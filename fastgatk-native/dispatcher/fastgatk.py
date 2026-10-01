#!/usr/bin/env python3
"""GATK-compatible dispatcher for the native fast-gatk tools.

This is intentionally a small, dependency-free launcher boundary.  It owns
launcher parsing and registry decisions, while a native tool owns tool
semantics.  Unknown tool parameters are explicit: by default this process
returns a usage error; with ``--fallback`` (or FASTGATK_ENABLE_FALLBACK=1) it
executes the configured GATK fallback command if it is available.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from urllib.parse import unquote, urlparse, urlsplit, urlunsplit
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence


EXIT_USAGE = 2
EXIT_BACKEND_UNAVAILABLE = 69
EXIT_FALLBACK_UNAVAILABLE = 70
EXIT_ARGS_FILE = 71
EXIT_REGISTRY = 72


class DispatcherError(RuntimeError):
    """An expected, user-facing dispatcher error."""

    def __init__(self, message: str, code: int = EXIT_USAGE,
                 category: str = "BAD_ARGUMENT") -> None:
        super().__init__(message)
        self.code = code
        self.category = category


class UnsupportedParameter(DispatcherError):
    def __init__(self, tool: str, option: str, argv: Sequence[str]) -> None:
        super().__init__(
            f"unsupported parameter for {tool}: {option}; "
            "no argument was discarded (use --fallback for explicit Java fallback)",
            EXIT_USAGE,
            "UNSUPPORTED_PARAMETER",
        )
        self.tool = tool
        self.option = option
        self.argv = list(argv)


@dataclass(frozen=True)
class RemoteOutput:
    uri: str
    local_path: Path


_REMOTE_INPUT_OPTIONS = {
    "-I", "--input", "-V", "--variant", "-R", "--reference",
    "-L", "--intervals", "--interval", "--region", "-XL", "--exclude-intervals",
    "--annotated-intervals", "--panel-of-normals", "--count-panel-of-normals",
    "--bqsr-recal-file", "--known-sites", "--resource", "--comparison",
    "--contamination-table", "--tumor-segmentation",
    "--orientation-bias-artifact-priors", "--ob-priors",
}
_OUTPUT_OPTIONS = {"-O", "--output", "--output-manifest", "--manifest",
                   "--output-vcf", "--output-bam", "--output-gvcf"}
_REMOTE_SCHEMES = {"http", "https", "s3", "gs", "gcs", "az", "abfs"}
_COMPANION_SUFFIXES = (".fai", ".dict", ".bai", ".crai", ".tbi", ".csi")


def _uri_scheme(value: str) -> str:
    parsed = urlparse(value)
    return parsed.scheme.lower() if parsed.scheme else ""


def _is_remote_uri(value: str) -> bool:
    return _uri_scheme(value) in _REMOTE_SCHEMES


def _with_suffix(uri: str, suffix: str) -> str:
    """Append a sidecar suffix to the URL path, before query/fragment data.

    Presigned URLs commonly carry authentication material in the query string;
    concatenating ``.fai`` to the complete URI would instead corrupt that
    query.  Keeping the query and fragment untouched also makes the helper
    safe for ordinary HTTP URLs.
    """
    parsed = urlsplit(uri)
    return urlunsplit((parsed.scheme, parsed.netloc, parsed.path + suffix,
                       parsed.query, parsed.fragment))


def _with_extension(uri: str, extension: str) -> str:
    """Replace the URL path extension while preserving auth query data."""
    parsed = urlsplit(uri)
    path = parsed.path
    dot = path.rfind(".")
    slash = path.rfind("/")
    if dot > slash:
        path = path[:dot] + extension
    else:
        path += extension
    return urlunsplit((parsed.scheme, parsed.netloc, path,
                       parsed.query, parsed.fragment))


def _file_uri_path(value: str) -> str | None:
    if _uri_scheme(value) != "file":
        return None
    parsed = urlparse(value)
    return unquote(parsed.path)


def _remote_cache_root() -> Path:
    configured = os.environ.get("FASTGATK_REMOTE_CACHE", "")
    if configured:
        root = Path(configured).expanduser()
    else:
        root = Path(tempfile.gettempdir()) / "fastgatk-remote-cache"
    root.mkdir(parents=True, exist_ok=True)
    return root


def _safe_remote_filename(uri: str, suffix: str = "") -> Path:
    parsed = urlparse(uri)
    name = Path(unquote(parsed.path)).name or "object"
    name = re.sub(r"[^A-Za-z0-9._-]", "_", name)[:160] or "object"
    digest = hashlib.sha256(uri.encode("utf-8")).hexdigest()[:24]
    return _remote_cache_root() / f"{digest}-{name}{suffix}"


def _download_remote(uri: str, destination: Path, required: bool = True) -> bool:
    if destination.is_file() and destination.stat().st_size > 0:
        return True
    temporary = destination.with_name(
        f".{destination.name}.part.{os.getpid()}"
    )
    command = [
        "curl", "--fail", "--silent", "--show-error", "--location",
        "--retry", "3", "--retry-delay", "1", "--connect-timeout", "30",
        "--output", str(temporary), uri,
    ]
    try:
        result = subprocess.run(command, check=False, text=True, capture_output=True)
    except OSError as exc:
        if required:
            raise DispatcherError(
                f"remote input requires curl: {exc}", EXIT_BACKEND_UNAVAILABLE,
                "BACKEND_UNAVAILABLE",
            ) from exc
        return False
    if result.returncode != 0 or not temporary.is_file() or temporary.stat().st_size == 0:
        try:
            temporary.unlink(missing_ok=True)
        except OSError:
            pass
        if required:
            detail = result.stderr.strip()[-1000:]
            raise DispatcherError(
                f"remote input download failed for {uri}{(': ' + detail) if detail else ''}",
                EXIT_BACKEND_UNAVAILABLE, "BACKEND_UNAVAILABLE",
            )
        return False
    os.replace(temporary, destination)
    return True


def _upload_remote(uri: str, source: Path) -> None:
    if _uri_scheme(uri) not in {"http", "https"}:
        raise DispatcherError(
            f"remote output upload currently requires HTTP(S)/presigned URI: {uri}",
            EXIT_BACKEND_UNAVAILABLE, "BACKEND_UNAVAILABLE",
        )
    command = [
        "curl", "--fail", "--silent", "--show-error", "--location",
        "--retry", "3", "--retry-delay", "1", "--connect-timeout", "30",
        "--upload-file", str(source), uri,
    ]
    try:
        result = subprocess.run(command, check=False, text=True, capture_output=True)
    except OSError as exc:
        raise DispatcherError(
            f"remote output requires curl: {exc}", EXIT_BACKEND_UNAVAILABLE,
            "BACKEND_UNAVAILABLE",
        ) from exc
    if result.returncode != 0:
        detail = result.stderr.strip()[-1000:]
        raise DispatcherError(
            f"remote output upload failed for {uri}{(': ' + detail) if detail else ''}",
            EXIT_BACKEND_UNAVAILABLE, "BACKEND_UNAVAILABLE",
        )


def _remote_output_local_path(uri: str) -> Path:
    destination = _safe_remote_filename(uri)
    output_root = _remote_cache_root() / "outputs"
    output_root.mkdir(parents=True, exist_ok=True)
    return output_root / destination.name


def _stage_remote_value(value: str, required_companions: bool = False) -> tuple[str, list[str]]:
    local_file = _file_uri_path(value)
    if local_file is not None:
        return local_file, []
    if not _is_remote_uri(value):
        return value, []
    destination = _safe_remote_filename(value)
    _download_remote(value, destination, required=True)
    staged = [value]
    # Local HTSlib/FAIDX readers discover sidecars by replacing the remote
    # suffix.  Cache companions under the same content-addressed stem.  A
    # reference requires .fai; BAM/VCF indexes are opportunistic because some
    # tools can scan unindexed inputs or use a caller-provided index.
    for companion in _COMPANION_SUFFIXES:
        # GATK dictionaries conventionally replace ``.fasta/.fa`` with
        # ``.dict`` while indexes append their suffix.  Try the canonical
        # dictionary URL first, then the appended spelling used by some
        # object stores; both are materialized at the local path expected by
        # HTSlib/GATK readers.
        if companion == ".dict":
            companion_uri = _with_extension(value, companion)
            companion_fallback_uri = _with_suffix(value, companion)
            companion_destination = destination.with_suffix(companion)
        else:
            companion_uri = _with_suffix(value, companion)
            companion_fallback_uri = ""
            companion_destination = destination.with_name(destination.name + companion)
        fetched = _download_remote(
            companion_uri, companion_destination,
            required=required_companions and companion == ".fai",
        )
        if not fetched and companion_fallback_uri:
            fetched = _download_remote(companion_fallback_uri, companion_destination,
                                       required=False)
            if fetched:
                companion_uri = companion_fallback_uri
        if fetched:
            staged.append(companion_uri)
    return str(destination), staged


def _remote_download_workers() -> int:
    """Return the bounded fan-out used for independent remote inputs.

    A dispatcher process is normally one task inside a SLURM/Nextflow
    allocation, so the default is deliberately small and does not consume a
    host-wide thread pool.  The environment override is useful for object
    stores with high per-request latency, but remains capped to keep a typo
    from creating an unbounded network storm.
    """
    raw = os.environ.get("FASTGATK_REMOTE_DOWNLOAD_THREADS", "").strip()
    if not raw:
        return min(4, max(1, os.cpu_count() or 1))
    try:
        value = int(raw)
    except ValueError as exc:
        raise DispatcherError(
            "FASTGATK_REMOTE_DOWNLOAD_THREADS must be an integer in [1, 64]",
            EXIT_USAGE, "BAD_ARGUMENT",
        ) from exc
    if value < 1 or value > 64:
        raise DispatcherError(
            "FASTGATK_REMOTE_DOWNLOAD_THREADS must be an integer in [1, 64]",
            EXIT_USAGE, "BAD_ARGUMENT",
        )
    return value


def _prefetch_remote_inputs(tool_args: Sequence[str]) -> tuple[
        dict[str, tuple[str, list[str]]], list[str]]:
    """Download independent input URIs concurrently, deterministically.

    The returned mapping is keyed by URI, so repeated ``-I/-V/-R`` values are
    fetched once.  Futures are consumed in first-seen order rather than
    completion order; this keeps manifest source ordering stable while still
    overlapping network latency.  Sidecar downloads for one input remain
    grouped inside ``_stage_remote_value`` and therefore retain the same
    required-reference semantics as the serial path.
    """
    requests: dict[str, bool] = {}
    index = 0
    while index < len(tool_args):
        token = tool_args[index]
        name = _option_name(token)
        resource_option = name.startswith("--resource:")
        if name not in _REMOTE_INPUT_OPTIONS and not resource_option:
            index += 1
            continue
        # VariantRecalibrator resource labels encode ``training=true`` and
        # similar attributes in the option token itself.  The following token
        # is always the resource URI; treating the first ``=`` as an inline
        # value would stage the label fragment instead of the file.
        inline = "=" in token and not resource_option
        value = token.split("=", 1)[1] if inline else (
            tool_args[index + 1] if index + 1 < len(tool_args) else ""
        )
        if value and _is_remote_uri(value) and _file_uri_path(value) is None:
            required = name in {"-R", "--reference"}
            requests[value] = requests.get(value, False) or required
        index += 1 if inline else 2
    if not requests:
        return {}, []

    ordered = list(requests.items())
    workers = min(_remote_download_workers(), len(ordered))
    if workers == 1:
        fetched = [
            _stage_remote_value(uri, required_companions=required)
            for uri, required in ordered
        ]
    else:
        with ThreadPoolExecutor(max_workers=workers,
                                thread_name_prefix="fastgatk-remote") as pool:
            futures = [pool.submit(_stage_remote_value, uri,
                                   required_companions=required)
                       for uri, required in ordered]
            # result() is intentionally read in input order.  If any request
            # fails, the exception is raised here and no native binary runs.
            fetched = [future.result() for future in futures]
    cache = {uri: result for (uri, _), result in zip(ordered, fetched)}
    sources = [source for _, result in zip(ordered, fetched)
               for source in result[1]]
    return cache, sources


def _remote_output_mode() -> str:
    return os.environ.get("FASTGATK_REMOTE_OUTPUT_MODE", "fail").strip().lower()


def _stage_remote_io(tool_args: Sequence[str], allow_remote_outputs: bool = False) -> tuple[
        list[str], list[str], list[RemoteOutput]]:
    """Stage remote input URI values without invoking a shell.

    The default is an input-only boundary.  When explicitly enabled by
    ``FASTGATK_REMOTE_OUTPUT_MODE=upload``, output URIs are rewritten to an
    atomic local staging path and returned for post-process upload.  The
    original values remain untouched for explicit Java fallback.
    """
    staged_args = list(tool_args)
    staged_cache, staged_sources = _prefetch_remote_inputs(tool_args)
    remote_outputs: list[RemoteOutput] = []
    index = 0
    while index < len(staged_args):
        token = staged_args[index]
        name = _option_name(token)
        if name in _OUTPUT_OPTIONS:
            value = token.split("=", 1)[1] if "=" in token else (
                staged_args[index + 1] if index + 1 < len(staged_args) else ""
            )
            local_output = _file_uri_path(value)
            if _is_remote_uri(value):
                if not allow_remote_outputs:
                    raise DispatcherError(
                        f"native remote output is not enabled: {value}; set "
                        "FASTGATK_REMOTE_OUTPUT_MODE=upload or use --fallback",
                        EXIT_BACKEND_UNAVAILABLE, "BACKEND_UNAVAILABLE",
                    )
                local_path = _remote_output_local_path(value)
                if "=" in token:
                    staged_args[index] = f"{name}={local_path}"
                elif index + 1 < len(staged_args):
                    staged_args[index + 1] = str(local_path)
                remote_outputs.append(RemoteOutput(value, local_path))
                index += 1 if "=" in token else 2
                continue
            if local_output is not None:
                if "=" in token:
                    staged_args[index] = f"{name}={local_output}"
                elif index + 1 < len(staged_args):
                    staged_args[index + 1] = local_output
            index += 1 if "=" in token else 2
            continue
        # VariantRecalibrator encodes a resource label in the option name
        # (``--resource:hapmap,training=true,...``), so it cannot be listed as
        # one finite registry option.  Its following token is still a normal
        # remote input and must pass through the same staging boundary.
        resource_option = name.startswith("--resource:")
        if name not in _REMOTE_INPUT_OPTIONS and not resource_option:
            index += 1
            continue
        # See the matching note in _prefetch_remote_inputs: resource labels
        # contain '=' but their value is the following token.
        inline = "=" in token and not resource_option
        value = token.split("=", 1)[1] if inline else (
            staged_args[index + 1] if index + 1 < len(staged_args) else ""
        )
        if value:
            staged_value, _sources = staged_cache.get(value, (None, None))
            if staged_value is None:
                staged_value, _sources = _stage_remote_value(
                    value, required_companions=name in {"-R", "--reference"}
                )
            if inline:
                staged_args[index] = f"{name}={staged_value}"
            elif staged_value != value:
                staged_args[index + 1] = staged_value
        index += 1 if inline else 2
    return staged_args, staged_sources, remote_outputs


def stage_remote_inputs(tool_args: Sequence[str]) -> tuple[list[str], list[str]]:
    """Backward-compatible input-only staging helper used by contracts."""
    staged, sources, _ = _stage_remote_io(tool_args, allow_remote_outputs=False)
    return staged, sources


def _rewrite_remote_manifests(outputs: Sequence[RemoteOutput]) -> None:
    replacements = [(str(output.local_path), output.uri) for output in outputs]
    for output in outputs:
        if not output.local_path.name.endswith(".manifest.json") or not output.local_path.is_file():
            continue
        text = output.local_path.read_text(encoding="utf-8")
        for local, remote in replacements:
            text = text.replace(local, remote)
        temporary = output.local_path.with_name(f".{output.local_path.name}.rewrite.{os.getpid()}")
        temporary.write_text(text, encoding="utf-8")
        os.replace(temporary, output.local_path)


def commit_remote_outputs(outputs: Sequence[RemoteOutput]) -> None:
    """Upload sidecars before primary objects and fail closed on any error."""
    if not outputs:
        return
    _rewrite_remote_manifests(outputs)
    sidecar_suffixes = (".fai", ".bai", ".crai", ".tbi", ".csi", ".md5")
    for output in outputs:
        if not output.local_path.is_file() or output.local_path.stat().st_size == 0:
            raise DispatcherError(
                f"native output was not produced: {output.local_path}",
                EXIT_BACKEND_UNAVAILABLE, "BACKEND_UNAVAILABLE",
            )
        for suffix in sidecar_suffixes:
            sidecar = output.local_path.with_name(output.local_path.name + suffix)
            if sidecar.is_file() and sidecar.stat().st_size > 0:
                _upload_remote(_with_suffix(output.uri, suffix), sidecar)
        # Commit the primary object last.  Object stores that expose the
        # presigned URL therefore never observe a complete primary before its
        # indexes have been uploaded.
        _upload_remote(output.uri, output.local_path)


@dataclass(frozen=True)
class Registry:
    path: Path
    payload: dict[str, Any]

    @property
    def version(self) -> str:
        return str(self.payload.get("registry_version", "unknown"))

    @property
    def schema_version(self) -> int:
        return int(self.payload.get("schema_version", 0))

    @property
    def tools(self) -> dict[str, dict[str, Any]]:
        tools = self.payload.get("tools")
        if not isinstance(tools, dict) or not tools:
            raise DispatcherError(
                f"registry has no tools: {self.path}", EXIT_REGISTRY, "REGISTRY_INVALID"
            )
        return tools

    def resolve_tool(self, name: str) -> tuple[str, dict[str, Any]]:
        lowered = name.lower()
        for canonical, entry in self.tools.items():
            aliases = entry.get("aliases", [])
            if canonical.lower() == lowered or any(str(alias).lower() == lowered for alias in aliases):
                if not isinstance(entry, dict):
                    raise DispatcherError(
                        f"registry entry is not an object: {canonical}",
                        EXIT_REGISTRY,
                        "REGISTRY_INVALID",
                    )
                return canonical, entry
        raise DispatcherError(
            f"unsupported GATK tool: {name}; use --list or explicit Java fallback",
            EXIT_USAGE,
            "UNSUPPORTED_TOOL",
        )


@dataclass(frozen=True)
class Parsed:
    launcher_args: tuple[str, ...]
    tool: str
    tool_args: tuple[str, ...]
    dry_run: bool
    fallback: bool
    launcher_separator: bool = False


def _default_registry_path() -> Path:
    return Path(__file__).resolve().with_name("tool_registry.json")


def load_registry(path: Path | None = None) -> Registry:
    registry_path = path or Path(os.environ.get("FASTGATK_TOOL_REGISTRY", _default_registry_path()))
    try:
        payload = json.loads(registry_path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise DispatcherError(
            f"tool registry not found: {registry_path}", EXIT_REGISTRY, "REGISTRY_UNAVAILABLE"
        ) from exc
    except json.JSONDecodeError as exc:
        raise DispatcherError(
            f"invalid tool registry {registry_path}: {exc}", EXIT_REGISTRY, "REGISTRY_INVALID"
        ) from exc
    if not isinstance(payload, dict):
        raise DispatcherError(
            f"invalid tool registry root: {registry_path}", EXIT_REGISTRY, "REGISTRY_INVALID"
        )
    return Registry(registry_path, payload)


def _split_args_file(text: str, source: Path) -> list[str]:
    try:
        # GATK argument files are whitespace-separated and support quoting.  We
        # intentionally do not enable shell comments: a '#' in a path/value is
        # data unless the caller explicitly quotes a different representation.
        return shlex.split(text, comments=False, posix=True)
    except ValueError as exc:
        raise DispatcherError(
            f"cannot parse argument file {source}: {exc}", EXIT_ARGS_FILE, "ARGS_FILE_INVALID"
        ) from exc


def expand_argument_files(argv: Sequence[str], cwd: Path | None = None,
                          max_depth: int = 16) -> list[str]:
    """Expand ``@file`` tokens without invoking a shell.

    Nested paths are resolved relative to the process working directory, which
    matches how GATK jobs generally stage argument files.  ``@@x`` escapes a
    literal ``@x`` token.  The expansion is deterministic and bounded.
    """

    root = cwd or Path.cwd()

    def expand(items: Iterable[str], depth: int) -> list[str]:
        if depth > max_depth:
            raise DispatcherError(
                f"argument file nesting exceeds {max_depth}", EXIT_ARGS_FILE, "ARGS_FILE_DEPTH"
            )
        result: list[str] = []
        for item in items:
            if item == "@@":
                result.append("@")
                continue
            if item.startswith("@@"):
                result.append(item[1:])
                continue
            if not item.startswith("@") or item == "@":
                result.append(item)
                continue
            path = Path(item[1:])
            if not path.is_absolute():
                path = root / path
            try:
                nested_text = path.read_text(encoding="utf-8")
            except OSError as exc:
                raise DispatcherError(
                    f"cannot read argument file {path}: {exc}", EXIT_ARGS_FILE, "ARGS_FILE_UNAVAILABLE"
                ) from exc
            result.extend(expand(_split_args_file(nested_text, path), depth + 1))
        return result

    return expand(list(argv), 0)


def _option_name(token: str) -> str:
    return token.split("=", 1)[0]


_OPTION_ALIASES: dict[str, set[str]] = {
    "-I": {"-I", "--input"},
    "-V": {"-V", "--variant"},
    "-R": {"-R", "--reference"},
    "-L": {"-L", "--intervals", "--interval", "--region"},
    "-O": {"-O", "--output"},
}

# Native read-filter masks have a finite, auditable class surface.  Most tools
# share this conservative HC/Mutect2 mask; read-QC tools can advertise a
# larger per-entry mask when their native implementation handles additional
# flag/coordinate predicates.  Dispatcher validation mirrors the selected
# surface so an arbitrary GATK filter class is never accepted by the native
# path and silently ignored.  With --fallback the full original argv is still
# forwarded to Java.
_NATIVE_READ_FILTERS = {
    "MappingQualityReadFilter",
    "NotDuplicateReadFilter",
    "MappedReadFilter",
    "NotSecondaryAlignmentReadFilter",
    "NotSupplementaryAlignmentReadFilter",
    "PassesVendorQualityCheckReadFilter",
    "GoodCigarReadFilter",
    "WellformedReadFilter",
    "NonZeroReferenceLengthAlignmentReadFilter",
    "ReadLengthReadFilter",
}


def _is_option(token: str, names: set[str]) -> bool:
    return _option_name(token) in names


def _consume_launcher(argv: Sequence[str]) -> Parsed:
    """Separate launcher options/tool from tool arguments.

    The parser accepts ``fastgatk HaplotypeCaller ...`` and, for compatibility
    with existing process wrappers, an implicit HC invocation beginning with a
    tool option (``fastgatk -I ...``).  Unknown launcher options are errors,
    never passed accidentally to a native binary.
    """

    launcher: list[str] = []
    tokens = list(argv)
    dry_run = False
    fallback = os.environ.get("FASTGATK_ENABLE_FALLBACK", "0") == "1"
    launcher_separator = False
    index = 0
    tool_index: int | None = None
    while index < len(tokens):
        token = tokens[index]
        if token == "--":
            if index + 1 >= len(tokens):
                raise DispatcherError("-- must be followed by a tool or tool arguments")
            launcher_separator = True
            tool_index = index + 1
            break
        if token in {"--dry-run", "--dryrun"}:
            dry_run = True
            launcher.append(token)
            index += 1
            continue
        if token in {"--fallback", "--java-fallback"}:
            fallback = True
            launcher.append(token)
            index += 1
            continue
        if token in {"--help", "-h", "--version", "--list"}:
            # Global actions are handled by main.  Keep them here only so a
            # tool name may appear before --help without becoming an unknown
            # launcher option.
            launcher.append(token)
            index += 1
            continue
        if token in {"--java-options", "-java-options", "--gatk-config-file"}:
            if index + 1 >= len(tokens):
                raise DispatcherError(f"missing value for launcher option {token}")
            launcher.extend([token, tokens[index + 1]])
            index += 2
            continue
        if (token.startswith("--java-options=") or token.startswith("-java-options=") or
                token.startswith("--gatk-config-file=")):
            launcher.append(token)
            index += 1
            continue
        if token.startswith("-"):
            # If the first non-launcher token is a tool option, allow the
            # default HaplotypeCaller process shape used by old wrappers.
            tool_index = index
            break
        tool_index = index
        break
    if tool_index is None:
        # Only global actions are allowed without a tool.  main handles them;
        # otherwise report a useful usage error.
        if launcher:
            return Parsed(tuple(launcher), "", tuple(), dry_run, fallback,
                          launcher_separator)
        raise DispatcherError("a GATK tool is required (for example HaplotypeCaller)")

    first = tokens[tool_index]
    if first.startswith("-"):
        tool = "HaplotypeCaller"
        raw_tool_args = list(tokens[tool_index:])
    else:
        tool = first
        raw_tool_args = list(tokens[tool_index + 1:])
    # Permit dispatcher controls after the tool as well.  They are removed
    # before registry validation, so they cannot accidentally reach a native
    # tool or be mistaken for a GATK option.  GATK's launcher options are also
    # accepted after the tool (the reference launcher searches the complete
    # argv), but are normalized into launcher_args so native binaries never
    # receive JVM/configuration controls they cannot honor.
    tool_args: list[str] = []
    index = 0
    after_tool_separator = False
    while index < len(raw_tool_args):
        token = raw_tool_args[index]
        if after_tool_separator:
            tool_args.append(token)
        elif token == "--":
            after_tool_separator = True
            tool_args.append(token)
        elif token in {"--dry-run", "--dryrun"}:
            dry_run = True
        elif token in {"--fallback", "--java-fallback"}:
            fallback = True
        elif token in {"--java-options", "-java-options", "--gatk-config-file"}:
            if index + 1 >= len(raw_tool_args):
                raise DispatcherError(f"missing value for launcher option {token}")
            launcher.extend([token, raw_tool_args[index + 1]])
            index += 1
        elif (token.startswith("--java-options=") or token.startswith("-java-options=") or
              token.startswith("--gatk-config-file=")):
            launcher.append(token)
        else:
            # A second -- belongs to GATK's tool/Spark argument split.  Keep
            # it and everything after it byte-for-byte for explicit fallback;
            # native validation rejects this boundary unless a fallback was
            # requested, so unsupported Spark tails cannot be silently lost.
            tool_args.append(token)
        index += 1
    return Parsed(tuple(launcher), tool, tuple(tool_args), dry_run, fallback,
                  launcher_separator)


def _split_option_value(tool: str, token: str, tokens: Sequence[str], index: int,
                        value_options: set[str], flag_options: set[str],
                        optional_boolean_options: set[str]) -> tuple[str, int]:
    name = _option_name(token)
    if name in flag_options:
        if "=" in token:
            if name not in optional_boolean_options:
                raise DispatcherError(f"flag option cannot have a value: {token}")
            value = token.split("=", 1)[1].lower()
            if value not in {"true", "false", "1", "0"}:
                raise DispatcherError(f"boolean option expects true or false: {token}")
            return token, index + 1
        if (name in optional_boolean_options and index + 1 < len(tokens) and
                tokens[index + 1].lower() in {"true", "false", "1", "0"}):
            return token, index + 2
        return token, index + 1
    if name not in value_options:
        raise UnsupportedParameter(tool, token, tokens)
    if "=" in token:
        value = token.split("=", 1)[1]
        if not value:
            raise DispatcherError(f"empty value for option {name}")
        return token, index + 1
    if index + 1 >= len(tokens):
        raise DispatcherError(f"missing value for option {token}")
    return token, index + 2


def validate_tool_args(tool: str, entry: dict[str, Any], tool_args: Sequence[str]) -> None:
    value_options = {str(x) for x in entry.get("value_options", [])}
    flag_options = {str(x) for x in entry.get("flag_options", [])}
    optional_boolean_options = {str(x) for x in entry.get("optional_boolean_options", [])}
    native_filter_table = entry.get("native_read_filters", _NATIVE_READ_FILTERS)
    if (not isinstance(entry.get("value_options", []), list) or
            not isinstance(entry.get("flag_options", []), list) or
            not isinstance(entry.get("optional_boolean_options", []), list) or
            not optional_boolean_options.issubset(flag_options) or
            not isinstance(native_filter_table, (list, set, tuple))):
        raise DispatcherError(f"invalid option table for {tool}", EXIT_REGISTRY, "REGISTRY_INVALID")
    index = 0
    seen: dict[str, int] = {}
    while index < len(tool_args):
        token = tool_args[index]
        if token == "--":
            raise DispatcherError(f"unexpected -- after tool {tool}")
        if not token.startswith("-"):
            raise DispatcherError(f"unexpected positional argument for {tool}: {token}")
        name = _option_name(token)
        if tool == "VariantRecalibrator" and name.startswith("--resource:"):
            # GATK resource labels are encoded in the option name, for example
            # --resource:hapmap,training=true,truth=true path.vcf.gz.  Keep the
            # complete labelled token and its following path intact instead of
            # forcing every release-specific label into the registry.
            seen[name] = seen.get(name, 0) + 1
            # Resource labels themselves contain '=' (training=true, etc.),
            # so an '=' never denotes an inline path in this GATK syntax.
            if index + 1 >= len(tool_args):
                raise DispatcherError(f"missing value for option {token}")
            index += 2
            continue
        if name not in value_options and name not in flag_options:
            raise UnsupportedParameter(tool, token, tool_args)
        seen[name] = seen.get(name, 0) + 1
        if name in {"-RF", "--read-filter", "-DF", "--disable-read-filter"}:
            filter_name = token.split("=", 1)[1] if "=" in token else (
                tool_args[index + 1] if index + 1 < len(tool_args) else ""
            )
            if not filter_name or filter_name.startswith("-"):
                raise DispatcherError(f"missing value for option {token}")
            native_filters = {str(value) for value in native_filter_table}
            if filter_name not in native_filters:
                raise UnsupportedParameter(tool, f"{name} {filter_name}", tool_args)
        _, index = _split_option_value(tool, token, tool_args, index, value_options,
                                       flag_options, optional_boolean_options)
    required_inputs = {str(x) for x in entry.get("required_inputs", [])}
    required_outputs = {str(x) for x in entry.get("required_outputs", [])}
    for required in required_inputs | required_outputs:
        # A registry requirement may name mutually exclusive alternatives as
        # ``-I|-V`` (for example, CheckReferenceCompatibility accepts either
        # an alignment input or a variant input).  Treat the whole expression
        # as one requirement while still honoring aliases for each option.
        alternatives = [candidate.strip() for candidate in required.split("|")
                        if candidate.strip()]
        aliases: set[str] = set()
        for candidate in alternatives or [required]:
            aliases.update(_OPTION_ALIASES.get(candidate, {candidate}))
        if not any(_option_name(token) in aliases for token in tool_args):
            raise DispatcherError(f"missing required option for {tool}: {required}")
    for group in entry.get("required_output_groups", []):
        aliases: set[str] = set()
        for required in group:
            aliases.update(_OPTION_ALIASES.get(str(required), {str(required)}))
        if not any(_option_name(token) in aliases for token in tool_args):
            rendered = ", ".join(str(value) for value in group)
            raise DispatcherError(f"one of the required output options is missing for {tool}: {rendered}")
    repeatable = {str(x) for x in entry.get("repeatable_options", [])}
    if (sum(seen.get(alias, 0) for alias in _OPTION_ALIASES["-I"]) > 1 and
            not any(alias in repeatable for alias in _OPTION_ALIASES["-I"])):
        raise DispatcherError(f"duplicate input option for {tool}")
    if (sum(seen.get(alias, 0) for alias in _OPTION_ALIASES["-O"]) > 1 and
            not any(alias in repeatable for alias in _OPTION_ALIASES["-O"])):
        raise DispatcherError(f"duplicate output option for {tool}")


def _project_root() -> Path:
    # dispatcher/ -> fastgatk-native/ -> workspace root
    return Path(__file__).resolve().parents[2]


def resolve_native_binary(canonical: str, entry: dict[str, Any]) -> Path:
    env_name = "FASTGATK_" + "".join(character if character.isalnum() else "_"
                                      for character in canonical.upper()) + "_BINARY"
    configured = os.environ.get(env_name)
    if canonical == "HaplotypeCaller" and not configured:
        # Preserve the original override name used by existing wrappers.
        configured = os.environ.get("FASTGATK_HC_BINARY")
    value = configured or str(entry.get("native_binary", ""))
    if not value:
        raise DispatcherError("registry has no native binary", EXIT_REGISTRY, "REGISTRY_INVALID")
    path = Path(value)
    if not path.is_absolute():
        configured_build = os.environ.get("FASTGATK_NATIVE_BUILD", "").strip()
        if configured_build:
            # Registry paths are rooted at fastgatk-native/build.  Preserve
            # the suffix below that build directory when a deployment uses a
            # custom artifact directory (for example a SLURM module cache).
            parts = path.parts
            try:
                build_index = parts.index("build")
            except ValueError:
                build_index = len(parts) - 1
            suffix = Path(*parts[build_index + 1:]) if build_index + 1 < len(parts) else Path(path.name)
            path = Path(configured_build).expanduser() / suffix
        else:
            path = _project_root() / path
            # A serial portability build is a valid native artifact too.  If
            # the conventional OpenMP directory is absent, probe the sibling
            # build-serial directory without changing registry semantics.
            if not path.is_file() and "/build/" in path.as_posix():
                serial = Path(path.as_posix().replace("/build/", "/build-serial/", 1))
                if serial.is_file():
                    path = serial
    return path


_LAUNCHER_VALUE_OPTIONS = {"--java-options", "-java-options", "--gatk-config-file"}


def _forwarded_launcher_args(launcher_args: Sequence[str]) -> list[str]:
    """Return only launcher options understood by the real GATK wrapper.

    Dispatcher-only controls such as ``--dry-run`` and ``--fallback`` must
    never leak into a Java invocation.  JVM/configuration options, however,
    are part of GATK's launcher contract and must be placed before the tool
    name even when the caller supplied them after the tool.
    """
    forwarded: list[str] = []
    index = 0
    while index < len(launcher_args):
        token = launcher_args[index]
        name = _option_name(token)
        if name in _LAUNCHER_VALUE_OPTIONS:
            forwarded.append(token)
            if "=" not in token:
                if index + 1 >= len(launcher_args):
                    raise DispatcherError(
                        f"missing value for launcher option {token}",
                        EXIT_USAGE, "BAD_ARGUMENT",
                    )
                forwarded.append(launcher_args[index + 1])
                index += 1
        index += 1
    return forwarded


def _launcher_option_values(launcher_args: Sequence[str], names: set[str]) -> list[str]:
    """Extract values from normalized launcher options without shell parsing."""
    values: list[str] = []
    index = 0
    while index < len(launcher_args):
        token = launcher_args[index]
        name = _option_name(token)
        if name in names:
            if "=" in token:
                value = token.split("=", 1)[1]
            elif index + 1 < len(launcher_args):
                value = launcher_args[index + 1]
                index += 1
            else:
                raise DispatcherError(
                    f"missing value for launcher option {token}",
                    EXIT_USAGE, "BAD_ARGUMENT",
                )
            if not value:
                raise DispatcherError(
                    f"empty value for launcher option {name}",
                    EXIT_USAGE, "BAD_ARGUMENT",
                )
            values.append(value)
        index += 1
    return values


def _fallback_argv(canonical: str, entry: dict[str, Any], tool_args: Sequence[str],
                   launcher_args: Sequence[str] = ()) -> list[str]:
    configured = os.environ.get("FASTGATK_GATK_BINARY")
    # FASTGATK_GATK_BINARY denotes the launcher, not a pre-bound tool
    # command.  Keep the canonical tool in the argv when an override is used
    # so a generic/unregistered fallback behaves exactly like a registry
    # fallback (`gatk Tool ...`).  Registry entries without an override keep
    # their explicit command verbatim for backwards compatibility.
    if configured:
        command = [configured, *_forwarded_launcher_args(launcher_args), canonical]
    else:
        configured_command = [str(x) for x in entry.get("fallback_command", [])]
        if not configured_command:
            raise DispatcherError("registry has no fallback command", EXIT_REGISTRY, "REGISTRY_INVALID")
        command = [configured_command[0], *_forwarded_launcher_args(launcher_args),
                   *configured_command[1:]]
    if not command:
        raise DispatcherError("registry has no fallback command", EXIT_REGISTRY, "REGISTRY_INVALID")
    return command + list(tool_args)


def _run_fallback(registry: Registry, canonical: str, entry: dict[str, Any],
                  tool_args: Sequence[str], reason: str, dry_run: bool,
                  launcher_args: Sequence[str] = ()) -> int:
    fallback = _fallback_argv(canonical, entry, tool_args, launcher_args)
    if dry_run:
        _print_json({
            "status": "dry-run",
            "execution_mode": "fallback",
            "fallback_reason": reason,
            "tool": canonical,
            "argv": fallback,
            "registry_version": registry.version,
        })
        return 0
    executable = shutil.which(fallback[0]) if not Path(fallback[0]).is_absolute() else fallback[0]
    if executable is None:
        raise DispatcherError(
            f"fallback executable not found: {fallback[0]}",
            EXIT_FALLBACK_UNAVAILABLE,
            "FALLBACK_UNAVAILABLE",
        )
    print(f"fastgatk: explicit fallback ({reason}) -> {' '.join(fallback)}", file=sys.stderr)
    return subprocess.run([executable] + fallback[1:], check=False).returncode


def _print_json(payload: dict[str, Any], stream: Any = sys.stdout) -> None:
    print(json.dumps(payload, sort_keys=True, ensure_ascii=False), file=stream)


def print_usage(registry: Registry | None = None) -> None:
    print("Usage: fastgatk [launcher-options] TOOL [tool-options]")
    print("       fastgatk [launcher-options] -I reads.bam -O calls.vcf")
    print("Launcher options: --help --version --list --dry-run --fallback "
          "--java-options VALUE --gatk-config-file PATH")
    print("Supported tools:")
    if registry is None:
        print("  HaplotypeCaller")
        return
    for canonical, entry in registry.tools.items():
        print(f"  {canonical}: status={entry.get('status', 'unknown')} aliases={','.join(entry.get('aliases', []))}")


def print_tool_help(tool: str, entry: dict[str, Any]) -> None:
    print(f"fastgatk {tool} ({entry.get('status', 'unknown')})")
    print("  -I, --input FILE                 input SAM/BAM/CRAM")
    print("  -R, --reference FILE             reference FASTA")
    print("  -L, --intervals REGION           contig:start-end")
    print("  -O, --output FILE                VCF or JSON summary")
    print("      --dry-run                    validate and print dispatch plan")
    print("      --fallback                   execute configured Java GATK fallback")
    print("  Registry-supported options:")
    print("    " + " ".join(entry.get("value_options", [])))
    print("    " + " ".join(entry.get("flag_options", [])))
    for note in entry.get("notes", []):
        print(f"  note: {note}")


def dispatch(registry: Registry, parsed: Parsed) -> int:
    if not parsed.tool:
        return 0
    try:
        canonical, entry = registry.resolve_tool(parsed.tool)
    except DispatcherError as error:
        # The execution plan deliberately keeps the long tail fallback-first:
        # a workflow using an unregistered GATK walker must still be able to
        # run unchanged when the caller explicitly opts into Java fallback.
        # Without --fallback we retain the fail-closed unknown-tool error so
        # a typo can never silently launch Java.
        if error.category != "UNSUPPORTED_TOOL" or not parsed.fallback:
            raise
        canonical = parsed.tool
        entry = {
            "status": "fallback-only",
            "aliases": [],
            "fallback_command": ["gatk", canonical],
            "fallback_policy": "explicit",
            "determinism": ["java"],
            "backends": ["GATK"],
            "value_options": [],
            "flag_options": [],
            "notes": ["Unregistered long-tail GATK tool; explicit Java fallback preserves the original argv."],
        }
    if any(token in {"--help", "-h"} for token in parsed.launcher_args) or any(
        token in {"--help", "-h"} for token in parsed.tool_args
    ):
        print_tool_help(canonical, entry)
        return 0

    # A fallback-only entry is intentionally explicit.  It preserves the
    # complete original argv for GATK instead of maintaining a partial native
    # allow-list that could reject a valid GATK release-specific parameter.
    if entry.get("status") == "fallback-only":
        if not parsed.fallback:
            raise DispatcherError(
                f"{canonical} is fallback-only; use --fallback or "
                "FASTGATK_ENABLE_FALLBACK=1",
                EXIT_BACKEND_UNAVAILABLE,
                "FALLBACK_REQUIRED",
            )
        return _run_fallback(
            registry, canonical, entry, parsed.tool_args,
            "registry status fallback-only", parsed.dry_run,
            parsed.launcher_args,
        )

    # A GATK properties file can alter codec, cloud, Spark, and tool defaults
    # that are not represented in the native registry.  Passing it through to
    # a native binary or adapter would silently ignore those settings, so
    # preserve the direct-replacement contract by routing this invocation to
    # Java.  The path/value itself remains in the fallback argv and is never
    # shell interpolated.
    config_files = _launcher_option_values(
        parsed.launcher_args, {"--gatk-config-file"}
    )
    if config_files:
        return _run_fallback(
            registry, canonical, entry, parsed.tool_args,
            "--gatk-config-file requires the Java GATK configuration boundary",
            parsed.dry_run, parsed.launcher_args,
        )

    # Adapter entries own the complete argv boundary and may forward
    # release-specific GATK parameters to an external backend.  Do not apply
    # the native allow-list here: the adapter must never discard an argument
    # merely because this registry version does not know it yet.
    if entry.get("status") == "adapter":
        native = resolve_native_binary(canonical, entry)
        plan = {
            "status": "dry-run" if parsed.dry_run else "dispatch",
            "execution_mode": "adapter",
            "tool": canonical,
            "binary": str(native),
            "binary_exists": native.is_file() and os.access(native, os.X_OK),
            "argv": [str(native), *parsed.tool_args],
            "launcher_args": list(parsed.launcher_args),
            "launcher_separator": parsed.launcher_separator,
            "registry": str(registry.path),
            "registry_version": registry.version,
            "determinism": os.environ.get("FASTGATK_DETERMINISM", "strict"),
        }
        if parsed.dry_run:
            if not plan["binary_exists"] and parsed.fallback:
                return _run_fallback(
                    registry, canonical, entry, parsed.tool_args,
                    f"adapter binary unavailable: {native}", True,
                    parsed.launcher_args,
                )
            _print_json(plan)
            return 0
        if not plan["binary_exists"]:
            if parsed.fallback:
                return _run_fallback(
                    registry, canonical, entry, parsed.tool_args,
                    f"adapter binary unavailable: {native}", False,
                    parsed.launcher_args,
                )
            raise DispatcherError(
                f"adapter binary unavailable: {native}; use --fallback for explicit Java fallback",
                EXIT_BACKEND_UNAVAILABLE,
                "BACKEND_UNAVAILABLE",
            )
        native_args, _, remote_outputs = _stage_remote_io(
            parsed.tool_args, allow_remote_outputs=_remote_output_mode() == "upload")
        return_code = subprocess.run([str(native), *native_args], check=False).returncode
        if return_code == 0:
            commit_remote_outputs(remote_outputs)
        return return_code

    # The post-tool separator is the GATK launcher boundary used for Spark
    # arguments.  Native non-Spark tools cannot consume that second argument
    # group; explicit fallback preserves it, while the default remains
    # fail-closed instead of dropping unknown arguments.
    if "--" in parsed.tool_args:
        if parsed.fallback:
            return _run_fallback(
                registry, canonical, entry, parsed.tool_args,
                "post-tool -- argument separator requires the Java launcher",
                parsed.dry_run, parsed.launcher_args,
            )
        raise DispatcherError(
            f"post-tool -- argument separator is unsupported for native {canonical}; "
            "use --fallback for the Java launcher",
            EXIT_USAGE, "UNSUPPORTED_PARAMETER",
        )

    try:
        validate_tool_args(canonical, entry, parsed.tool_args)
    except UnsupportedParameter as exc:
        if not parsed.fallback:
            raise
        return _run_fallback(registry, canonical, entry, parsed.tool_args, str(exc),
                             parsed.dry_run, parsed.launcher_args)

    native = resolve_native_binary(canonical, entry)
    plan = {
        "status": "dry-run" if parsed.dry_run else "dispatch",
        "execution_mode": "native",
        "tool": canonical,
        "binary": str(native),
        "binary_exists": native.is_file() and os.access(native, os.X_OK),
        "argv": [str(native), *parsed.tool_args],
        "launcher_args": list(parsed.launcher_args),
        "launcher_separator": parsed.launcher_separator,
        "registry": str(registry.path),
        "registry_version": registry.version,
        "determinism": os.environ.get("FASTGATK_DETERMINISM", "strict"),
    }
    if parsed.dry_run:
        if not plan["binary_exists"] and parsed.fallback:
            return _run_fallback(
                registry, canonical, entry, parsed.tool_args,
                f"native binary unavailable: {native}", True,
                parsed.launcher_args,
            )
        _print_json(plan)
        return 0
    if not plan["binary_exists"]:
        if parsed.fallback:
            return _run_fallback(
                registry, canonical, entry, parsed.tool_args,
                f"native binary unavailable: {native}", False,
                parsed.launcher_args,
            )
        raise DispatcherError(
            f"native binary unavailable: {native}; use --fallback for explicit Java fallback",
            EXIT_BACKEND_UNAVAILABLE,
            "BACKEND_UNAVAILABLE",
        )
    native_args, _, remote_outputs = _stage_remote_io(
        parsed.tool_args, allow_remote_outputs=_remote_output_mode() == "upload")
    return_code = subprocess.run([str(native), *native_args], check=False).returncode
    if return_code == 0:
        commit_remote_outputs(remote_outputs)
    return return_code


def main(argv: Sequence[str] | None = None) -> int:
    raw = list(sys.argv[1:] if argv is None else argv)
    # OpenMP spin budget: libgomp's default spin count burns CPU at
    # barrier-heavy kernels (DREAM somatic: 33 CPU-s vs 19 at spin5000, for
    # ~10% wall).  libgomp parses GOMP_SPINCOUNT at startup, so the launcher
    # must set it before exec; an explicit user setting always wins.  See
    # fastgatk-native/docs/resource-parity.md for the measured frontier.
    os.environ.setdefault("GOMP_SPINCOUNT", "5000")
    try:
        expanded = expand_argument_files(raw)
        registry = load_registry()
        if not expanded or expanded in (["--help"], ["-h"]):
            print_usage(registry)
            return 0
        if expanded in (["--version"],):
            print(f"fastgatk {registry.version} (registry schema {registry.schema_version})")
            return 0
        if expanded in (["--list"],):
            for canonical, entry in registry.tools.items():
                print(f"{canonical}\t{entry.get('status', 'unknown')}\t{','.join(entry.get('aliases', []))}")
            return 0
        parsed = _consume_launcher(expanded)
        if not parsed.tool and any(token in {"--help", "-h"}
                                   for token in parsed.launcher_args):
            print_usage(registry)
            return 0
        if not parsed.tool and any(token == "--version"
                                   for token in parsed.launcher_args):
            print(f"fastgatk {registry.version} (registry schema {registry.schema_version})")
            return 0
        if parsed.launcher_args and "--list" in parsed.launcher_args and not parsed.tool:
            for canonical, entry in registry.tools.items():
                print(f"{canonical}\t{entry.get('status', 'unknown')}\t{','.join(entry.get('aliases', []))}")
            return 0
        return dispatch(registry, parsed)
    except DispatcherError as exc:
        _print_json({
            "status": "error",
            "category": exc.category,
            "message": str(exc),
            "tool": getattr(exc, "tool", None),
            "registry": str(_default_registry_path()),
        }, sys.stderr)
        return exc.code
    except KeyboardInterrupt:
        _print_json({"status": "error", "category": "INTERRUPTED", "message": "interrupted"}, sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
