#!/usr/bin/env python3
"""Locked HG002 acquisition; no pipeline execution or biological eligibility claims.

fetch_assets defaults to small files. Callers must enforce the resource gate before
opting into large=True. Every reuse hashes the bytes again. Unpinned publisher
manifests have an explicit HTTPS trust-on-first-use boundary and immutable local
receipts; dependent files additionally match the checked-in MD5 pins.
"""
from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import http.client
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import tarfile
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import uuid

LOCK_PATH = Path(__file__).with_suffix(".json")
CHUNK_BYTES = 1024 * 1024
MAX_METADATA_BYTES = 8 * 1024 * 1024
MAX_LINE_BYTES = 65536
MAX_MANIFEST_ENTRIES = 100000
MAX_ARCHIVE_BYTES = 128 * 1024**3
MAX_MEMBER_BYTES = 16 * 1024**3
MANIFEST_IDS = ("stratification_manifest", "truth_manifest", "stratification_member_manifest")


class AssetError(RuntimeError):
    """Acquisition gate failure with machine-readable evidence."""

    def __init__(self, code: str, message: str, **evidence):
        super().__init__(message)
        self.code = code
        self.evidence = evidence

    def as_dict(self) -> dict:
        return {"code": self.code, "message": str(self), "evidence": self.evidence}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _safe_name(value: str) -> str:
    if not isinstance(value, str) or not value or "\\" in value or ":" in value:
        raise AssetError("unsafe_path", "Invalid relative path", path=value)
    while value.startswith("./"):
        value = value[2:]
    parts = value.rstrip("/").split("/")
    if any(p in ("", ".", "..") for p in parts) or any(ord(c) < 32 for c in value):
        raise AssetError("unsafe_path", "Traversal or control character in path", path=value)
    if PurePosixPath(value).is_absolute():
        raise AssetError("unsafe_path", "Absolute path is forbidden", path=value)
    return "/".join(parts)


def _regular(path: Path) -> None:
    if not stat.S_ISREG(path.lstat().st_mode):
        raise AssetError("unsafe_file", "Expected a regular non-symlink file", path=str(path))


def _directory(path: Path) -> None:
    if path.is_symlink():
        raise AssetError("unsafe_directory", "Symlink directory is forbidden", path=str(path))
    path.mkdir(parents=True, exist_ok=True)
    if not path.is_dir():
        raise AssetError("unsafe_directory", "Expected directory", path=str(path))


def file_digest(path: Path) -> dict:
    """Return bytes, MD5, SHA-256 using constant-size reads (not cached metadata)."""
    path = Path(path)
    _regular(path)
    before = path.stat()
    md5, sha256, size = hashlib.md5(), hashlib.sha256(), 0
    with path.open("rb") as handle:
        while chunk := handle.read(CHUNK_BYTES):
            md5.update(chunk)
            sha256.update(chunk)
            size += len(chunk)
    after = path.stat()
    if (before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (
        after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns
    ) or size != after.st_size:
        raise AssetError("file_changed", "File changed while hashing", path=str(path))
    return {"bytes": size, "md5": md5.hexdigest(), "sha256": sha256.hexdigest()}


def _read_json(path: Path) -> dict:
    _regular(path)
    if path.stat().st_size > MAX_METADATA_BYTES:
        raise AssetError("metadata_limit", "Metadata exceeds bounded size", path=str(path))
    with path.open(encoding="utf-8") as handle:
        result = json.load(handle)
    if not isinstance(result, dict):
        raise AssetError("invalid_metadata", "Expected JSON object", path=str(path))
    return result


def _atomic_json(path: Path, value: dict) -> None:
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".part", dir=path.parent)
    temporary = Path(name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        _sync_directory(path.parent)
    finally:
        temporary.unlink(missing_ok=True)


def _sync_directory(path: Path) -> None:
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def load_lock() -> dict:
    """Load and check the adjacent versioned asset lock, never resolve new pins."""
    lock = _read_json(LOCK_PATH)
    if (lock.get("schema_version") != 1 or lock.get("sample") != "HG002"
            or lock.get("reference") != "GRCh38_no_alt"):
        raise AssetError("invalid_lock", "Unsupported lock schema/sample/reference")
    assets = lock.get("assets", {})
    if not isinstance(assets, dict):
        raise AssetError("invalid_lock", "Assets must be a mapping")
    required = {"bam", "bam_index", "reference", "truth", "truth_index", "benchmark_bed",
                "dbsnp", "dbsnp_index", "mills", "mills_index", "stratifications", *MANIFEST_IDS}
    if not required <= assets.keys():
        raise AssetError("invalid_lock", "Missing required asset ids", missing=sorted(required - assets.keys()))
    filenames = set()
    for asset_id, asset in assets.items():
        if not isinstance(asset, dict):
            raise AssetError("invalid_lock", "Asset must be an object", asset=asset_id)
        filename = _safe_name(asset["filename"])
        if "/" in filename or filename in filenames:
            raise AssetError("invalid_lock", "Asset filenames must be unique basenames", asset=asset_id)
        filenames.add(filename)
        url = urllib.parse.urlsplit(asset["url"])
        if url.scheme != "https" or not url.netloc or url.username or url.password or url.fragment:
            raise AssetError("invalid_lock", "Assets require public HTTPS URLs", asset=asset_id)
        md5 = asset.get("md5")
        if md5 is None:
            if asset_id not in MANIFEST_IDS or asset.get("checksum_policy") != "freeze_sha256_at_first_fetch":
                raise AssetError("invalid_lock", "Only official manifests may bootstrap a checksum", asset=asset_id)
        elif not re.fullmatch(r"[0-9a-f]{32}", md5):
            raise AssetError("invalid_lock", "Invalid locked MD5", asset=asset_id)
        if not isinstance(asset.get("large"), bool) or not asset.get("role"):
            raise AssetError("invalid_lock", "Missing role/large classification", asset=asset_id)
        if "size_bytes" in asset and (type(asset["size_bytes"]) is not int or asset["size_bytes"] <= 0):
            raise AssetError("invalid_lock", "Invalid expected object length", asset=asset_id)
        if "generation" in asset and not re.fullmatch(r"[0-9]+", str(asset["generation"])):
            raise AssetError("invalid_lock", "Invalid GCS generation", asset=asset_id)
        if "checksum_manifest" in asset:
            if asset["checksum_manifest"] not in MANIFEST_IDS:
                raise AssetError("invalid_lock", "Unknown checksum authority", asset=asset_id)
            _safe_name(asset["manifest_entry"])
    for tool in ("samtools", "bcftools", "happy"):
        if not re.fullmatch(r"[^\s]+@sha256:[0-9a-f]{64}", lock["images"][tool]):
            raise AssetError("invalid_lock", "Unpinned container image", tool=tool)
    return lock


def _lines(path: Path, *, compressed: bool = False):
    opener = gzip.open if compressed else open
    total = 0
    with opener(path, "rb") as handle:
        number = 0
        while raw := handle.readline(MAX_LINE_BYTES + 1):
            number += 1
            total += len(raw)
            if len(raw) > MAX_LINE_BYTES or total > MAX_MEMBER_BYTES:
                raise AssetError("text_limit", "Text line/member exceeds safety bound", path=str(path), line=number)
            yield number, raw.decode("utf-8").rstrip("\r\n")


def _parse_md5(path: Path) -> dict[str, str]:
    if path.stat().st_size > MAX_METADATA_BYTES:
        raise AssetError("manifest_limit", "Checksum manifest is too large", path=str(path))
    entries = {}
    for number, line in _lines(path):
        if not line or line.startswith("#"):
            continue
        match = re.fullmatch(r"([0-9a-fA-F]{32}) [ *](.+)", line)
        if not match:
            raise AssetError("invalid_manifest", "Malformed MD5 manifest entry", path=str(path), line=number)
        name = _safe_name(match[2])
        if name in entries:
            raise AssetError("duplicate_manifest_entry", "Duplicate checksum path", path=str(path), entry=name)
        entries[name] = match[1].lower()
        if len(entries) > MAX_MANIFEST_ENTRIES:
            raise AssetError("manifest_limit", "Too many checksum entries", path=str(path))
    if not entries:
        raise AssetError("invalid_manifest", "Empty checksum manifest", path=str(path))
    return entries


def _resolved_asset(asset: dict, manifests: dict) -> dict:
    resolved = dict(asset)
    if "checksum_manifest" in asset:
        authority = asset["checksum_manifest"]
        entry = _safe_name(asset["manifest_entry"])
        expected = manifests.get(authority, {}).get(entry)
        if expected is None:
            raise AssetError("manifest_entry_missing", "No official checksum for asset", manifest=authority, entry=entry)
        if asset.get("md5") is not None and asset["md5"] != expected:
            raise AssetError("publisher_checksum_changed", "Official manifest disagrees with lock", entry=entry,
                             locked_md5=asset["md5"], manifest_md5=expected)
        resolved["md5"] = expected
    return resolved


def _request_url(asset: dict) -> str:
    parts = urllib.parse.urlsplit(asset["url"])
    query = urllib.parse.parse_qsl(parts.query, keep_blank_values=True)
    if "generation" in asset:
        generation = str(asset["generation"])
        if any(key == "generation" and value != generation for key, value in query):
            raise AssetError("generation_conflict", "URL conflicts with locked generation")
        query = [(key, value) for key, value in query if key != "generation"]
        query.append(("generation", generation))
    return urllib.parse.urlunsplit(parts._replace(query=urllib.parse.urlencode(query)))


class _SameOriginRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        old, new = urllib.parse.urlsplit(req.full_url), urllib.parse.urlsplit(newurl)
        if new.scheme != "https" or new.netloc != old.netloc:
            raise AssetError("unsafe_redirect", "Refusing cross-origin or non-HTTPS redirect", url=newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _open_url(request):
    return urllib.request.build_opener(_SameOriginRedirect()).open(request, timeout=120)


def _check_digest(path: Path, digest: dict, asset: dict, receipt: dict | None = None) -> None:
    for expected_key, actual_key in (("size_bytes", "bytes"), ("md5", "md5")):
        expected = asset.get(expected_key)
        if expected is not None and digest[actual_key] != expected:
            raise AssetError("checksum_mismatch", "Asset does not match lock/manifest", path=str(path),
                             field=actual_key, expected=expected, observed=digest[actual_key])
    if receipt is not None:
        if receipt.get("url") != _request_url(asset):
            raise AssetError("receipt_identity_changed", "Receipt does not match requested object", path=str(path))
        for key in ("bytes", "md5", "sha256"):
            if receipt.get(key) != digest[key]:
                raise AssetError("frozen_content_changed", "Bytes differ from frozen acquisition receipt",
                                 path=str(path), field=key, expected=receipt.get(key), observed=digest[key])
    elif asset.get("md5") is None:
        raise AssetError("unfrozen_manifest", "Unpinned existing manifest has no frozen receipt", path=str(path))


def _receipt_path(path: Path) -> Path:
    return path.with_name(path.name + ".receipt.json")


def _verify_one(path: Path, asset: dict) -> dict:
    receipt_path = _receipt_path(path)
    receipt = _read_json(receipt_path) if receipt_path.exists() else None
    digest = file_digest(path)
    _check_digest(path, digest, asset, receipt)
    return {"status": "verified", "path": str(path), **digest, "url": _request_url(asset),
            "verified_at": _now(), "receipt": str(receipt_path) if receipt else None,
            "acquisition_provenance": receipt or {"origin": "externally_supplied_locked_bytes", "download_not_observed": True}}


def _reject_partial(part: Path, metadata: Path, error: AssetError) -> None:
    # Quarantine prevents another invocation from accepting or resuming bad bytes.
    suffix = ".rejected." + uuid.uuid4().hex
    rejected = []
    for path in (part, metadata):
        if path.exists() or path.is_symlink():
            destination = path.with_name(path.name + suffix)
            os.replace(path, destination)
            rejected.append(str(destination))
    if rejected:
        error.evidence["quarantined"] = rejected
        _sync_directory(part.parent)


def _transfer(path: Path, asset: dict, *, manifest: bool, small_limit: int) -> dict:
    if path.exists() or path.is_symlink():
        return _verify_one(path, asset)
    part = path.with_name(path.name + ".part")
    metadata_path = path.with_name(path.name + ".part.json")
    receipt_path = _receipt_path(path)
    url = _request_url(asset)
    expected = {"url": url, "md5": asset.get("md5"), "size_bytes": asset.get("size_bytes")}
    limit = MAX_METADATA_BYTES if manifest else (None if asset["large"] else small_limit)
    old_receipt = _read_json(receipt_path) if receipt_path.exists() else None
    offset = 0
    state = None
    if part.exists() or part.is_symlink():
        _regular(part)
        if not metadata_path.exists():
            raise AssetError("partial_identity_missing", "Cannot resume without object identity", path=str(part))
        state = _read_json(metadata_path)
        if state.get("expected") != expected:
            raise AssetError("partial_identity_changed", "Partial file belongs to different pinned object", path=str(part))
        offset = part.stat().st_size
        total = state.get("identity", {}).get("total_bytes")
        if type(total) is not int or offset > total:
            raise AssetError("invalid_partial_size", "Partial length exceeds or lacks object length", path=str(part))
    elif metadata_path.exists():
        raise AssetError("orphan_partial_metadata", "Partial metadata exists without bytes", path=str(metadata_path))

    started_at, started = _now(), time.monotonic()
    downloaded = 0
    if state is None or offset < state["identity"]["total_bytes"]:
        headers = {"Accept-Encoding": "identity", "User-Agent": "fastgatk-giab-assets/1"}
        if offset:
            identity = state["identity"]
            etag = identity.get("etag")
            if not etag and not asset.get("generation"):
                raise AssetError("resume_identity_unavailable", "Resume requires strong ETag or pinned GCS generation", path=str(part))
            headers["Range"] = f"bytes={offset}-"
            if etag:
                headers["If-Match"] = etag
        with _open_url(urllib.request.Request(url, headers=headers)) as response:
            status = response.status
            if response.headers.get("Content-Encoding", "identity").lower() != "identity":
                raise AssetError("encoded_transfer", "Server applied content encoding to asset", url=url)
            raw_etag = response.headers.get("ETag")
            etag = raw_etag if raw_etag and not raw_etag.startswith("W/") else None
            generation = response.headers.get("x-goog-generation")
            if asset.get("generation") and generation != str(asset["generation"]):
                raise AssetError("generation_changed", "Response did not prove locked GCS generation", url=url,
                                 expected=str(asset["generation"]), observed=generation)
            length = response.headers.get("Content-Length")
            if length is None or not re.fullmatch(r"[0-9]+", length):
                raise AssetError("length_missing", "Response requires an exact Content-Length", url=url)
            length = int(length)
            if offset:
                match = re.fullmatch(r"bytes ([0-9]+)-([0-9]+)/([0-9]+)", response.headers.get("Content-Range", ""))
                if status != 206 or match is None:
                    raise AssetError("invalid_range", "Resume requires 206 and exact Content-Range", url=url, status=status)
                first, last, total = map(int, match.groups())
                if first != offset or last != total - 1 or length != total - offset:
                    raise AssetError("invalid_range", "Response range is not the requested remaining bytes", url=url,
                                     content_range=match[0], offset=offset, length=length)
                identity = state["identity"]
                if total != identity["total_bytes"] or response.geturl() != identity["response_url"]:
                    raise AssetError("object_changed", "Object size or resolved URL changed during resume", url=url)
                if identity.get("etag") and identity["etag"] != etag:
                    raise AssetError("object_changed", "Strong ETag changed during resume", url=url)
                if identity.get("generation") != generation:
                    raise AssetError("object_changed", "Generation changed during resume", url=url)
            else:
                if status != 200 or response.headers.get("Content-Range") is not None:
                    raise AssetError("unexpected_partial_response", "Initial download requires complete 200 response", url=url, status=status)
                total = length
                if state is not None:
                    # Even a zero-length interrupted part must retain its identity.
                    identity = state["identity"]
                    if (identity.get("etag") != etag or identity.get("generation") != generation
                            or identity["total_bytes"] != total or identity["response_url"] != response.geturl()):
                        raise AssetError("object_changed", "Object changed since initial response", url=url)
            if total <= 0 or (asset.get("size_bytes") is not None and total != asset["size_bytes"]):
                raise AssetError("length_mismatch", "Response length disagrees with lock", url=url,
                                 expected=asset.get("size_bytes"), observed=total)
            if limit is not None and total > limit:
                raise AssetError("download_limit", "Object exceeds selected acquisition size class", url=url, bytes=total, limit=limit)
            if state is None:
                state = {"schema_version": 1, "expected": expected, "started_at": started_at,
                         "identity": {"etag": etag, "generation": generation, "total_bytes": total,
                                      "response_url": response.geturl(), "last_modified": response.headers.get("Last-Modified")}}
                _atomic_json(metadata_path, state)
            with part.open("ab" if offset else "wb") as handle:
                try:
                    while chunk := response.read(CHUNK_BYTES):
                        downloaded += len(chunk)
                        if downloaded > length:
                            raise AssetError("response_overrun", "Body exceeds declared response length", url=url)
                        handle.write(chunk)
                finally:
                    handle.flush()
                    os.fsync(handle.fileno())
            if downloaded != length:
                # A cleanly truncated HTTP body is resumable, but never publishable.
                raise OSError(f"Incomplete download of {url}: received {downloaded} of {length} response bytes")
    digest = file_digest(part)
    if digest["bytes"] != state["identity"]["total_bytes"]:
        raise AssetError("length_mismatch", "Complete partial differs from advertised length", path=str(part))
    if limit is not None and digest["bytes"] > limit:
        raise AssetError("download_limit", "Partial exceeds selected acquisition size class", path=str(part))
    # Bootstrap is only permitted for a freshly downloaded official manifest.
    if asset.get("md5") is None and old_receipt is None:
        if not manifest:
            raise AssetError("unlocked_asset", "Data assets require publisher checksums", path=str(part))
    else:
        _check_digest(part, digest, asset, old_receipt)
    if manifest:
        _parse_md5(part)
    receipt = {"schema_version": 1, "url": url, **digest, "origin": "https_download",
               "checksum_authority": asset.get("checksum_manifest", "lock" if asset.get("md5") else "HTTPS_first_fetch"),
               "identity": state["identity"], "started_at": state["started_at"], "completed_at": _now(),
               "last_attempt_started_at": started_at, "last_attempt_seconds": time.monotonic() - started,
               "last_attempt_downloaded_bytes": downloaded, "resumed_from_bytes": offset}
    # Receipt first: a crash cannot leave an unpinned manifest silently adoptable.
    _atomic_json(receipt_path, receipt)
    os.replace(part, path)
    metadata_path.unlink()
    _sync_directory(path.parent)
    return {"status": "verified", "path": str(path), **digest, "url": url,
            "verified_at": _now(), "receipt": str(receipt_path), "acquisition_provenance": receipt}


def _download_asset(path: Path, asset: dict, *, manifest: bool, small_limit: int) -> dict:
    try:
        return _transfer(path, asset, manifest=manifest, small_limit=small_limit)
    except urllib.error.HTTPError as error:
        if error.code in (412, 416):
            failure = AssetError("resume_rejected", "Server rejected object identity or requested range", status=error.code, url=error.url)
            _reject_partial(path.with_name(path.name + ".part"), path.with_name(path.name + ".part.json"), failure)
            raise failure from error
        raise
    except AssetError as error:
        _reject_partial(path.with_name(path.name + ".part"), path.with_name(path.name + ".part.json"), error)
        raise


def _stratification_checksums(lock: dict, manifests: dict) -> dict[str, str]:
    root = lock["stratification"]["archive_root"]
    checksums = {name: md5 for name, md5 in manifests["stratification_manifest"].items() if name.startswith(root + "/")}
    for name, md5 in manifests["stratification_member_manifest"].items():
        key = root + "/" + name
        if checksums.get(key) != md5:
            raise AssetError("stratification_manifest_conflict", "GRCh38 member and root manifests disagree", entry=key)
    entry = lock["stratification"]["entrypoint"]
    if checksums.get(entry) != lock["stratification"]["entrypoint_md5"]:
        raise AssetError("stratification_manifest_conflict", "Entrypoint differs from locked TSV", entry=entry)
    return checksums


def _bed_stats(path: Path) -> dict:
    records, chromosomes = 0, set()
    for number, line in _lines(path, compressed=path.name.endswith(".gz")):
        if not line or line.startswith(("#", "track ", "browser ")):
            continue
        fields = line.split("\t")
        if len(fields) < 3 or not fields[0] or not fields[1].isdigit() or not fields[2].isdigit():
            raise AssetError("invalid_stratification_bed", "Malformed BED record", path=str(path), line=number)
        start, end = int(fields[1]), int(fields[2])
        if start >= end:
            raise AssetError("invalid_stratification_bed", "BED requires 0 <= start < end", path=str(path), line=number)
        chromosomes.add(fields[0])
        if len(chromosomes) > MAX_MANIFEST_ENTRIES:
            raise AssetError("bed_limit", "Too many BED contigs", path=str(path))
        records += 1
    return {"records": records, "contigs": sorted(chromosomes), "reference_bounds": "not_checked"}


def _audit_stratifications(directory: Path, lock: dict, manifests: dict, *, frozen: dict | None = None) -> dict:
    checksums = _stratification_checksums(lock, manifests)
    members = {}
    for parent, dirs, files in os.walk(directory, followlinks=False):
        for name in dirs:
            if (Path(parent) / name).is_symlink():
                raise AssetError("unsafe_archive", "Extracted tree contains a symlink", path=str(Path(parent) / name))
        for name in files:
            path = Path(parent) / name
            relative = path.relative_to(directory).as_posix()
            expected = checksums.get(relative)
            if expected is None:
                raise AssetError("unlisted_archive_member", "Archive file has no official checksum", entry=relative)
            digest = file_digest(path)
            if digest["md5"] != expected:
                raise AssetError("archive_member_checksum", "Extracted member differs from official checksum", entry=relative,
                                 expected=expected, observed=digest["md5"])
            if frozen is not None and frozen.get("members", {}).get(relative, {}).get("sha256") != digest["sha256"]:
                raise AssetError("frozen_member_changed", "Extracted member differs from frozen SHA-256", entry=relative)
            if name.endswith((".bed", ".bed.gz")):
                digest["bed"] = _bed_stats(path)
            members[relative] = digest
            if len(members) > MAX_MANIFEST_ENTRIES:
                raise AssetError("archive_limit", "Too many extracted files")
    prefix = lock["stratification"]["archive_root"] + "/"
    required = {prefix + name for name in manifests["stratification_member_manifest"]}
    required.add(lock["stratification"]["entrypoint"])
    missing = required - members.keys()
    if missing:
        raise AssetError("archive_members_missing", "Archive lacks required official members", entries=sorted(missing))
    if frozen is not None and set(frozen.get("members", {})) != members.keys():
        raise AssetError("frozen_member_set_changed", "Extracted archive member set changed")
    entrypoint = directory / lock["stratification"]["entrypoint"]
    layers = []
    labels = set()
    for number, line in _lines(entrypoint):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 2 or not fields[0] or fields[0] in labels:
            raise AssetError("invalid_stratification_tsv", "Expected unique label and BED path", path=str(entrypoint), line=number)
        label, relative = fields[0], _safe_name(fields[1])
        key = prefix + relative
        if key not in members or not relative.endswith((".bed", ".bed.gz")):
            raise AssetError("invalid_stratification_tsv", "TSV refers to absent/non-BED member", entry=relative, line=number)
        if relative not in manifests["stratification_member_manifest"]:
            raise AssetError("unlisted_stratification", "TSV BED lacks a GRCh38 official checksum", entry=relative)
        labels.add(label)
        layers.append({"label": label, "member": key})
    if not layers:
        raise AssetError("invalid_stratification_tsv", "No stratifications in official TSV")
    return {"status": "verified", "members": members, "layers": layers,
            "bed_validation": "checksums, gzip integrity and coordinate syntax; reference bounds require biological preflight"}


class _BoundedTarInfo(tarfile.TarInfo):
    """Bound extension headers before tarfile allocates their declared bodies."""

    def _proc_member(self, archive):
        if self.type == tarfile.GNUTYPE_SPARSE:
            raise AssetError("unsafe_archive", "Sparse archive members are forbidden", entry=self.name)
        if self.type in (b"L", b"K", b"x", b"g", b"X") and not 0 <= self.size <= MAX_LINE_BYTES:
            raise AssetError("archive_limit", "Tar extension header exceeds bound", entry=self.name, bytes=self.size)
        return super()._proc_member(archive)

    def _proc_gnusparse_10(self, next_member, pax_headers, archive):
        # This PAX variant reads its sparse map before yielding a TarInfo.
        raise AssetError("unsafe_archive", "PAX sparse archive members are forbidden", entry=self.name)


def _prepare_stratifications(data_root: Path, lock: dict, manifests: dict, assets: dict, *, extract: bool) -> dict:
    archive = Path(assets["stratifications"]["path"])
    prepared = data_root / "prepared"
    directory = prepared / "stratifications"
    receipt_path = prepared / "stratifications-manifest.json"
    authority = {asset_id: assets[asset_id]["sha256"] for asset_id in ("stratifications", "stratification_manifest", "stratification_member_manifest")}
    started = time.monotonic()
    if directory.exists() or directory.is_symlink():
        if directory.is_symlink() or not directory.is_dir():
            raise AssetError("unsafe_archive", "Expected regular extracted directory", path=str(directory))
        if not receipt_path.exists():
            raise AssetError("unfrozen_archive", "Extracted tree has no acquisition audit", path=str(directory))
        frozen = _read_json(receipt_path)
        if frozen.get("source_sha256") != authority:
            raise AssetError("archive_source_changed", "Extracted tree belongs to different frozen assets")
        audit = _audit_stratifications(directory, lock, manifests, frozen=frozen)
    elif not extract:
        return {"status": "missing", "path": str(directory), "error": "Stratifications have not been safely extracted"}
    else:
        _directory(prepared)
        staging = Path(tempfile.mkdtemp(prefix="stratifications.part.", dir=prepared))
        try:
            checksums = _stratification_checksums(lock, manifests)
            prefix = lock["stratification"]["archive_root"]
            seen, expanded = set(), 0
            with tarfile.open(archive, "r|gz", tarinfo=_BoundedTarInfo) as handle:
                for member in handle:
                    name = _safe_name(member.name)
                    if name != prefix and not name.startswith(prefix + "/"):
                        raise AssetError("unsafe_archive", "Member lies outside declared archive root", entry=name)
                    if name in seen or len(seen) >= MAX_MANIFEST_ENTRIES:
                        raise AssetError("unsafe_archive", "Duplicate member or member-count limit", entry=name)
                    seen.add(name)
                    destination = staging / name
                    if member.isdir():
                        destination.mkdir(parents=True, exist_ok=True)
                    elif member.isfile() and not member.issparse():
                        if member.size < 0 or member.size > MAX_MEMBER_BYTES:
                            raise AssetError("archive_limit", "Member exceeds extraction bound", entry=name, bytes=member.size)
                        expanded += member.size
                        if expanded > MAX_ARCHIVE_BYTES:
                            raise AssetError("archive_limit", "Archive exceeds total extraction bound", bytes=expanded)
                        if name not in checksums:
                            raise AssetError("unlisted_archive_member", "Archive file has no official checksum", entry=name)
                        destination.parent.mkdir(parents=True, exist_ok=True)
                        source = handle.extractfile(member)
                        if source is None:
                            raise AssetError("invalid_archive", "Regular member has no body", entry=name)
                        with source, destination.open("xb") as output:
                            shutil.copyfileobj(source, output, CHUNK_BYTES)
                            output.flush()
                            os.fsync(output.fileno())
                        if destination.stat().st_size != member.size:
                            raise AssetError("invalid_archive", "Truncated member", entry=name)
                    else:
                        raise AssetError("unsafe_archive", "Links, devices, FIFOs and sparse files are forbidden", entry=name)
                    # tarfile otherwise retains every TarInfo, even in stream mode.
                    handle.members.clear()
            audit = _audit_stratifications(staging, lock, manifests)
            audit.update({"schema_version": 1, "source_sha256": authority, "created_at": _now(),
                          "tsv": str(directory / lock["stratification"]["entrypoint"]), "path": str(directory)})
            _atomic_json(receipt_path, audit)
            os.replace(staging, directory)
            _sync_directory(prepared)
        finally:
            if staging.exists():
                shutil.rmtree(staging)
    return {"status": "verified", "path": str(directory), "tsv": str(directory / lock["stratification"]["entrypoint"]),
            "manifest": str(receipt_path), "source_sha256": authority, "members": len(audit["members"]),
            "layers": len(audit["layers"]), "verified_at": _now(), "elapsed_seconds": time.monotonic() - started,
            "reference_eligibility": "pending"}


@contextmanager
def _acquisition_lock(data_root: Path):
    _directory(data_root)
    lock_path = data_root / ".giab-acquisition.lock"
    fd = os.open(lock_path, os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise AssetError("acquisition_busy", "Another fetch/verify owns this data root", path=str(data_root)) from error
        yield
    finally:
        os.close(fd)


def _error_dict(error: Exception) -> dict:
    if isinstance(error, AssetError):
        return error.as_dict()
    return {"code": "io_error" if isinstance(error, OSError) else "invalid_asset",
            "message": str(error), "exception": type(error).__name__}


def _event(data_root: Path, value: dict) -> None:
    path = data_root / "acquisition-events.jsonl"
    fd = os.open(path, os.O_WRONLY | os.O_APPEND | os.O_CREAT | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, "a", encoding="utf-8") as handle:
        handle.write(json.dumps(value, sort_keys=True, allow_nan=False) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def _acquire(data_root: Path, *, fetch: bool, large: bool, only: set[str] | None = None) -> dict:
    data_root = Path(data_root).absolute()
    started_at, started = _now(), time.monotonic()
    report = {"schema_version": 1, "operation": "fetch" if fetch else "verify", "large_enabled": large,
              "data_root": str(data_root), "started_at": started_at, "status": "blocked", "requested_status": "failed",
              "assets": {}, "errors": [], "stratifications": {"status": "not_evaluated"},
              "reference_eligibility": "pending"}
    try:
        lock = load_lock()
        report.update({"sample": lock["sample"], "reference": lock["reference"],
                       "lock": {"path": str(LOCK_PATH), **file_digest(LOCK_PATH)}})
        selected = set(lock["assets"]) if only is None else set(only)
        if selected - lock["assets"].keys():
            raise AssetError("unknown_asset_selection", "Unknown selected asset ids",
                             assets=sorted(selected - lock["assets"].keys()))
        selected.update(MANIFEST_IDS)
        report["requested_assets"] = sorted(selected)
        with _acquisition_lock(data_root):
            raw = data_root / "raw"
            _directory(raw)
            manifests = {}
            entries = lock["assets"]
            ordered = list(MANIFEST_IDS) + [key for key, asset in entries.items() if key not in MANIFEST_IDS and not asset["large"]]
            # Materialize/reference-check assets before the expensive full BAM.
            ordered += [key for key, asset in entries.items() if asset["large"] and key != "bam"]
            ordered.append("bam")
            for asset_id in ordered:
                asset = entries[asset_id]
                path = raw / asset["filename"]
                if asset_id not in selected:
                    report["assets"][asset_id] = {"status": "not_selected", "path": str(path)}
                    continue
                if fetch and asset["large"] and not large:
                    report["assets"][asset_id] = {"status": "deferred", "path": str(path), "reason": "requires_explicit_large_fetch"}
                    continue
                if fetch and report["errors"]:
                    report["assets"][asset_id] = {"status": "not_attempted", "path": str(path), "reason": "earlier_acquisition_gate_failed"}
                    continue
                item_started = time.monotonic()
                try:
                    asset = _resolved_asset(asset, manifests)
                    if not fetch and not path.exists():
                        raise AssetError("asset_missing", "Required raw asset is absent", asset=asset_id, path=str(path))
                    item = (_download_asset(path, asset, manifest=asset_id in MANIFEST_IDS,
                                            small_limit=lock["acquisition"]["small_asset_limit_bytes"])
                            if fetch else _verify_one(path, asset))
                    if asset_id in MANIFEST_IDS:
                        manifests[asset_id] = _parse_md5(path)
                    if asset_id == MANIFEST_IDS[-1]:
                        # Resolve every locked checksum before downloading payloads.
                        for dependent in entries.values():
                            _resolved_asset(dependent, manifests)
                        _stratification_checksums(lock, manifests)
                    item["elapsed_seconds"] = time.monotonic() - item_started
                    report["assets"][asset_id] = item
                    if asset_id == "stratifications":
                        report["stratifications"] = _prepare_stratifications(
                            data_root, lock, manifests, report["assets"], extract=fetch)
                        if report["stratifications"]["status"] != "verified":
                            raise AssetError("stratifications_missing", "Archive also requires audited extracted members")
                except (AssetError, OSError, ValueError, KeyError, TypeError, http.client.HTTPException, tarfile.TarError, EOFError) as error:
                    detail = {"asset": asset_id, **_error_dict(error)}
                    report["errors"].append(detail)
                    report["assets"][asset_id] = {"status": "missing" if detail["code"] == "asset_missing" else "failed",
                                                "path": str(path), "error": detail, "elapsed_seconds": time.monotonic() - item_started}
                    if asset_id == "stratifications":
                        report["stratifications"] = {"status": "failed", "error": detail}
                _event(data_root, {"at": _now(), "operation": report["operation"], "asset": asset_id, **report["assets"][asset_id]})
            if not report["errors"]:
                report["requested_status"] = "verified"
                report["status"] = "partial" if any(item["status"] in ("deferred", "not_selected")
                                                  for item in report["assets"].values()) else "complete"
            report.update({"completed_at": _now(), "elapsed_seconds": time.monotonic() - started,
                           "manifest": str(data_root / "asset-manifest.json")})
            _atomic_json(data_root / "asset-manifest.json", report)
            _event(data_root, {"at": _now(), "operation": report["operation"], "status": report["status"],
                               "elapsed_seconds": report["elapsed_seconds"], "errors": report["errors"]})
    except (AssetError, OSError, ValueError, KeyError, TypeError) as error:
        report["errors"].append(_error_dict(error))
        report.update({"status": "blocked", "requested_status": "failed", "completed_at": _now(),
                       "elapsed_seconds": time.monotonic() - started})
    return report


def fetch_assets(data_root: Path, *, large: bool = False, only: set[str] | None = None) -> dict:
    """Fetch manifests, small resources, then explicitly authorized large assets.

    only selects payload ids; official manifests are always acquired first.
    Large order is reference, known sites, stratifications (including audit), BAM.
    status=partial with requested_status=verified is successful small/selected fetch;
    status=complete requires all raw assets and audited extracted stratifications.
    Any failure is status=blocked with structured errors, never a success fallback.
    """
    return _acquire(data_root, fetch=True, large=large, only=only)


def verify_assets(data_root: Path) -> dict:
    """Offline rehash of every raw/stratification asset; never downloads or repairs."""
    return _acquire(data_root, fetch=False, large=True)
