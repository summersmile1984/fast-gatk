#!/usr/bin/env python3
"""Verify that the pinned oracle/benchmark fixtures have not drifted.

The digest gate binds GATK/native comparisons to exact input bytes.  It does
not claim that a digest proves algorithmic equivalence; that remains the job
of each tool-specific oracle.  The manifest is deliberately checked from the
source tree so CTest and standalone invocations use the same corpus.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "fastgatk-native/tests/pinned_fixture_digests.sha256"


def read_manifest() -> list[tuple[str, str]]:
    if not MANIFEST.is_file():
        raise SystemExit(f"missing fixture digest manifest: {MANIFEST}")
    entries: list[tuple[str, str]] = []
    seen: set[str] = set()
    for line_number, raw in enumerate(MANIFEST.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split(maxsplit=1)
        if len(fields) != 2 or len(fields[0]) != 64:
            raise SystemExit(f"invalid digest manifest line {line_number}: {raw!r}")
        digest, relative = fields
        relative = relative.lstrip("*").strip()
        path = Path(relative)
        if path.is_absolute() or ".." in path.parts:
            raise SystemExit(f"manifest path escapes repository at line {line_number}: {relative}")
        if relative in seen:
            raise SystemExit(f"duplicate fixture path at line {line_number}: {relative}")
        seen.add(relative)
        entries.append((digest.lower(), relative))
    if not entries:
        raise SystemExit("fixture digest manifest is empty")
    return entries


def digest(path: Path) -> tuple[str, int]:
    hasher = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(chunk)
            size += len(chunk)
    return hasher.hexdigest(), size


def main() -> int:
    results = []
    for expected, relative in read_manifest():
        path = ROOT / relative
        if not path.is_file():
            raise SystemExit(f"missing pinned fixture: {relative}")
        actual, size = digest(path)
        if actual != expected:
            raise SystemExit(
                f"fixture digest mismatch: {relative}: expected {expected}, got {actual}"
            )
        results.append({"path": relative, "bytes": size, "sha256": actual})
    print(json.dumps({
        "schema_version": 1,
        "status": "pass",
        "manifest": str(MANIFEST.relative_to(ROOT)),
        "files": results,
        "count": len(results),
        "total_bytes": sum(item["bytes"] for item in results),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
