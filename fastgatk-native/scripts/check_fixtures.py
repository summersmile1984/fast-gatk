#!/usr/bin/env python3
"""Fail-loud fixture presence check for the gitignored corpora.

The tracked small fixtures are already pinned by
``pinned_fixture_digests.sha256`` and verified by ``verify_fixture_digests.py``.
This script covers the *untracked* corpora (``fixtures/``, ``testdata/``,
``gatk-source/src/test/resources/``) that the oracle and benchmark scripts
reference as repository-relative paths.  It exists so a fresh machine fails
with an explicit missing-file list instead of a confusing downstream error.

Usage
-----
    python3 fastgatk-native/scripts/check_fixtures.py            # exit 0 ok, 1 missing
    FASTGATK_REQUIRE_GATK_ORACLE=1 python3 .../check_fixtures.py  # (same; always fail-loud)
"""

from __future__ import annotations

import ast
import glob
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "fastgatk-native" / "scripts"

_EXTENSIONS = (
    ".bam", ".bai", ".cram", ".crai", ".sam", ".fa", ".fasta", ".fai", ".dict",
    ".vcf", ".vcf.gz", ".tbi", ".bed", ".bed.gz", ".interval_list", ".tsv",
    ".tranches", ".csv", ".stats", ".pon", ".g.vcf", ".g.vcf.gz", ".hd5", ".hdf5",
    ".gz", ".idx",
)
_PREFIXES = ("fixtures/", "testdata/", "gatk-source/src/test/resources/")


def _string_literals(path: Path):
    try:
        tree = ast.parse(path.read_text(encoding="utf-8", errors="replace"))
    except Exception:
        return
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            yield node.value


def referenced_fixtures() -> dict[str, set[str]]:
    refs: dict[str, set[str]] = {}
    for script in sorted(glob.glob(str(SCRIPTS / "verify_*.py"))) + \
            sorted(glob.glob(str(SCRIPTS / "benchmark_*.py"))):
        for literal in _string_literals(Path(script)):
            for prefix in _PREFIXES:
                # repo-relative fixture literals start at the beginning of the
                # string (e.g. "fixtures/chr20/mnp.bam"); a match deeper inside
                # "fastgatk-native/tests/fixtures/..." is a bundled (tracked)
                # fixture, not a gitignored corpus path.
                if not literal.startswith(prefix):
                    continue
                candidate = literal.split()[0].rstrip('",\';()')
                if not candidate.lower().endswith(_EXTENSIONS):
                    continue
                # skip path navigation and source-code citations
                if ".." in candidate or ":" in candidate.split("/")[-1]:
                    continue
                refs.setdefault(prefix.rstrip("/"), set()).add(candidate)
    return refs


def main() -> int:
    missing: dict[str, list[str]] = {}
    present = 0
    for group, paths in sorted(referenced_fixtures().items()):
        for rel in sorted(paths):
            if (ROOT / rel).exists():
                present += 1
            else:
                missing.setdefault(group, []).append(rel)
    print(f"fixture check: {present} present, "
          f"{sum(len(v) for v in missing.values())} missing")
    if missing:
        for group, paths in sorted(missing.items()):
            print(f"[{group}] missing:")
            for rel in paths:
                print(f"    {rel}")
        print("\nRestore these corpora (see FIXTURE_CORPUS.md and each dir's "
              "MANIFEST.md / fetch_real_testdata.sh) before running the suite.")
        return 1
    print("all referenced gitignored fixtures are present")
    return 0


if __name__ == "__main__":
    sys.exit(main())
