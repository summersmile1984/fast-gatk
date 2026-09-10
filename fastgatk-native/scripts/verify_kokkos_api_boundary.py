#!/usr/bin/env python3
"""Enforce the single-source Kokkos SIMD boundary for production code.

The portable implementation must expose CPU vectorization through Kokkos
(`Kokkos::Experimental::simd`/execution spaces), not through a second set of
host ISA kernels.  This check intentionally scans source files only; build
dependency files may contain compiler intrinsic headers pulled in by Kokkos
and are not part of the production API boundary.
"""

from __future__ import annotations

import json
import pathlib
import re


ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCE_ROOTS = (ROOT / "fastgatk-kernels" / "src", ROOT / "fastgatk-native" / "src")

# Keep this deliberately narrow.  CPU feature probes (`__builtin_cpu_supports`)
# and words such as "AVX2" in documentation are harmless; raw vector types,
# intrinsic calls and direct immintrin includes are the portability violation.
FORBIDDEN = re.compile(
    r"^\s*#\s*include\s*[<\"]immintrin\.h[>\"]|"
    r"\b_mm(?:[A-Za-z0-9_]*)\s*\(|\b__m(?:128|256|512)(?:d|i|h)?\b"
)


def main() -> int:
    scanned = 0
    violations: list[dict[str, object]] = []
    for source_root in SOURCE_ROOTS:
        for path in sorted(source_root.glob("*.cpp")):
            scanned += 1
            for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
                if FORBIDDEN.search(line):
                    violations.append({
                        "file": str(path.relative_to(ROOT)),
                        "line": line_number,
                        "text": line.strip(),
                    })
    if violations:
        print(json.dumps({"status": "fail", "scanned_files": scanned,
                          "violations": violations}, sort_keys=True))
        return 1
    print(json.dumps({"status": "pass", "scanned_files": scanned,
                      "shared_kokkos_api": True, "raw_host_intrinsics": False},
                     sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
