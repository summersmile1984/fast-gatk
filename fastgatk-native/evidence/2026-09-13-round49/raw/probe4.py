#!/usr/bin/env python3
"""Round-49 batch 4: controls for the '*' ownership/placement fixture and a
DP=0 variant of it."""
from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import probe  # noqa: E402
import probe3  # noqa: E402

DENSE = probe.DENSE
HET2, HOME_STAR = probe.HET2, probe.HOM_STAR

CASES: list[dict] = []


def case(name, why, body, args=None, compare_positions=None):
    CASES.append(dict(name=name, why=why, body=body, reference=probe.AAA,
                      args=args or [], header=probe.HEADER, native_args=[],
                      compare_positions=compare_positions))


case("g-star-nomid-control",
     "control for g-star-chain-dense: the intermediate '*' record at 4 is "
     "removed, so the only recorded deletion is record 1's [2,6] span, which "
     "does NOT reach position 7.  The '*' at 7 must then be treated as a "
     "spurious spanning deletion by BOTH tools (covering the negative case of "
     "the coverage test that keeps row 7 alive in the chain fixture)",
     "chr1\t2\t.\tAAAAAA\tAA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HET2 + "\n"
     "chr1\t7\t.\tAAAA\t*,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HOME_STAR + "\n",
     args=[DENSE], compare_positions=[2, 7])

case("g-star-mid-only-control",
     "control: only record 1 and the intermediate '*' record at 4; the row at 4 "
     "must survive on both sides",
     "chr1\t2\t.\tAAAAAA\tAA,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HET2 + "\n"
     "chr1\t4\t.\tAAAAA\t*,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HOME_STAR + "\n",
     args=[DENSE], compare_positions=[2, 4])

case("dp0-star-ownership-dense",
     "the DP=0 reach gap combined with '*' ownership: with INFO/DP=0 GATK takes "
     "the `result = originalVC` arm, so GenotypingEngine.recordDeletions() is "
     "never reached and the downstream '*' has no covering deletion at all",
     "chr1\t2\t.\tAAAA\tAA,<NON_REF>\t.\tPASS\tDP=0\tGT:DP:AD:PL\t" + HET2 + "\n"
     "chr1\t4\t.\tAA\t*,<NON_REF>\t.\tPASS\tDP=20\tGT:DP:AD:PL\t" + HOME_STAR + "\n",
     args=[DENSE], compare_positions=[2, 4])


def main() -> int:
    selected = [c for c in CASES if not sys.argv[1:] or c["name"] in sys.argv[1:]]
    results = {}
    for c in selected:
        saved = probe3.CASES
        probe3.CASES = [c]
        results[c["name"]] = probe3.run_case(c)
        probe3.CASES = saved
    print("=" * 100)
    for name, ok in results.items():
        print(f"{'EQUAL ' if ok else 'DIFFER'}  {name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
