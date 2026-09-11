"""Oracle-presence guard: make "no GATK oracle" impossible to mistake for parity.

Why this module exists
----------------------
Many `verify_*.py` scripts wrap their GATK comparison in a presence guard shaped
like::

    java = root / "third_party/jdk17/bin/java"
    jar  = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    if java.exists() and jar.exists():
        ... run GATK, compare against native output ...
    # (implicit else: nothing happens)

and the dedicated `*_gatk_oracle.py` scripts do the mirror image::

    if not JAVA.exists() or not GATK.exists():
        print("skip: bundled GATK/JDK oracle not present")
        return 0

Neither form ever *asserts* that the oracle exists.  When the jar or the JDK
disappears, the test keeps exiting 0 and CTest keeps printing ``Passed``: the
script has silently degraded into "native output compared against native
expectations recorded earlier".  Green then no longer means matching GATK.

Policy implemented here
-----------------------
* oracle present                          -> the caller behaves exactly as before;
* oracle missing, ``FASTGATK_REQUIRE_GATK_ORACLE`` set
                                          -> raise :class:`GATKOracleMissing`
                                             (non-zero exit) naming the missing path;
* oracle missing, that variable unset     -> print an unmistakable
                                             "NOT verified against GATK" notice and
                                             let the caller keep its old skip path.

Usage
-----
Positive guard (the comparison only runs when the oracle is there)::

    from oracle_guard import oracle_ready
    if oracle_ready("verify_filter_intervals.py", java, jar):
        ... run GATK, compare ...

Negative skip branch (the test already decided to bail out)::

    from oracle_guard import oracle_not_verified
    if not NATIVE.exists() or not JAVA.exists() or not GATK.exists():
        oracle_not_verified("verify_x.py", JAVA, GATK)   # no-op if only NATIVE is missing
        print("skip: ...")
        return 0

Both helpers ignore paths they were not given: ``oracle_not_verified`` stays
silent when every oracle path it was handed exists, so a skip caused by a missing
native binary or fixture is not mislabelled as a missing oracle.

Note on visibility: CTest runs here with ``--output-on-failure``, which discards
the output of *passing* tests.  The notice therefore reaches the regression log
only when a single test is re-run verbosely (``ctest -V -R <name>``) or the script
is run by hand.  The notice is never suppressed, only captured.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Iterable, List

REQUIRE_ENV = "FASTGATK_REQUIRE_GATK_ORACLE"

#: Values of REQUIRE_ENV that mean "do not require" (anything else counts as set).
_DISABLED = {"", "0", "false", "no", "off"}

_RULE = "=" * 78


class GATKOracleMissing(RuntimeError):
    """Raised when the oracle is required but absent.  Uncaught -> exit code 1."""


def requires_oracle() -> bool:
    """True when ``FASTGATK_REQUIRE_GATK_ORACLE`` is set to a truthy value."""
    return os.environ.get(REQUIRE_ENV, "").strip().lower() not in _DISABLED


def missing_paths(paths: Iterable[object]) -> List[str]:
    """Return the entries of ``paths`` that are not existing regular files."""
    out: List[str] = []
    for p in paths:
        try:
            candidate = Path(p)
        except TypeError:
            continue
        if not candidate.is_file():
            out.append(str(candidate))
    return out


def _reported(paths: Iterable[object]) -> List[str]:
    """The subset of ``paths`` that is missing.

    Callers are responsible for passing only *oracle* paths (the JDK binary and
    the pinned GATK jar); every path handed in is therefore reported verbatim.
    No name-based filtering happens here on purpose: a ``JAVA=/bogus/java``
    override at run time yields a path that looks nothing like the vendored JDK,
    yet it is still the oracle slot that failed to resolve.
    """
    return missing_paths(paths)


def _raise(test: str, missing: List[str]) -> None:
    raise GATKOracleMissing(
        "\n".join(
            [
                _RULE,
                f"[ORACLE REQUIRED] {test}: {REQUIRE_ENV}={os.environ.get(REQUIRE_ENV, '')!r}",
                "The GATK oracle is required for this test but was not found:",
                *[f"    missing: {m}" for m in missing],
                "",
                "The GATK comparison is the whole point of this test, so it fails loudly",
                "instead of passing without it.  Unset the variable to fall back to the old",
                "skip-and-pass behaviour (which prints a NOT-VERIFIED notice).",
                _RULE,
            ]
        )
    )


_NOTIFIED: set = set()


def _notice(test: str, missing: List[str]) -> None:
    # A script may guard several comparison sites; one notice per script is
    # enough to make "this run was not GATK-verified" impossible to miss.
    if test in _NOTIFIED:
        return
    _NOTIFIED.add(test)
    text = "\n".join(
        [
            _RULE,
            f"[NOT VERIFIED AGAINST GATK] {test}",
            "The GATK oracle is absent, so the expectations asserted by this test were",
            "NOT verified against GATK 4.6.2.0.  A PASS from this script is NOT evidence",
            "of GATK parity -- it only means native output matched native expectations",
            "recorded earlier.",
            "",
            *[f"    missing: {m}" for m in missing],
            "",
            f"Set {REQUIRE_ENV}=1 to turn this into a hard failure.",
            _RULE,
        ]
    )
    sys.stderr.write(text + "\n")
    sys.stderr.flush()


def oracle_ready(test: str, *paths: object) -> bool:
    """Oracle-presence test to use in place of ``java.exists() and jar.exists()``.

    Returns True when every path given exists (caller proceeds exactly as before).
    Otherwise raises (env var set) or prints the NOT-VERIFIED notice and returns
    False (env var unset).
    """
    missing = missing_paths(paths)
    if not missing:
        return True
    if requires_oracle():
        _raise(test, missing)
    _notice(test, missing)
    return False


def oracle_not_verified(test: str, *paths: object) -> None:
    """Annotate an already-taken "skip because something is missing" branch.

    Silent no-op unless one of the given *oracle* paths is really absent, so a
    branch entered because the native binary or a fixture is missing keeps its old
    behaviour.
    """
    reported = _reported(paths)
    if not reported:
        return
    if requires_oracle():
        _raise(test, reported)
    _notice(test, reported)


if __name__ == "__main__":  # manual self-check: python3 oracle_guard.py <jar> <java>
    _args = sys.argv[1:]
    if not _args:
        _args = [
            str(Path(__file__).resolve().parents[2] / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"),
            str(Path(__file__).resolve().parents[2] / "third_party/jdk17/bin/java"),
        ]
    print(f"{REQUIRE_ENV} set: {requires_oracle()}")
    for _p in _args:
        print(f"  {_p}: exists={Path(_p).is_file()}")
    print(f"oracle_ready -> {oracle_ready('oracle_guard.py (self-check)', *_args)}")
