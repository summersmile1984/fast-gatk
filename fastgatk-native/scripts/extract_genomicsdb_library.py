#!/usr/bin/env python3
"""Extract only libtiledbgenomicsdb.so from the pinned GATK package.

The GATK distribution is a large ZIP/JAR.  Keeping this extraction explicit
and atomic avoids unpacking the whole package during a native build and never
leaves a partially written shared object that could be linked by a parallel
build.
"""

from __future__ import annotations

import argparse
import os
import shutil
import tempfile
import zipfile
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", required=True, type=Path)
    parser.add_argument("--member", default="libtiledbgenomicsdb.so")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    if not args.package.is_file():
        raise SystemExit(f"GATK package is missing: {args.package}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.package) as archive:
        try:
            source = archive.open(args.member, "r")
        except KeyError as exc:
            raise SystemExit(f"GATK package does not contain {args.member}") from exc
        fd, temporary_name = tempfile.mkstemp(
            prefix=f".{args.output.name}.", dir=args.output.parent
        )
        try:
            with os.fdopen(fd, "wb") as destination, source:
                shutil.copyfileobj(source, destination)
                destination.flush()
                os.fsync(destination.fileno())
            os.chmod(temporary_name, 0o755)
            os.replace(temporary_name, args.output)
            # The packaged object carries SONAME libtiledbgenomicsdb.so.1,
            # so the runtime loader needs the conventional sibling symlink
            # even though the ZIP member is named libtiledbgenomicsdb.so.
            soname = args.output.with_name(args.output.name + ".1")
            try:
                soname.unlink()
            except FileNotFoundError:
                pass
            soname.symlink_to(args.output.name)
        finally:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
