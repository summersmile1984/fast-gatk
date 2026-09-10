#!/usr/bin/env python3
"""Pinned GATK oracle for CollectAllelicCounts read-filter resolution.

This slice is intentionally narrower than full CNV parity.  It proves that
the native pileup applies the same resolved default list as GATK 4.6.2.0 and
that ``--disable-tool-default-read-filters`` removes that list before explicit
filters are added.  Secondary/supplementary/QC-fail flags are included because
they are *not* CollectAllelicCounts defaults in this release.
"""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
BUILD = pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native/build"))
NATIVE = pathlib.Path(os.environ.get(
    "FASTGATK_COLLECT_ALLELIC_COUNTS_BINARY",
    str(BUILD / "fastgatk-collect-allelic-counts"),
))
JAVA = ROOT / "third_party/jdk17/bin/java"
JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


def rows(path: pathlib.Path) -> list[tuple[str, int, int, int, str, str]]:
    result: list[tuple[str, int, int, int, str, str]] = []
    for line in path.read_text(encoding="ascii").splitlines():
        if not line or line.startswith("@") or line.startswith("CONTIG"):
            continue
        fields = line.split("\t")
        result.append((fields[0], int(fields[1]), int(fields[2]), int(fields[3]),
                       fields[4], fields[5]))
    return result


def write_reference(path: pathlib.Path) -> None:
    sequence = ("ACGT" * 30)[:100]
    path.write_text(f">chr1\n{sequence}\n", encoding="ascii")
    path.with_suffix(path.suffix + ".fai").write_text(
        f"chr1\t{len(sequence)}\t6\t{len(sequence)}\t{len(sequence) + 1}\n",
        encoding="ascii",
    )


def write_sam(path: pathlib.Path) -> None:
    bases = "TTTT"
    quality = "IIII"
    records = [
        ("normal", 0, 60),
        ("duplicate", 1024, 60),
        ("low-mapq", 0, 20),
        ("secondary", 256, 60),
        ("supplementary", 2048, 60),
        ("qcfail", 512, 60),
    ]
    lines = [
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chr1\tLN:100",
        "@RG\tID:rg1\tSM:SAMPLE",
    ]
    for name, flags, mapq in records:
        lines.append(
            f"{name}\t{flags}\tchr1\t1\t{mapq}\t4M\t*\t0\t0\t{bases}\t{quality}\tRG:Z:rg1"
        )
    # A spliced read is rejected by GATK's default WellformedReadFilter
    # (CigarContainsNoNOperator), but contributes its M bases when all tool
    # defaults are disabled.  Keep it in the same indexed traversal so this
    # checks the filter boundary rather than a separate input mode.
    lines.append(
        "spliced\t0\tchr1\t1\t60\t2M2N2M\t*\t0\t0\tTTTT\tIIII\tRG:Z:rg1"
    )
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def run_native(extra: list[str], bam: pathlib.Path, reference: pathlib.Path,
               work: pathlib.Path, label: str) -> tuple[list[tuple], dict]:
    output = work / f"native-{label}.tsv"
    manifest = work / f"native-{label}.json"
    command = [str(NATIVE), "-I", str(bam), "-R", str(reference), "-L", "chr1:1-4",
               "-O", str(output), "--sample", "SAMPLE", "--output-manifest", str(manifest),
               "--batch-records", "2"] + extra
    result = subprocess.run(command, text=True, capture_output=True, env=os.environ.copy())
    if result.returncode != 0:
        raise AssertionError(f"native {label} failed:\n{result.stdout}\n{result.stderr}")
    summary = json.loads(result.stdout.splitlines()[-1])
    return rows(output), json.loads(manifest.read_text(encoding="ascii"))


def run_java(extra: list[str], bam: pathlib.Path, reference: pathlib.Path,
             work: pathlib.Path, label: str) -> list[tuple]:
    output = work / f"java-{label}.tsv"
    command = [str(JAVA), "-jar", str(JAR), "CollectAllelicCounts", "-I", str(bam),
               "-R", str(reference), "-L", "chr1:1-4", "-O", str(output)] + extra
    result = subprocess.run(command, text=True, capture_output=True, env=os.environ.copy())
    if result.returncode != 0:
        raise AssertionError(f"GATK {label} failed:\n{result.stdout}\n{result.stderr[-5000:]}")
    return rows(output)


def compare_case(label: str, native_extra: list[str], java_extra: list[str],
                 bam: pathlib.Path, reference: pathlib.Path, work: pathlib.Path,
                 expected_used: int, expected_mapping_quality: bool) -> None:
    native_rows, manifest = run_native(native_extra, bam, reference, work, label)
    java_rows = run_java(java_extra, bam, reference, work, label)
    if native_rows != java_rows:
        raise AssertionError({"label": label, "native": native_rows, "java": java_rows})
    if manifest["reads_used"] != expected_used:
        raise AssertionError({"label": label, "manifest": manifest})
    if manifest["read_filter_policy"]["mapping_quality"] != expected_mapping_quality:
        raise AssertionError({"label": label, "manifest": manifest})


def main() -> int:
    if not NATIVE.exists():
        raise SystemExit(f"missing native binary: {NATIVE}")
    if not JAVA.exists() or not JAR.exists():
        raise SystemExit("missing pinned GATK 4.6.2.0 runtime")
    with tempfile.TemporaryDirectory(prefix="fastgatk-collect-allelic-gatk-") as directory:
        work = pathlib.Path(directory)
        reference = work / "reference.fasta"
        sam = work / "reads.sam"
        write_reference(reference)
        dict_result = subprocess.run([
            str(JAVA), "-jar", str(JAR), "CreateSequenceDictionary",
            f"R={reference}", f"O={reference.with_suffix('.dict')}",
        ], text=True, capture_output=True, env=os.environ.copy())
        if dict_result.returncode != 0:
            raise AssertionError("Picard CreateSequenceDictionary failed:\n" + dict_result.stderr[-5000:])
        write_sam(sam)
        bam = work / "reads.bam"
        sort_result = subprocess.run([
            str(JAVA), "-jar", str(JAR), "SortSam", f"I={sam}", f"O={bam}",
            "SORT_ORDER=coordinate", "CREATE_INDEX=true", "VALIDATION_STRINGENCY=LENIENT",
        ], text=True, capture_output=True, env=os.environ.copy())
        if sort_result.returncode != 0:
            raise AssertionError("Picard SortSam failed:\n" + sort_result.stderr[-5000:])
        bam_index = pathlib.Path(str(bam) + ".bai")
        if not bam_index.exists():
            bam_index = bam.with_suffix(".bai")
        if not bam_index.exists():
            raise AssertionError("Picard SortSam did not create BAM index")

        # GATK's LocusWalker defaults are Wellformed + Mapped, and this tool
        # adds Mapped + NonZeroReferenceLength + NotDuplicate + MAPQ>=30.
        # Secondary, supplementary and QC-fail reads are retained.
        compare_case("default", [], [], bam, reference, work, expected_used=4,
                     expected_mapping_quality=True)

        # The advanced Boolean accepts both the explicit GATK form and native
        # optional-Boolean spelling.  Every default filter is removed here;
        # the interval traversal still exposes the six mapped records.
        compare_case("no-defaults", ["--disable-tool-default-read-filters", "true"],
                     ["--disable-tool-default-read-filters", "true"],
                     bam, reference, work, expected_used=7,
                     expected_mapping_quality=False)

        # Explicit read filters are resolved after disable-tool-defaults,
        # independent of command-line order.  Re-enable only duplicate and
        # mapping-quality filters and compare exact pileup rows.
        compare_case(
            "explicit",
            ["--disable-tool-default-read-filters", "true", "--read-filter",
             "NotDuplicateReadFilter", "--read-filter", "MappingQualityReadFilter"],
            ["--disable-tool-default-read-filters", "true", "--read-filter",
             "NotDuplicateReadFilter", "--read-filter", "MappingQualityReadFilter",
             "--minimum-mapping-quality", "30"],
            bam, reference, work, expected_used=5, expected_mapping_quality=True,
        )

        # -XL applies to the locus traversal itself.  The native path must
        # remove these sites from its device locus table while still decoding
        # a read that crosses the excluded span; treating it as a whole-read
        # filter would change neighboring pileups.
        excluded_native, excluded_manifest = run_native(
            ["-XL", "chr1:2-3"], bam, reference, work, "excluded")
        excluded_java = run_java(["-XL", "chr1:2-3"], bam, reference, work, "excluded")
        if excluded_native != excluded_java:
            raise AssertionError({"label": "excluded", "native": excluded_native,
                                  "java": excluded_java})
        if [row[1] for row in excluded_native] != [1, 4]:
            raise AssertionError({"label": "excluded", "rows": excluded_native})
        if excluded_manifest["excluded_intervals"] != 1:
            raise AssertionError({"label": "excluded", "manifest": excluded_manifest})

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "default_filter_resolution_exact": True,
            "secondary_supplementary_qcfail_retained_by_default": True,
            "wellformed_cigar_n_filter_exact": True,
            "disable_tool_defaults_exact": True,
            "explicit_filter_reenable_exact": True,
            "exclude_interval_locus_boundary_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
