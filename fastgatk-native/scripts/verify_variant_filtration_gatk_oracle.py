#!/usr/bin/env python3
"""Focused GATK 4.6.2.0 VariantFiltration semantic oracle.

This oracle is intentionally separate from ``verify_variant_filtration.py``.
It keeps a small, release-pinned fixture and compares the records emitted by
the native implementation with GATK's Java writer for the common production
slice: site QUAL/INFO predicates, genotype FORMAT predicates and FILTER/FT
materialization.  A second invocation covers ``--invert-filter-expression``.
Headers, compression and index bytes are not compared; all semantic record
columns (including QUAL, INFO, FORMAT and genotype values) are compared after
normalizing only FILTER label order and INFO key order.

JEXL outside this explicitly tested subset (arbitrary annotation objects,
annotation engines, cloud codecs and full Java expression semantics) remains a
documented fallback rather than being silently treated as supported.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import subprocess
import tempfile
from decimal import Decimal, InvalidOperation
from pathlib import Path


VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##INFO=<ID=QD,Number=1,Type=Float,Description=Quality by depth>
##INFO=<ID=AS_SCORE,Number=A,Type=Integer,Description=Per-ALT score>
##FILTER=<ID=LowExisting,Description=An input filter retained by the oracle fixture>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=FT,Number=1,Type=String,Description=Genotype filter>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t10\trsLow\tA\tG\t25\tPASS\tDP=5;QD=1.0;AS_SCORE=3\tGT:GQ:DP:AD\t0/1:10:5:3,2\t1/1:50:20:0,20
chr1\t20\trsGood\tC\tT\t60\tPASS\tDP=20;QD=3.0;AS_SCORE=8\tGT:GQ:DP:AD\t0/0:50:20:20,0\t0/1:5:12:8,4
chr1\t30\trsExisting\tG\tA\t60\tLowExisting\tDP=30;QD=2.5;AS_SCORE=9\tGT:GQ:DP:AD:FT\t0/1:40:30:15,15:Prior\t./.:.:.:.:PASS
"""


def invoke(command: list[str], label: str) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"{label} failed ({result.returncode}): {result.stderr[-4000:]}")
    return result


def canonical_numeric(value: str) -> str:
    """Canonicalize Java/HTSlib textual spellings without using binary float."""
    try:
        number = Decimal(value)
    except InvalidOperation:
        return value
    if not number.is_finite():
        return value
    rendered = format(number.normalize(), "f")
    if "." in rendered:
        rendered = rendered.rstrip("0").rstrip(".")
    return rendered or "0"


def read_vcf(path: Path) -> list[tuple[object, ...]]:
    opener = gzip.open if path.name.endswith(".gz") else open
    rows: list[tuple[object, ...]] = []
    info_types: dict[str, str] = {}
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("##INFO=<"):
                match = re.search(r"ID=([^,>]+).*Type=(Integer|Float|Flag)(?:,|>)", line)
                if match:
                    info_types[match.group(1)] = match.group(2)
                continue
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                raise AssertionError(f"unexpected VCF row in {path}: {line!r}")
            filters = tuple(sorted(name for name in fields[6].split(";")
                                   if name and name != "PASS"))
            info_items = []
            if fields[7] not in ("", "."):
                for item in fields[7].split(";"):
                    key, separator, value = item.partition("=")
                    if separator and info_types.get(key) in ("Integer", "Float"):
                        value = ",".join(canonical_numeric(token) for token in value.split(","))
                    info_items.append((key, value if separator else None))
            format_keys = fields[8].split(":")
            samples = []
            for sample in fields[9:]:
                values = sample.split(":")
                # HTSJDK omits trailing missing FORMAT values, while HTSlib
                # commonly writes explicit dots.  They denote the same VCF
                # value, so pad only omitted trailing values before comparing.
                if len(values) > len(format_keys):
                    raise AssertionError(f"FORMAT/sample width mismatch in {path}: {line!r}")
                values.extend("." for _ in range(len(format_keys) - len(values)))
                # HTSJDK may emit FORMAT keys in a different canonical order
                # from HTSlib. Compare key/value semantics, not that writer
                # implementation detail; every value is retained exactly.
                samples.append(tuple(sorted(zip(format_keys, values))))
            rows.append((fields[0], fields[1], fields[2], fields[3], fields[4],
                         fields[5], filters, tuple(sorted(info_items)),
                         tuple(sorted(format_keys)), *samples))
    return rows


def header_line(path: Path, prefix: str) -> str | None:
    """Return one metadata line without writer-specific quoting differences."""
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith(prefix):
                return line.rstrip("\n")
    return None


def run_native(binary: Path, source: Path, output: Path, args: list[str]) -> dict:
    result = invoke([str(binary), "-V", str(source), "-O", str(output),
                     "--create-output-variant-index=false", *args], "native VariantFiltration")
    try:
        return json.loads(result.stdout.splitlines()[-1])
    except (IndexError, json.JSONDecodeError) as error:
        raise AssertionError(f"native did not emit a JSON summary: {result.stdout!r}") from error


def run_gatk(java: Path, jar: Path, source: Path, output: Path, args: list[str]) -> None:
    invoke([str(java), "-jar", str(jar), "VariantFiltration", "-V", str(source),
            "-O", str(output), *args], "GATK VariantFiltration")


def assert_exact(native: Path, gatk: Path, label: str) -> int:
    native_rows = read_vcf(native)
    gatk_rows = read_vcf(gatk)
    if native_rows != gatk_rows:
        raise AssertionError(
            f"{label} native/GATK semantic rows differ:\n"
            f"native={native_rows!r}\nGATK={gatk_rows!r}")
    return len(native_rows)


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_VARIANT_FILTRATION_BINARY",
        root / "fastgatk-native/build/fastgatk-variant-filtration"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    required = (binary, java, jar)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("VariantFiltration Java oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-filtration-gatk-oracle-") as directory:
        work = Path(directory)
        source = work / "input.vcf"
        source.write_text(VCF, encoding="utf-8")

        common = [
            "--filter-expression", "QUAL < 30", "--filter-name", "LowQual",
            "--filter-expression", 'vc.getAttribute("DP") < 10', "--filter-name", "LowDP",
            "--genotype-filter-expression", "GQ < 20", "--genotype-filter-name", "LowGQ",
        ]
        gatk = work / "gatk.vcf"
        native = work / "native.vcf"
        run_gatk(java, jar, source, gatk, common)
        native_summary = run_native(binary, source, native, common)
        record_count = assert_exact(native, gatk, "site+genotype predicates")
        if record_count != 3 or native_summary.get("input_records") != 3:
            raise AssertionError(f"unexpected VariantFiltration summary: {native_summary}")

        invert_args = ["--filter-expression", "QUAL < 30", "--filter-name", "NotLowQual",
                       "--invert-filter-expression"]
        gatk_invert = work / "gatk-invert.vcf"
        native_invert = work / "native-invert.vcf"
        run_gatk(java, jar, source, gatk_invert, invert_args)
        invert_summary = run_native(binary, source, native_invert, invert_args)
        invert_count = assert_exact(native_invert, gatk_invert, "inverted site predicate")
        if invert_count != 3 or invert_summary.get("input_records") != 3:
            raise AssertionError(f"unexpected inverted summary: {invert_summary}")

        # The first record exercises both site predicates and the S1
        # genotype predicates; this guards against a summary-only smoke test.
        rows = read_vcf(native)
        first = rows[0]
        if first[6] != ("LowDP", "LowQual"):
            raise AssertionError(f"unexpected native FILTER set: {first[6]!r}")
        if "FT" not in str(first[8]):
            raise AssertionError(f"native genotype FT column missing: {first!r}")
        if "LowGQ" not in str(first[9]):
            raise AssertionError(f"native genotype FT missing expected labels: {first!r}")
        invert_rows = read_vcf(native_invert)
        if invert_rows[0][6] != () or invert_rows[1][6] != ("NotLowQual",):
            raise AssertionError(f"unexpected invert FILTER semantics: {invert_rows!r}")

        # GATK's invalidate-previous-filters resets only site-level FILTER;
        # genotype FT history remains available for downstream chained passes.
        invalidate_args = [
            "--filter-expression", "QUAL < 30", "--filter-name", "NewSite",
            "--genotype-filter-expression", "GQ < 20", "--genotype-filter-name", "NewGT",
            "--invalidate-previous-filters", "true",
        ]
        gatk_invalidate = work / "gatk-invalidate.vcf"
        native_invalidate = work / "native-invalidate.vcf"
        invalidate_manifest = work / "native-invalidate.manifest.json"
        run_gatk(java, jar, source, gatk_invalidate, invalidate_args)
        invalidate_summary = run_native(
            binary, source, native_invalidate,
            [*invalidate_args, "--output-manifest", str(invalidate_manifest)])
        invalidate_count = assert_exact(
            native_invalidate, gatk_invalidate, "invalidated previous site filters")
        if invalidate_count != 3 or invalidate_summary.get("input_records") != 3:
            raise AssertionError(f"unexpected invalidate summary: {invalidate_summary}")
        invalidated_rows = read_vcf(native_invalidate)
        if invalidated_rows[0][6] != ("NewSite",):
            raise AssertionError(f"site FILTER history was not reset: {invalidated_rows!r}")
        if "NewGT" not in str(invalidated_rows[0][9]) or "Prior" not in str(invalidated_rows[2][9]):
            raise AssertionError(f"genotype FT history was unexpectedly removed: {invalidated_rows!r}")
        invalidate_metadata = json.loads(invalidate_manifest.read_text(encoding="utf-8"))
        if invalidate_metadata["compatibility"]["invalidate_previous_filters"] is not True:
            raise AssertionError(f"manifest missed invalidate_previous_filters: {invalidate_metadata}")

        # VariantFiltration's mask is a FeatureInput and GATK requires a
        # random-access index for it.  Keep the mask fixture indexed by the
        # pinned Java IndexFeatureFile command, then compare both normal
        # overlap semantics and GATK's reverse --filter-not-in-mask mode.
        # The input's pre-existing LowExisting filter makes this a chained
        # filtering check as well: a mask must append, never erase, history.
        mask = work / "mask.vcf"
        mask.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=100>\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            "chr1\t20\tmask\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/1\n",
            encoding="utf-8",
        )
        invoke([str(java), "-jar", str(jar), "IndexFeatureFile", "-I", str(mask)],
               "GATK IndexFeatureFile for VariantFiltration mask")

        mask_args = ["--mask", str(mask)]
        gatk_mask = work / "gatk-mask.vcf"
        native_mask = work / "native-mask.vcf"
        run_gatk(java, jar, source, gatk_mask, mask_args)
        native_mask_summary = run_native(binary, source, native_mask, mask_args)
        mask_count = assert_exact(native_mask, gatk_mask, "mask overlap")
        if mask_count != 3 or native_mask_summary.get("masked_records") != 1:
            raise AssertionError(f"unexpected mask summary: {native_mask_summary}")
        mask_rows = read_vcf(native_mask)
        if mask_rows[0][6] != () or mask_rows[1][6] != ("Mask",) or mask_rows[2][6] != ("LowExisting",):
            raise AssertionError(f"unexpected mask FILTER semantics: {mask_rows!r}")
        mask_metadata = json.loads(Path(str(native_mask) + ".manifest.json").read_text(encoding="utf-8"))
        if mask_metadata["compatibility"]["mask_filter"] is not True or \
                mask_metadata["compatibility"]["filter_not_in_mask"] is not False:
            raise AssertionError(f"mask compatibility metadata missing: {mask_metadata}")

        reverse_args = ["--mask", str(mask), "--mask-name", "NotMask", "--filter-not-in-mask"]
        gatk_reverse = work / "gatk-mask-reverse.vcf"
        native_reverse = work / "native-mask-reverse.vcf"
        run_gatk(java, jar, source, gatk_reverse, reverse_args)
        native_reverse_summary = run_native(binary, source, native_reverse, reverse_args)
        reverse_count = assert_exact(native_reverse, gatk_reverse, "reverse mask")
        if reverse_count != 3 or native_reverse_summary.get("masked_records") != 2 or \
                native_reverse_summary.get("filter_not_in_mask") is not True:
            raise AssertionError(f"unexpected reverse mask summary: {native_reverse_summary}")
        reverse_rows = read_vcf(native_reverse)
        if reverse_rows[0][6] != ("NotMask",) or reverse_rows[1][6] != () or \
                reverse_rows[2][6] != ("LowExisting", "NotMask"):
            raise AssertionError(f"unexpected reverse mask FILTER semantics: {reverse_rows!r}")

        # A custom mask description is a header-only option, but it must be
        # accepted by both implementations and retained in the FILTER
        # declaration.  Quote normalization is intentionally ignored here:
        # htsjdk quotes descriptions while HTSlib leaves plain text unquoted.
        described_args = ["--mask", str(mask), "--mask-name", "CustomMask",
                          "--mask-description", "my mask 100%"]
        gatk_described = work / "gatk-mask-described.vcf"
        native_described = work / "native-mask-described.vcf"
        run_gatk(java, jar, source, gatk_described, described_args)
        native_described_summary = run_native(binary, source, native_described, described_args)
        assert_exact(native_described, gatk_described, "custom mask description")
        native_mask_header = header_line(native_described, "##FILTER=<ID=CustomMask,")
        gatk_mask_header = header_line(gatk_described, "##FILTER=<ID=CustomMask,")
        if native_mask_header is None or gatk_mask_header is None or \
                "Description=my mask 100%" not in native_mask_header or \
                'Description="my mask 100%"' not in gatk_mask_header:
            raise AssertionError(
                f"custom mask description missing: native={native_mask_header!r}, gatk={gatk_mask_header!r}")
        described_metadata = json.loads(
            Path(str(native_described) + ".manifest.json").read_text(encoding="utf-8"))
        if described_metadata["compatibility"]["mask_description"] != "my mask 100%":
            raise AssertionError(f"mask description metadata missing: {described_metadata}")

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "site_genotype_records": record_count,
            "inverted_records": invert_count,
            "qual_info_format_semantics_exact": True,
            "filter_inversion_exact": True,
            "invalidate_previous_filters_exact": True,
            "mask_overlap_exact": True,
            "mask_reverse_exact": True,
            "mask_description_exact": True,
            "headers_and_index_bytes_compared": False,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
