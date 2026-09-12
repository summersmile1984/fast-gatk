#!/usr/bin/env python3
"""Strict pinned-GATK oracle for the **ORDER** of the GenotypeGVCFs output header.

The gap this gate pins
----------------------
Every other GenotypeGVCFs gate compares header *content*: the ``##FILTER``
group in order, every other content line as a multiset.  None of them pins the
position of a header line, so a whole divergence class was invisible: native's
GATK-compatibility header preserves the *input's* line order plus its own INFO
rank list, while GATK's output is written by htsjdk in **sorted** order.  The
consequences measured on these fixtures are the ``##contig`` position (native
1, GATK 27), the position of ``INFO/AD`` (native last, GATK right after ``AC``)
and the order of the two ``##INFO=<ID=DP,...>`` lines.

htsjdk's comparator (read out of the pinned jar, not inferred)
--------------------------------------------------------------
``javap -p -c`` on ``third_party/gatk-package/gatk-4.6.2.0/
gatk-package-4.6.2.0-local.jar`` gives, verbatim:

    // htsjdk.variant.vcf.VCFHeaderLine  (every header line type except contigs)
    public int compareTo(java.lang.Object);
       0: aload_0
       1: invokevirtual #115  // Method toString:()Ljava/lang/String;
       4: aload_1
       5: invokevirtual #116  // Method java/lang/Object.toString:()Ljava/lang/String;
       8: invokevirtual #119  // Method java/lang/String.compareTo:(Ljava/lang/String;)I
      11: ireturn

    // htsjdk.variant.vcf.VCFContigHeaderLine  (overrides it)
    public int compareTo(java.lang.Object);
       0: aload_1
       1: instanceof    #2    // class htsjdk/variant/vcf/VCFContigHeaderLine
       4: ifeq          22
       7: aload_0
       8: getfield      #39   // Field contigIndex:Ljava/lang/Integer;
      11: aload_1
      12: checkcast     #2
      15: getfield      #39   // Field contigIndex:Ljava/lang/Integer;
      18: invokevirtual #130  // Method java/lang/Integer.compareTo:(Ljava/lang/Integer;)I
      21: ireturn
      22: aload_0
      23: aload_1
      24: invokespecial #132  // Method htsjdk/variant/vcf/VCFSimpleHeaderLine.compareTo:(Ljava/lang/Object;)I
      27: ireturn

    // htsjdk.variant.vcf.VCFHeader.getMetaDataInSortedOrder()
    public java.util.Set<VCFHeaderLine> getMetaDataInSortedOrder();
       0: aload_0
       1: new           #461  // class java/util/TreeSet
       4: dup
       5: aload_0
       6: getfield      #103  // Field mMetaData:Set;
       9: invokespecial #462  // Method java/util/TreeSet."<init>":(Collection;)V
      12: invokevirtual #458  // Method makeGetMetaDataSet:(Set;)Set
      15: areturn

    // htsjdk.variant.variantcontext.writer.VCFWriter
    //   .writeHeader(VCFHeader, Writer, String versionString, String streamName)
    //   writes "##fileformat=<version>\n", then iterates
    //   header.getMetaDataInSortedOrder(), SKIPS any line whose key satisfies
    //   VCFHeaderVersion.isFormatString(key), and writes "##" + line + "\n" for
    //   each; then "#" + HEADER_FIELDS + samples + "\n".

So the comparator is, precisely:

    1. the ``##fileformat`` line is written first, unconditionally, out of the
       header version -- it is not part of the sort;
    2. every other metadata line is ordered by ``toString()``, i.e. the **full
       line text without the leading ``##``** (``VCFHeaderLine.toString()`` is
       ``key + "=" + value``), compared as Java ``String.compareTo``
       (UTF-16 code unit order).  Because every line starts with its key and
       ``"="``, this is "by key, then by the rest of the line";
    3. ``VCFContigHeaderLine`` overrides ``compareTo``: two contig lines compare
       by ``contigIndex`` (the index the line was given while the *input*
       header was parsed, i.e. the input's declaration order -- measured: an
       input declaring ``chr1, chr10, chr2`` against a reference dictionary in
       the order ``chr2, chr10, chr1`` comes back ``chr1, chr10, chr2``, and
       vice versa); a contig line against a non-contig line falls back to
       ``toString().compareTo(...)``, i.e. plain text.  Contig lines therefore
       keep the input's order but move as a block to the position where
       ``##contig=...`` sorts among the rest: after every ``##INFO=`` and
       before ``##source=``;
    4. the sort is a ``TreeSet``, so two lines with equal ``toString()`` collapse
       to one -- but the header *set* de-duplicates by ``equals()``
       (class + key + value), so this only matters for two different line types
       spelling the same text.

Where the producer's identity lines land
----------------------------------------
``##GATKCommandLine`` and ``##source`` are ordinary header lines and sort like
any other.  Measured on these fixtures GATK writes ``##GATKCommandLine=...``
between ``##FORMAT=`` and ``##INFO=`` (``"GATKCommandLine" > "FORMAT"`` and
``< "INFO"``) and ``##source=GenotypeGVCFs`` last, after the contigs.  Both are
**producer identity** -- they name the program, its command line and its
wall-clock date -- so they are normalised away before the comparison, exactly as
in the sibling gates; everything else (``##fileformat``, ``##ALT``,
``##FILTER``, ``##FORMAT``, ``##INFO``, ``##contig``, ``#CHROM``) is a statement
about the file's contents and must match byte for byte **and position for
position**.

What this gate asserts
----------------------
For each fixture:
1. GATK's own output header still satisfies the comparator above (the gate's
   model of htsjdk is re-derived from GATK's bytes, so the gate cannot pass
   because the model drifted);
2. GATK's and native's identity-stripped header lines are byte-identical
   **position by position**, including ``##fileformat`` first and ``#CHROM``
   last;
3. the contig lines appear in the input header's declaration order on both
   sides (the comparator's item 3, asserted rather than assumed).
"""

from __future__ import annotations

import argparse
import bisect
import json
import os
import pathlib
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import oracle_guard  # noqa: E402

# --------------------------------------------------------------------------
# Fixtures.  The declarations are byte-identical across the three inputs; only
# their ORDER differs, so the fixtures isolate the ordering rule from the
# content rule (which the sibling gates already pin).
# --------------------------------------------------------------------------
INFO_DP = "##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>"
INFO_AD = "##INFO=<ID=AD,Number=R,Type=Integer,Description=Allele depths>"
ALT_NON_REF = ("##ALT=<ID=NON_REF,Description="
               "Represents any possible alternate allele>")
FORMAT_LINES = [
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>",
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>",
    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>",
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>",
]
CONTIG = "##contig=<ID=chr1,length=100>"
CHROM = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTAR"
FILEFORMAT = "##fileformat=VCFv4.2"

# The registered fixture shape (verify_genotype_gvcf_spandel_gatk_oracle.py
# HEADER): contig first, INFO before FORMAT -- clearly NOT htsjdk's order.
CONTIG_FIRST_BODY = [CONTIG, ALT_NON_REF, INFO_DP, INFO_AD, *FORMAT_LINES]
# Control: already in htsjdk's text order (ALT < FILTER < FORMAT < INFO <
# contig), i.e. the shape a GATK-written gVCF's header has.
SORTED_BODY = [ALT_NON_REF, *FORMAT_LINES, INFO_AD, INFO_DP, CONTIG]
# Third order: the reverse of the first, same content.
REVERSED_BODY = list(reversed(CONTIG_FIRST_BODY))

# verify_genotype_gvcf.py:378-381, byte for byte.
STAR_RECORD = (
    "chr1\t2\t.\tA\t*,G,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:12,8,0,0:0,0,100,100,100,100,100,100,100,100\n"
)

CASES = [
    {
        "case": "header-order-star-contig-first",
        "why": "the registered fixture's header shape (##contig first, ##INFO "
               "before ##FORMAT): native preserves the input order, GATK emits "
               "htsjdk's sorted order",
        "body": CONTIG_FIRST_BODY,
        "args": [],
        "gatk_contig_order": ["chr1"],
    },
    {
        "case": "header-order-star-input-already-sorted",
        "why": "control: the input header is already in htsjdk's text order, so "
               "only the positions native moves (##INFO pulled to the end, the "
               "##INFO/DP pair, ##contig) can still diverge",
        "body": SORTED_BODY,
        "args": [],
        "gatk_contig_order": ["chr1"],
    },
    {
        "case": "header-order-star-input-reversed",
        "why": "third distinct input order (the reverse of the first), same "
               "declarations: proves the rule is a sort of the header's own "
               "lines and not a permutation native could inherit",
        "body": REVERSED_BODY,
        "args": [],
        "gatk_contig_order": ["chr1"],
    },
    {
        "case": "header-order-star-dense",
        "why": "the same ordering contract in the other traversal "
               "(--include-non-variant-sites), which reaches the aggregate "
               "writer rather than the streaming one",
        "body": CONTIG_FIRST_BODY,
        "args": ["--include-non-variant-sites"],
        "gatk_contig_order": ["chr1"],
    },
]

# Producer identity: the program that wrote the file, its command line and its
# wall-clock date.  Native is a different program and has no
# --add-output-vcf-command-line equivalent on this surface, so these lines are
# normalised away before the comparison (they are a *record of the producer*,
# not a declaration about the file's contents).
IDENTITY_KEYS = frozenset({"GATKCommandLine", "source", "fileDate"})


def is_identity(line: str) -> bool:
    if not line.startswith("##"):
        return False
    key = line[2:].split("=", 1)[0]
    if key in IDENTITY_KEYS:
        return True
    # Belt and braces for a key that carries the producer's own run metadata
    # under another name (e.g. a future ##fileDate spelling inside a compound
    # line): a Date= or CommandLine= field is identity wherever it appears.
    return 'CommandLine="' in line or 'Date="' in line


def key_of(line: str) -> str:
    return line[2:].split("=", 1)[0] if line.startswith("##") else ""


def htsjdk_order(lines: list[str]) -> list[str]:
    """The order htsjdk's VCFWriter.writeHeader() would emit ``lines`` in.

    Implements exactly the comparator documented in this file's docstring: the
    ``##fileformat`` line first, then everything else by full line text, except
    that contiguous runs of ``##contig=`` lines keep their input order while
    moving to the text position of the ``##contig=`` prefix, then ``#CHROM``.
    """
    fileformat = [line for line in lines if key_of(line) == "fileformat"]
    chrom = [line for line in lines if line.startswith("#CHROM")]
    rest = [line for line in lines
            if key_of(line) != "fileformat" and not line.startswith("#CHROM")]
    contigs = [line for line in rest if line.startswith("##contig=")]
    others = sorted(line for line in rest if not line.startswith("##contig="))
    # "##contig=" is a strict prefix of every contig line, and no non-contig
    # line shares that prefix, so lowering against the prefix gives exactly the
    # position where the contig block sorts.
    at = bisect.bisect_left(others, "##contig=")
    return fileformat + others[:at] + contigs + others[at:] + chrom


def header_lines(path: pathlib.Path) -> list[str]:
    opener = None
    if path.name.endswith(".gz"):
        import gzip
        opener = gzip.open
    else:
        opener = open
    with opener(path, "rt", encoding="utf-8") as stream:
        return [line.rstrip("\n") for line in stream if line.startswith("#")]


def comparable(lines: list[str]) -> list[str]:
    return [line for line in lines if not is_identity(line)]


def position_diff(gatk: list[str], native: list[str], limit: int = 60) -> str:
    rows = []
    for index in range(max(len(gatk), len(native))):
        left = gatk[index] if index < len(gatk) else "<no line>"
        right = native[index] if index < len(native) else "<no line>"
        if left != right:
            rows.append(f"      [{index:2d}] GATK   {left}\n"
                        f"           NATIVE {right}")
    if len(rows) > limit:
        rows = rows[:limit] + [f"      ... {len(rows) - limit} more"]
    return "\n".join(rows)


def contig_names(lines: list[str]) -> list[str]:
    names = []
    for line in lines:
        if line.startswith("##contig=<ID="):
            names.append(line[len("##contig=<ID="):].split(",", 1)[0].split(">", 1)[0])
    return names


def write_reference(work: pathlib.Path) -> pathlib.Path:
    """A 100 bp chr1 whose bases match every fixture REF allele."""
    reference = work / "reference.fa"
    reference.write_text(">chr1\n" + "A" * 100 + "\n", encoding="utf-8")
    reference.with_name(reference.name + ".fai").write_text(
        "chr1\t100\t6\t100\t101\n", encoding="utf-8")
    reference.with_suffix(".dict").write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n", encoding="utf-8")
    return reference


def invoke(command: list[str], label: str, timeout: int) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True,
                            check=False, timeout=timeout)
    if result.returncode != 0:
        print(f"# {label} exited {result.returncode}")
        print(result.stderr[-4000:])
    return result


def run_case(case: dict, work: pathlib.Path, reference: pathlib.Path,
             java: pathlib.Path, jar: pathlib.Path, native: pathlib.Path,
             timeout: int) -> dict:
    source = work / f"{case['case']}.g.vcf"
    source.write_text(
        FILEFORMAT + "\n" + "\n".join(case["body"]) + "\n" + CHROM + "\n"
        + STAR_RECORD, encoding="utf-8")
    index_result = invoke([str(java), "-Xmx1g", "-jar", str(jar),
                           "IndexFeatureFile", "-I", str(source)],
                          f"GATK IndexFeatureFile [{case['case']}]", timeout)
    gatk_out = work / f"{case['case']}-gatk.vcf"
    native_out = work / f"{case['case']}-native.vcf"
    common = ["-R", str(reference), "-V", str(source), *case["args"]]
    gatk_result = invoke([str(java), "-Xmx1g", "-jar", str(jar),
                          "GenotypeGVCFs", *common, "-O", str(gatk_out),
                          "--create-output-variant-index", "false"],
                         f"GATK GenotypeGVCFs [{case['case']}]", timeout)
    native_result = invoke([str(native), *common,
                            "--gatk-compatible-annotations",
                            "-O", str(native_out)],
                           f"native GenotypeGVCFs [{case['case']}]", timeout)

    gatk_headers = header_lines(gatk_out) if gatk_out.exists() else []
    native_headers = header_lines(native_out) if native_out.exists() else []
    gatk_compared = comparable(gatk_headers)
    native_compared = comparable(native_headers)
    predicted = htsjdk_order(gatk_headers) if gatk_headers else []

    result = {
        "case": case["case"],
        "why": case["why"],
        "args": common,
        "gatk_exit": gatk_result.returncode,
        "native_exit": native_result.returncode,
        "gatk_header_lines": len(gatk_headers),
        "native_header_lines": len(native_headers),
        "gatk_compared_lines": len(gatk_compared),
        "native_compared_lines": len(native_compared),
        "gatk_header": gatk_headers,
        "native_header": native_headers,
        "comparator_reproduces_gatk": predicted == gatk_headers,
        "gatk_contig_order": contig_names(gatk_headers),
        "native_contig_order": contig_names(native_headers),
        "mismatched_positions": sum(
            1 for index in range(max(len(gatk_compared), len(native_compared)))
            if (gatk_compared[index] if index < len(gatk_compared) else None)
            != (native_compared[index] if index < len(native_compared) else None)),
        "violations": [],
    }
    if index_result.returncode != 0:
        result["violations"].append(
            f"GATK IndexFeatureFile exited {index_result.returncode}")
    if gatk_result.returncode != 0:
        result["violations"].append(f"GATK exited {gatk_result.returncode}")
    if native_result.returncode != 0:
        result["violations"].append(f"native exited {native_result.returncode}")
    # The gate's model of htsjdk must still describe pinned GATK's own output;
    # otherwise the comparison below could pass because the model drifted.
    if not result["comparator_reproduces_gatk"]:
        result["violations"].append(
            "GATK moved away from the pinned htsjdk comparator: its own header "
            "is no longer '##fileformat' first + the rest sorted by line text "
            "(contigs by input index) + '#CHROM' last")
    # Producer-identity lines are normalised away, so the comparison must still
    # be a comparison of a real header.
    if result["gatk_compared_lines"] < 20:
        result["violations"].append(
            "the identity-stripped GATK header is too short to be a header: "
            f"{result['gatk_compared_lines']} lines")
    if result["gatk_contig_order"] != case["gatk_contig_order"]:
        result["violations"].append(
            "GATK moved away from the measured contig order: "
            f"{result['gatk_contig_order']} expected {case['gatk_contig_order']}")
    if result["native_contig_order"] != result["gatk_contig_order"]:
        result["violations"].append(
            "contig lines are not in the same order: "
            f"GATK={result['gatk_contig_order']} "
            f"NATIVE={result['native_contig_order']}")
    if len(gatk_compared) != len(native_compared):
        result["violations"].append(
            "identity-stripped header line counts differ: "
            f"GATK={len(gatk_compared)} NATIVE={len(native_compared)}")
    if gatk_compared != native_compared:
        result["violations"].append(
            "header lines are not byte-identical in the same positions "
            f"({result['mismatched_positions']} position(s) differ):\n"
            + position_diff(gatk_compared, native_compared))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict oracle for the ORDER of the GenotypeGVCFs output "
                    "header against pinned GATK 4.6.2.0 (htsjdk's comparator).")
    parser.add_argument(
        "--native", default=os.environ.get("FASTGATK_GENOTYPE_BINARY"),
        help="native GenotypeGVCFs binary (default: $FASTGATK_GENOTYPE_BINARY "
             "or fastgatk-native/build/fastgatk-genotype-gvcf)")
    parser.add_argument("--case", action="append", default=None,
                        help="run only the named case(s)")
    parser.add_argument("--timeout", type=int, default=300,
                        help="per-process timeout in seconds (default 300)")
    arguments = parser.parse_args()

    root = pathlib.Path(__file__).resolve().parents[2]
    native = pathlib.Path(arguments.native) if arguments.native else (
        pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD",
                                    root / "fastgatk-native" / "build"))
        / "fastgatk-genotype-gvcf")
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = pathlib.Path(os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java")))

    assets = [native, java, jar]
    if not all(path.is_file() for path in assets):
        oracle_guard.oracle_not_verified(
            "verify_genotype_gvcf_header_order_gatk_oracle.py", java, jar)
        missing = sorted(str(path) for path in assets if not path.is_file())
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing oracle assets: {missing}")
        print(json.dumps({"status": "skipped",
                          "reason": "missing native binary or pinned GATK assets",
                          "missing": missing}, sort_keys=True))
        return 0

    selected = [case for case in CASES
                if arguments.case is None or case["case"] in arguments.case]
    if not selected:
        raise SystemExit(f"no such case: {arguments.case}")

    results: list[dict] = []
    with tempfile.TemporaryDirectory(
            prefix="fastgatk-genotype-header-order-oracle-") as directory:
        work = pathlib.Path(directory)
        reference = write_reference(work)
        for case in selected:
            results.append(run_case(case, work, reference, java, jar, native,
                                    arguments.timeout))

    violations: list[str] = []
    for result in results:
        violations.extend(f"[{result['case']}] {item}"
                          for item in result["violations"])

    print(f"# {pathlib.Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# GATK rule under test: htsjdk's VCFWriter.writeHeader() writes "
          "##fileformat first, then every metadata line in VCFHeaderLine"
          ".compareTo order -- the FULL LINE TEXT (VCFContigHeaderLine instead "
          "compares by the input parse order of the contigs, and against a "
          "non-contig line by text) -- then #CHROM.  ##GATKCommandLine and "
          "##source are producer identity and are normalised away.")
    print(f"# cases={[case['case'] for case in selected]}")
    for result in results:
        print(f"[{result['case']}]")
        print(f"    why: {result['why']}")
        print(f"    args: {' '.join(result['args']) if result['args'] else '(none)'}")
        print(f"    gatk_exit={result['gatk_exit']} native_exit={result['native_exit']}")
        print(f"    header lines: GATK={result['gatk_header_lines']} "
              f"NATIVE={result['native_header_lines']}; compared "
              f"(identity stripped): GATK={result['gatk_compared_lines']} "
              f"NATIVE={result['native_compared_lines']}")
        print(f"    comparator reproduces GATK's own order: "
              f"{result['comparator_reproduces_gatk']}")
        print(f"    contig order: GATK={result['gatk_contig_order']} "
              f"NATIVE={result['native_contig_order']}")
        print(f"    mismatched positions: {result['mismatched_positions']}")
        if result["mismatched_positions"]:
            print("    position-by-position diff:")
            print(position_diff(comparable(result["gatk_header"]),
                                comparable(result["native_header"])))
        for violation in result["violations"]:
            print(f"    VIOLATION: {violation}")

    payload = {
        "status": "pass" if not violations else "divergence",
        "gatk_version": "4.6.2.0",
        "cases": results,
        "violations": violations,
        "header_compared": "every header line of both outputs, in order, after "
                           "normalising away the producer-identity lines "
                           "(##GATKCommandLine, ##source, ##fileDate and any "
                           "CommandLine=/Date= field); includes ##fileformat "
                           "first and #CHROM last",
        "htsjdk_comparator": "VCFHeaderLine.compareTo = toString().compareTo(); "
                             "VCFContigHeaderLine.compareTo = contigIndex order "
                             "against another contig, text otherwise; "
                             "VCFHeader.getMetaDataInSortedOrder() = TreeSet; "
                             "VCFWriter.writeHeader writes ##fileformat first "
                             "and skips isFormatString keys",
    }
    print(json.dumps(payload, sort_keys=True))
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
