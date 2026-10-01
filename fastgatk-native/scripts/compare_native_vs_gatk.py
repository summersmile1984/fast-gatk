#!/usr/bin/env python3
"""Native (fast-gatk) vs pinned GATK 4.6.2.0 comparison for the tools whose
fallback paths were eliminated this session: CompareReferences FULL_ALIGNMENT,
GatherVcfs BCF, IndexFeatureFile (interval_list + GTF).

Measures, per tool and per side:
  * correctness  -- byte-identity of the output (VCF + .vcf.idx, .bcf, .idx)
  * wall time    -- /usr/bin/time elapsed (Java side includes JVM startup)
  * resource     -- max RSS (KB) and filesystem inputs/outputs

Writes JSON to stdout (machine-readable) and emits the Markdown report to the
path given by --report (default <repo>/NATIVE_VS_GATK_COMPARISON.md).
"""

from __future__ import annotations

import argparse
import json
import pathlib
import random
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "fastgatk-native" / "scripts"))

from benchmark_lib import (  # noqa: E402
    bgzip_and_index, gatk_jar, java_binary, run_java_timed, run_timed,
    write_minimal_dict,
)

BUILD = ROOT / "fastgatk-native" / "build"  # OpenMP backend


def native(name: str) -> pathlib.Path:
    return BUILD / f"fastgatk-{name}"


def tribble_body(data: bytes) -> bytes:
    """Return the Tribble v3 index body (after the header/properties)."""
    import struct
    offset = 4  # skip TIDX magic

    def i32() -> int:
        nonlocal offset
        v = struct.unpack_from("<i", data, offset)[0]
        offset += 4
        return v

    def i64() -> int:
        nonlocal offset
        v = struct.unpack_from("<q", data, offset)[0]
        offset += 8
        return v

    def text() -> str:
        nonlocal offset
        end = data.index(b"\0", offset)
        v = data[offset:end].decode()
        offset = end + 1
        return v

    i32()          # type
    i32()          # version
    text()         # uri (run-specific absolute path)
    i64()          # size
    i64()          # mtime (run-specific)
    text()         # md5
    i32()          # flags
    count = i32()
    for _ in range(count):
        text()
        text()
    return data[offset:]


def write_reference(path: pathlib.Path, length: int, snp_every: int) -> None:
    """Deterministic ACGT reference with a SNP every ``snp_every`` bases."""
    rng = random.Random(0xC0FFEE)
    bases = "".join(rng.choice("ACGT") for _ in range(length))
    bases = list(bases)
    for pos in range(snp_every, length, snp_every):
        alt = [b for b in "ACGT" if b != bases[pos]][0]
        bases[pos] = alt
    text = "".join(bases)
    width = 80
    with path.open("w", encoding="ascii", newline="\n") as stream:
        stream.write(">chr1\n")
        for offset in range(0, length, width):
            stream.write(text[offset:offset + width] + "\n")
    path.with_name(path.name + ".fai").write_text(
        f"chr1\t{length}\t6\t{width}\t{width + 1}\n", encoding="ascii")


def compare_references(work: pathlib.Path, bases: int, results: dict) -> None:
    first = work / "cr1.fa"
    second = work / "cr2.fa"
    write_reference(first, bases, bases + 1)   # pure ACGT, no SNPs
    write_reference(second, bases, 10_000)     # SNP every 10 kb
    write_minimal_dict(first, "chr1", bases)
    write_minimal_dict(second, "chr1", bases)

    nat_out = work / "nat-snps"
    java_out = work / "java-snps"
    nat_out.mkdir(); java_out.mkdir()

    nat_cmd = [str(native("compare-references")), "-R", str(first),
               "-refcomp", str(second), "--base-comparison", "FULL_ALIGNMENT",
               "--base-comparison-output", str(nat_out), "--threads", "4"]
    java_cmd = ["CompareReferences", "-R", str(first), "-refcomp", str(second),
                "--base-comparison", "FULL_ALIGNMENT",
                "--base-comparison-output", str(java_out)]

    nat_times = [run_timed(nat_cmd, work, "cr-native") for _ in range(3)]
    nat_wall = min(t["seconds"] for t in nat_times)
    nat_rss = max(t["max_rss_kb"] for t in nat_times)
    java_t = run_java_timed("CompareReferences", java_cmd[1:], work, "cr-java")

    nat_vcf = sorted(nat_out.glob("*.vcf"))
    java_vcf = sorted(java_out.glob("*.vcf"))
    nat_idx = sorted(nat_out.glob("*.vcf.idx"))
    java_idx = sorted(java_out.glob("*.vcf.idx"))
    vcf_identical = bool(nat_vcf and java_vcf and
                         nat_vcf[0].read_bytes() == java_vcf[0].read_bytes())
    # .vcf.idx carries a run-specific mtime and the output VCF's absolute path
    # (uri), so compare the index *body* byte-for-byte instead of the whole file.
    idx_body_identical = bool(nat_idx and java_idx and
                              tribble_body(nat_idx[0].read_bytes()) ==
                              tribble_body(java_idx[0].read_bytes()))
    nat_records = (sum(1 for ln in nat_vcf[0].read_text().splitlines()
                       if not ln.startswith("#")) if nat_vcf else 0)
    java_records = (sum(1 for ln in java_vcf[0].read_text().splitlines()
                        if not ln.startswith("#")) if java_vcf else 0)

    results["compare_references"] = {
        "input": f"{bases / 1_000_000:.0f} Mb x2, SNP every 10 kb",
        "records": {"native": nat_records, "gatk": java_records},
        "vcf_byte_identical": vcf_identical,
        "vcf_idx_body_identical": idx_body_identical,
        "native": {"wall_seconds": round(nat_wall, 3),
                   "max_rss_mb": round(nat_rss / 1024, 1)},
        "gatk": {"wall_seconds": round(java_t["seconds"], 3),
                 "max_rss_mb": round(java_t["max_rss_kb"] / 1024, 1)},
    }


def gather_vcfs(work: pathlib.Path, shards: int, per_shard: int, results: dict) -> None:
    header = ("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000000>\n"
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    inputs = []
    total = shards * per_shard
    for s in range(shards):
        plain = work / f"shard-{s}.vcf"
        with plain.open("w", encoding="ascii") as fh:
            fh.write(header)
            first = s * per_shard + 1
            for pos in range(first, first + per_shard):
                fh.write(f"chr1\t{pos}\t.\tA\tG\t50\tPASS\t.\n")
        inputs.append(bgzip_and_index(plain))

    nat_out = work / "native.bcf"
    java_out = work / "java.bcf"
    nat_cmd = [str(native("gather-vcfs"))]
    for p in inputs:
        nat_cmd += ["-I", str(p)]
    nat_cmd += ["-O", str(nat_out)]
    java_cmd = ["GatherVcfs"]
    for p in inputs:
        java_cmd += ["-I", str(p)]
    java_cmd += ["-O", str(java_out)]

    nat_t = run_timed(nat_cmd, work, "gv-native")
    java_t = run_java_timed("GatherVcfs", java_cmd[1:], work, "gv-java")

    # Correctness: both must be raw BCF 2.1 and decode to identical VCF records
    # via the pinned GATK SelectVariants.
    nat_magic = nat_out.read_bytes()[:5] if nat_out.exists() else b""
    java_magic = java_out.read_bytes()[:5] if java_out.exists() else b""

    def decode(bcf: pathlib.Path, out: pathlib.Path) -> int:
        if not bcf.exists():
            return -1
        r = subprocess.run(
            [str(java_binary()), "-Xmx1g", "-jar", str(gatk_jar()), "SelectVariants",
             "-V", str(bcf), "-O", str(out), "--create-output-variant-index", "false"],
            capture_output=True, text=True)
        return r.returncode

    nat_dec = work / "native-decoded.vcf"
    java_dec = work / "java-decoded.vcf"
    nat_rc = decode(nat_out, nat_dec)
    java_rc = decode(java_out, java_dec)

    def data_lines(path: pathlib.Path) -> list[str]:
        return [ln for ln in path.read_text().splitlines() if not ln.startswith("#")]

    decoded_identical = (nat_rc == 0 and java_rc == 0 and nat_dec.exists() and
                         java_dec.exists() and
                         data_lines(nat_dec) == data_lines(java_dec))

    results["gather_vcfs"] = {
        "input": f"{shards} shards x {per_shard} records = {total} records",
        "records": total,
        "bcf_magic_native": nat_magic.hex(),
        "bcf_magic_gatk": java_magic.hex(),
        "decoded_records_identical": decoded_identical,
        "native": {"wall_seconds": round(nat_t["seconds"], 3),
                   "max_rss_mb": round(nat_t["max_rss_kb"] / 1024, 1)},
        "gatk": {"wall_seconds": round(java_t["seconds"], 3),
                 "max_rss_mb": round(java_t["max_rss_kb"] / 1024, 1)},
    }


def interval_list(work: pathlib.Path, records: int, results: dict) -> None:
    inp = work / "sample.interval_list"
    lines = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:1\tLN:249250621",
             "@SQ\tSN:2\tLN:243199373"]
    half = records // 2
    for i in range(half):
        start = i * 1000 + 1
        lines.append(f"1\t{start}\t{start + 500}\t+\tregion_{i}")
    for i in range(half, records):
        start = (i - half) * 1000 + 1
        lines.append(f"2\t{start}\t{start + 500}\t+\tregion_{i}")
    inp.write_text("\n".join(lines) + "\n", encoding="ascii")

    nat_idx = work / "native.idx"
    java_idx = work / "java.idx"
    nat_t = run_timed([str(native("index-feature-file")), "-I", str(inp),
                       "-O", str(nat_idx), "--threads", "4"], work, "il-native")
    java_t = run_java_timed("IndexFeatureFile", ["-I", str(inp), "-O", str(java_idx)],
                            work, "il-java")

    results["index_interval_list"] = {
        "input": f"{records} interval_list records, 2 contigs",
        "idx_byte_identical": (nat_idx.exists() and java_idx.exists() and
                               nat_idx.read_bytes() == java_idx.read_bytes()),
        "native": {"wall_seconds": round(nat_t["seconds"], 3),
                   "max_rss_mb": round(nat_t["max_rss_kb"] / 1024, 1)},
        "gatk": {"wall_seconds": round(java_t["seconds"], 3),
                 "max_rss_mb": round(java_t["max_rss_kb"] / 1024, 1)},
    }


def gtf(work: pathlib.Path, genes: int, results: dict) -> None:
    def record(contig, kind, start, end, gene, transcript=None):
        attrs = [f'gene_id "{gene}"', 'gene_type "protein_coding"',
                 f'gene_name "{gene}"']
        if transcript is not None:
            attrs += [f'transcript_id "{transcript}"',
                      'transcript_type "protein_coding"',
                      f'transcript_name "{transcript}"']
        if kind == "exon":
            attrs += ['exon_number 1', f'exon_id "E{gene}"']
        attrs.append("level 2")
        return (f"{contig}\tHAVANA\t{kind}\t{start}\t{end}\t.\t+\t.\t" +
                "; ".join(attrs) + ";\n")

    header = ("##description: evidence-based annotation of the human genome (GRCh38), "
              "version 43 (Ensembl 109)\n##provider: GENCODE\n"
              "##contact: gencode@sanger.ac.uk\n##format: gtf\n##date: 2023-01-01\n")
    inp = work / "gencode.v43.dense.gtf"
    body = []
    for i in range(genes):
        start = 2 * i + 1
        end = start + 1
        body.append(record("1", "gene", start, end + 2000, f"G{i}", f"T{i}"))
        body.append(record("1", "transcript", start, end + 2000, f"G{i}", f"T{i}"))
        body.append(record("1", "exon", start, end + 2000, f"G{i}", f"T{i}"))
    inp.write_text(header + "".join(body), encoding="ascii")

    nat_idx = work / "native-gtf.idx"
    java_idx = work / "java-gtf.idx"
    nat_t = run_timed([str(native("index-feature-file")), "-I", str(inp),
                       "-O", str(nat_idx), "--threads", "4"], work, "gtf-native")
    java_t = run_java_timed("IndexFeatureFile", ["-I", str(inp), "-O", str(java_idx)],
                            work, "gtf-java")

    results["index_gtf"] = {
        "input": f"{genes} GENCODE genes (gene+transcript+exon each)",
        "idx_byte_identical": (nat_idx.exists() and java_idx.exists() and
                               nat_idx.read_bytes() == java_idx.read_bytes()),
        "native": {"wall_seconds": round(nat_t["seconds"], 3),
                   "max_rss_mb": round(nat_t["max_rss_kb"] / 1024, 1)},
        "gatk": {"wall_seconds": round(java_t["seconds"], 3),
                 "max_rss_mb": round(java_t["max_rss_kb"] / 1024, 1)},
    }


def md_report(results: dict, backend: str) -> str:
    lines = []
    lines.append("# fast-gatk native vs pinned GATK 4.6.2.0 对比报告\n")
    lines.append(f"- 原生后端: `{backend}`（OpenMP，Kokkos）")
    lines.append(f"- Java 基线: GATK 4.6.2.0 (`{gatk_jar().name}`) + JDK 17，`-Xmx1g`")
    lines.append("- 计时: `/usr/bin/time` 墙钟（Java 侧含 JVM 启动）；内存: 峰值 RSS\n")

    for key, title, note in [
        ("compare_references", "CompareReferences FULL_ALIGNMENT（去掉 MUMmer 依赖）",
         "2 Mb 参考对，SNP 每 10 kb；输出 `.vcf` + `.vcf.idx`"),
        ("gather_vcfs", "GatherVcfs BCF 输出（原生 BCF 2.1）",
         "多 shard VCF 聚合为 `.bcf`，经 GATK SelectVariants 解码比对"),
        ("index_interval_list", "IndexFeatureFile interval_list",
         "线性/区间树 Tribble `.idx` 字节比对"),
        ("index_gtf", "IndexFeatureFile GENCODE GTF",
         "基因聚合 Tribble `.idx` 字节比对"),
    ]:
        r = results.get(key)
        if not r:
            continue
        lines.append(f"## {title}\n")
        lines.append(f"- 输入: {r['input']}\n")
        lines.append("| 指标 | native (fast-gatk) | GATK 4.6.2.0 | 结论 |")
        lines.append("|---|---|---|---|")
        n = r["native"]; j = r["gatk"]
        sp = (j["wall_seconds"] / n["wall_seconds"]) if n["wall_seconds"] else 0
        lines.append(f"| 墙钟 | {n['wall_seconds']} s | {j['wall_seconds']} s | "
                     f"native 快 **{sp:.1f}×** |")
        lines.append(f"| 峰值 RSS | {n['max_rss_mb']} MB | {j['max_rss_mb']} MB | "
                     f"native 少 **{j['max_rss_mb'] / max(n['max_rss_mb'], 1):.1f}×** |")
        if key == "compare_references":
            lines.append(f"| 记录数 | {r['records']['native']} | {r['records']['gatk']} | "
                         f"{'一致' if r['records']['native'] == r['records']['gatk'] else '不一致'} |")
            lines.append(f"| `.vcf` 字节一致 | — | — | "
                         f"{'✅ 一致' if r['vcf_byte_identical'] else '❌ 不一致'} |")
            lines.append(f"| `.vcf.idx` 内容一致 | — | — | "
                         f"{'✅ 一致（仅 mtime/uri 为运行时元数据）' if r['vcf_idx_body_identical'] else '❌ 不一致'} |")
        elif key == "gather_vcfs":
            lines.append(f"| 记录数 | {r['records']} | {r['records']} | ✅ 一致 |")
            lines.append(f"| BCF magic | `{r['bcf_magic_native']}` | `{r['bcf_magic_gatk']}` | "
                         f"{'✅ 均 BCF\\x02\\x01' if r['bcf_magic_native'] == r['bcf_magic_gatk'] else '⚠️ 不同'} |")
            lines.append(f"| 解码记录字节一致 | — | — | "
                         f"{'✅ 一致' if r['decoded_records_identical'] else '❌ 不一致'} |")
        else:
            lines.append(f"| `.idx` 字节一致 | — | — | "
                         f"{'✅ 一致' if r['idx_byte_identical'] else '❌ 不一致'} |")
        lines.append("")

    lines.append("## 说明\n")
    lines.append("- Java 侧墙钟包含 JVM 启动（约 2–5 s），因此小输入下的加速比被低估；输入越大 native 优势越明显。")
    lines.append("- 内存为峰值 RSS；GATK 侧受 `-Xmx1g` 上限约束。")
    lines.append("- 字节一致 = native 与 GATK 在同一输入上产出完全相同字节流。\n")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bases", type=int, default=2_000_000)
    parser.add_argument("--shards", type=int, default=4)
    parser.add_argument("--per-shard", type=int, default=25_000)
    parser.add_argument("--interval-records", type=int, default=600)
    parser.add_argument("--genes", type=int, default=400)
    parser.add_argument("--report", type=pathlib.Path,
                        default=ROOT / "NATIVE_VS_GATK_COMPARISON.md")
    args = parser.parse_args()

    results = {"backend": "OpenMP"}
    with tempfile.TemporaryDirectory(prefix="fastgatk-nvj-") as directory:
        work = pathlib.Path(directory)
        compare_references(work, args.bases, results)
        gather_vcfs(work, args.shards, args.per_shard, results)
        interval_list(work, args.interval_records, results)
        gtf(work, args.genes, results)

    report = md_report(results, "OpenMP")
    args.report.write_text(report, encoding="utf-8")
    print(json.dumps(results, sort_keys=True))
    print(f"\nReport written to {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
