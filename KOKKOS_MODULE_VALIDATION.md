# Kokkos 四类代表模块验证报告

## 结论

统一的 Kokkos 数据/执行模式已经在四类负载上跑通：规则 DP、统计归约、
不规则图和 region 流水线。generic、ZEN3/AVX2、ZEN4/AVX-512 三个 CPU
变体使用同一份源码，所有验证指标跨变体一致。

这证明的是架构模式可复用，不是四个 GATK 工具已经完成兼容。每个模块的
GATK 语义缺口见文末。

## 正确性和确定性

命令：

```bash
bash kokkos-modules-demo/scripts/build_variants.sh
python3 kokkos-modules-demo/scripts/verify_modules.py --records=512 --threads=4
```

结果：`status=pass`，12 个组合（3 个 CPU 变体 × 4 个模块）全部通过，
并且跨变体输出签名一致。验证对象如下：

| 模块 | 验证内容 |
|---|---|
| Smith-Waterman | affine-gap 最大分数 checksum；match=10、mismatch=-15、gap-open=-30、gap-extend=-5 |
| BQSR | observations、errors、quality sum、quality×cycle×context bin checksum |
| local assembly graph | k-mer edge 数、branching node 数、入/出度 checksum |
| region stream | 覆盖位置、错配数、quality/mapQ 总量；同时验证文件输入适配 |

完整机器可读结果：[verify-modules.json](</home/turing-agents/Documents/fast-gatk/kokkos-modules-demo/results/verify-modules.json>)。

## CPU 基准

命令：

```bash
python3 kokkos-modules-demo/scripts/benchmark_modules.py \
  --variant=avx512 --records=2048 --iterations=5 --repeats=3 \
  --cores=1,2,4,8,16 \
  --output=kokkos-modules-demo/results/benchmark-avx512.json
```

主机为 AMD Ryzen 9 7945HX；进程固定到相同的物理 CPU，三次运行取按时间排序的中位数。
吞吐量包含 5 次 kernel iteration；region 的 I/O 单独计时。

| cores | SW (M cell updates/s) | BQSR (M obs/s) | graph (M edges/s) | region (M records/s) |
|---:|---:|---:|---:|---:|
| 1 | 354.7 | 217.7 | 307.0 | 80.0 |
| 2 | 712.9 | 223.4 | 143.6 | 64.7 |
| 4 | 1,404.3 | 244.9 | 175.9 | 81.4 |
| 8 | 2,676.6 | 270.2 | 230.5 | 81.7 |
| 16 | 4,938.9 | 91.7 | 95.1 | 27.1 |

解释：Smith-Waterman 的规则 DP 在本机接近线性扩展；BQSR 的共享 bin
原子更新、图的竞争度数更新和很小的 region batch 在高核数受到同步与调度
开销限制。这些结果支持生产实现使用局部归约、长度/工作量分桶、较大 batch
和持久化 workspace，而不是把 PairHMM 的粒度照搬到所有算法。

完整结果：[benchmark-avx512.json](</home/turing-agents/Documents/fast-gatk/kokkos-modules-demo/results/benchmark-avx512.json>)。

## 当前边界

- Smith-Waterman 目前只验证最大分数；GATK 的 backtrack、CIGAR、overhang 和 tie-breaking 尚未接入。
- 这里的 BQSR 四模块原型仍只验证 Kokkos 统计/分桶；native 路径已另外提供
  GATKReport v1.1/RecalTable0/1/2、可合并 covariate sidecar，以及
  exact→read-group/quality→context/quality→cycle/quality→global quality
  的 ApplyBQSR 层级回退，但完整 GATK recalibration model/rounding 仍未完成。
- assembly 是 k-mer edge/degree/branching 原型，尚未实现完整 ReadThreadingGraph、路径搜索和 haplotype 生成。
- region 输入是 `position base_code quality mapq` 的窄适配格式，用来验证 region 分块和 I/O/compute 解耦；当前环境没有 htslib/samtools，因此尚未测真实 BAM/CRAM 解码。
- 当前只在 CPU OpenMP 后端实测。CUDA/HIP/SYCL 需要可用设备和编译器；GPU 版应为这些模块重新组织 Team/Vector/tile kernel，不能把 CPU 结果当成 GPU 性能结论。native CMake 现在提供 `fastgatk-kokkos-backend-config` 配置门禁：它在不创建 GPU context 的前提下检查 `Kokkos_ENABLE_*`、`Kokkos_ARCH_*`、生成的 `KokkosConfigCommon.cmake` 和编译器 provenance；因此 CPU-only CI 可以发现 stale/mixed backend cache，而不会把“没有设备”误报成 GPU 性能结果。

## 2026-08-21 native 组合回归快照

本机为 AMD Ryzen 9 7945HX（16 physical/32 logical，单 NUMA），CPU flags
包含 `avx512f/avx512dq/avx512bw/avx512vl`；当前 native build 使用
Kokkos 5.2.0 + OpenMP。以下门禁在同一次构建中通过：

| 门禁 | 结果 |
|---|---|
| CTest | 17/17 passed |
| dispatcher registry | 72/72 checks passed |
| PairHMM Strict oracle | scalar/AVX2/AVX-512，4096/4096 raw-bit diff=0 |
| SW GATK overhang oracle | SOFTCLIP/INDEL/LEADING_INDEL/IGNORE，Java tie priority、CIGAR 和 Kokkos/Host score 一致 |
| BQSR report/apply fallback | GATKReport v1.1、sidecar Gather 合并、无 sidecar 的 RecalTable1/2 回退通过 |
| GATK call-set oracle | sentinel 通过，Jaccard=1.0；`bit_identical=false`（按设计） |
| Nextflow/SLURM | local DSL2-equivalent + resource simulation passed |

文件边界基准（包含 HTSlib decode、Kokkos prepare/execute、VCF 压缩和 Tabix）
在 `17:69000-70000` 上测得：HC gVCF wall time 约 5.57 s、Mutect2 wall time
约 11.87 s。HC telemetry 中 PairHMM 为 5228 pairs，prepare 约 1.62 ms，
execute 约 1.19 s；这些数字是当前 native 原型的可重复基线，不代表相对 Java
GATK 的端到端加速结论，后者必须在相同线程、内存、压缩和输入下另行对照。
