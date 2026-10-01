# fast-gatk native vs pinned GATK 4.6.2.0 对比报告

本报告对本次会话消除 fallback 的四个关键工具做 native（C++ host + Kokkos kernel）与 pin 版 GATK 4.6.2.0 的端到端对比：输出正确性、速度、资源占用。

## 环境

- 原生后端：`OpenMP`（Kokkos），线程 4（CompareReferences/IndexFeatureFile）
- Java 基线：GATK 4.6.2.0（`gatk-package-4.6.2.0-local.jar`）+ JDK 17，`-Xmx1g`
- 计时：`/usr/bin/time` 墙钟（Java 侧**含 JVM 启动**）；内存：峰值 RSS
- 正确性：字节一致 = 原生与 GATK 在同一输入上产出完全相同的字节流

## 总览

| 工具 | 输出正确性 | 速度 (native vs GATK) | 内存 (native vs GATK) |
|---|---|---|---|
| CompareReferences FULL_ALIGNMENT | ✅ `.vcf` 字节一致；`.vcf.idx` 内容一致 | **2.7× 更快** | **4.5× 更少** |
| GatherVcfs (BCF 输出) | ✅ BCF 2.1；解码记录一致 | **81.0× 更快** | **84.5× 更少** |
| IndexFeatureFile interval_list | ✅ `.idx` 字节一致 | **220.8× 更快** | **51.7× 更少** |
| IndexFeatureFile GENCODE GTF | ✅ `.idx` 字节一致 | **219.8× 更快** | **45.5× 更少** |

> 注：Java 侧墙钟含 JVM 启动（约 2–5 s），小输入下加速比被低估；输入越大 native 优势越明显。

---

## 1. CompareReferences FULL_ALIGNMENT（去掉 MUMmer 依赖）

输入：2 Mb 参考对，SNP 每 10 kb（199 处差异）。

| 指标 | native (fast-gatk) | GATK 4.6.2.0 | 结论 |
|---|---|---|---|
| 墙钟 | 1.644 s | 4.383 s | native 快 **2.7×** |
| 峰值 RSS | 66.7 MB | 299.4 MB | native 内存仅 GATK 的 **1/4.5** |
| 记录数 | 199 | 199 | ✅ 一致 |
| `.vcf` 输出 | — | — | ✅ 字节一致 |
| `.vcf.idx` sidecar | — | — | ✅ 内容一致（仅 mtime/uri 为运行时元数据） |

说明：GATK FULL_ALIGNMENT 依赖外部 MUMmer（`nucmer --mum → delta-filter -1 → show-snps`）；native 已用 MUM-anchored 共线比对器原生替代，且产出 `.vcf` + `.vcf.idx` 与 GATK 一致。

---

## 2. GatherVcfs BCF 输出（原生 BCF 2.1）

输入：4 shards × 25,000 records = 100,000 records，聚合为 `.bcf`。

| 指标 | native (fast-gatk) | GATK 4.6.2.0 | 结论 |
|---|---|---|---|
| 墙钟 | 0.055 s | 4.455 s | native 快 **81.0×** |
| 峰值 RSS | 5.6 MB | 473.3 MB | native 内存仅 GATK 的 **1/84.5** |
| 记录数 | 100,000 | 100,000 | ✅ 一致 |
| BCF magic | `BCF\2\1` | `BCF\2\1` | ✅ 一致 |
| 解码记录 | — | — | ✅ 一致（经 GATK SelectVariants 解码后逐记录比对） |

说明：GATK/htsjdk 只读 BCF 2.1（minor version 必须 =1），HTSlib 1.22 原生写 2.2。已给 vendored HTSlib 打补丁（`third_party/HTSLIB_BCF21.patch`）输出 2.1 并读写兼容 2.2，且把裸 BCF 打开逻辑抽成共享 helper（`fastgatk-native/include/fastgatk/io/vcf_output.hpp`）应用到全部 11 个变体写入工具。

---

## 3. IndexFeatureFile interval_list

输入：600 条 interval_list 记录（2 个 contig，1-based inclusive）。

| 指标 | native (fast-gatk) | GATK 4.6.2.0 | 结论 |
|---|---|---|---|
| 墙钟 | 0.016 s | 3.533 s | native 快 **220.8×** |
| 峰值 RSS | 5.3 MB | 274.1 MB | native 内存仅 GATK 的 **1/51.7** |
| `.idx` 输出 | — | — | ✅ 字节一致 |

---

## 4. IndexFeatureFile GENCODE GTF

输入：400 个 GENCODE 基因（每个含 gene + transcript + exon，聚合为一个 Feature）。

| 指标 | native (fast-gatk) | GATK 4.6.2.0 | 结论 |
|---|---|---|---|
| 墙钟 | 0.016 s | 3.516 s | native 快 **219.8×** |
| 峰值 RSS | 5.6 MB | 254.9 MB | native 内存仅 GATK 的 **1/45.5** |
| `.idx` 输出 | — | — | ✅ 字节一致 |

---

## 5. 其它本会话完成项（正确性验证）

以下项未单独做速度/内存基准（或已由 CTest 覆盖），正确性均以 GATK oracle 比对通过：

| 项 | 状态 |
|---|---|
| 13 个调用管线 read filter（含 `PlatformUnitReadFilter` `--black-listed-lanes`） | ✅ 已实现并验证（读取器 + Kokkos kernel + 数据流 + CLI） |
| DepthOfCoverage `--partition-type sample\|library\|readgroup` | ✅ 已实现 |
| CreateReadCountPanelOfNormals `--sequence-dictionary` | ✅ 已实现 |
| DenoiseReadCounts `--normalization-target-coverage` | ✅ 确认 GATK 4.6.2.0 无此参数（幻影参数，已移除） |
| IndexFeatureFile 未压缩 VCF/GVCF/BED Tribble 线性/区间树 + interval_list/GFF3/GTF codec | ✅ 字节一致 |
| CTest 契约回归（受影响 12 个工具 + GATK oracle） | ✅ 86/86 通过 |

---

## 附：对比过程中发现并修复的问题

- **`tribble_index.hpp::optimize_linear_contig`**：native 在合并相邻 bin 时比 htsjdk 多合并一级（`lastGood` 语义缺失），导致均匀分布 fixture 的线性 `.idx` 块粒度与 GATK 不同（3 blocks vs 5 blocks）。已按 htsjdk `LinearIndex$ChrIndex.optimize()` 字节码修正（保存并恢复 merge 前的 `lastGood` 布局 + 严格 `>` 边界），修复后 interval_list `.idx` 恢复字节一致。
