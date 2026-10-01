> **范围裁决（2026-09-28）**：医学/肿瘤应用面（gCNV 种系 CNV 轴、Funcotator 报告/整合余项、
> STR 位点字节 pin、A1/cloud-URI 散项）暂不纳入项目目标；本清单中对应条目降级为"范围外存档"。
> 通用基因计算对齐目标已达成（50 工具、语料全绿）。

# GATK 范围内尚未完成的重要算法

日期：2026-09-16  
对照：pinned GATK **4.6.2.0**（`gatk-source` + `gatk-package-4.6.2.0-local.jar`）  
口径：只谈 GATK 工具箱里的算法，不谈 BWA/STAR/VEP 等 GATK 之外的生态。  
状态快照：native registry **48 / 210** 命令（47 contract-compatible + 1 GenomicsDB adapter）；运行时 **无 JVM**；Java GATK 仅作 CI oracle。

---

## 1. 结论

在 GATK 范围内，**germline/somatic 短读 calling 的主干数值核已经有 native 实现**：局部 haplotype 装配、PairHMM、SW、BQSR、二倍体基因型、VQSR VBEM、Mutect2 似然公式、FilterMutectCalls 的 Java-parity 学习路径。  
**还没完成、且会改变生产结果的，主要不是「缺一个 CLI」，而是下面几类算法：**

1. **HaplotypeCaller 装配图的剩余状态机**（完整 SeqGraph / K-best / 多 AssemblyRegion）。
2. **gVCF 联合基因分型语义**（`*` spanning deletion、ReblockGVCF、GenomicsDB 原生 import）。
3. **体细胞 joint 模型的默认路径与线粒体重复区**（公式已齐，输入图和默认学习路径未齐）。
4. **CNV 的完整概率模型**（ModelSegments MCMC/Gibbs、GermlineCNVCaller 整条 HMM）。
5. **整族未注册的 GATK 算法**：Funcotator、VariantAnnotator 插件、CNNScoreVariants、DragSTR、SV、Spark/PathSeq/BWA。

一句话：当前覆盖的是 **GATK Best Practices 短读 SNP/INDEL 流水线的可运行子集**；同一流水线里会改 call set 的剩余算法，集中在装配图、联合分型、CNV 概率模型和注释/SV/CNN 整族缺席。

---

## 2. 对照：GATK 里已经落地的算法

便于后面「未完成」不和「已有骨架」混淆。

| GATK 算法 | Native 落点 | 完成度（对 4.6.2.0） |
|---|---|---|
| ActivityProfile / AssemblyRegion 切窗 | `activity_profile` + Host | 主路径有；任意窗口未证 |
| ReadThreading 局部 de Bruijn、k 循环、dangling recover-all | `kmer_graph` + HC | chr17 VCF / recover-all 双叉已对 GATK；SeqGraph 拓扑未 1:1 |
| Affine-gap Smith–Waterman（haplotype / dangling） | Kokkos score + Host CIGAR | 有 oracle |
| PairHMM（regular；少量 flow） | `pairhmm_kokkos` | 常规路径有；GKL flow 全语料未齐 |
| 贝叶斯基因型 PL/GQ、HW prior | `genotype` + HC | 二倍体常用路径有 |
| BQSR covariate 计数 / ApplyBQSR | `bqsr` | 生产路径可用 |
| Picard MarkDuplicates / SortSam | native 工具 | 合同级可用 |
| VQSR VBEM–GMM + ApplyVQSR scatter/gather | VariantRecalibrator | pinned 2-Gaussian **VQSLOD 位级**；`culprit` 字段仍错 |
| Mutect2 SomaticLikelihoodsEngine 公式 | Mutect2 kernel | 公式机器级；真实窗 TLOD 仍受装配上下文差 ~3.5 |
| FilterMutectCalls SomaticClusteringModel（env Java 路径） | FMC + kernels | env 路径阈值 Δ≈3e-12；**默认路径仍是 native 近似学习** |
| 区间 CNV 计数 / PoN 去噪 / 点估计分段 | CollectReadCounts 等 | 骨架有；Gibbs/MCMC 未齐 |

以上说明：**「工具在 registry」≠「算法 1:1」。** 下文只列仍会改生物学输出或挡住 Best Practices 替换的缺口。

---

## 3. 未完成算法（按 GATK Best Practices 影响排序）

影响定义为：缺了它，同一条 GATK 推荐流水线会得到不同 call set、不同 FILTER、或根本跑不起来。

### 3.1 P0 — 会改 germline/somatic call set 的核心状态机

#### A. 完整 ReadThreadingAssembler / SeqGraph（HC/Mutect2 共用）

GATK 源：`ReadThreadingAssembler`、`ReadThreadingGraph`、`SeqGraph.simplifyGraph`、`GraphBasedKBestHaplotypeFinder`。

| 子算法 | GATK 做什么 | Native 现状 | 为什么重要 |
|---|---|---|---|
| SeqGraph 线性链/钻石/tail 压缩 | 把 k-mer 图收成序列图再 K-best | 有压缩计数和路径语言；**不承诺与 GATK DOT 拓扑字节一致** | 改变进入 PairHMM 的 haplotype 集合 |
| K-best 寻路终止于 reference sink | 按边权重找源到汇路径 | 有界 DFS/路径枚举；与 GATK K-best 未在任意图上证等价 | 漏/多重 haplotype → 漏/假 call |
| recover-all 后 `hasCycles()` 丢弃该 k | recover-all 引入环则该 k 作废、换更大 k | 环检测在建图阶段；**未证明与「恢复后再丢弃」同一时刻** | 换 k 失败会换一套图 |
| 多 AssemblyRegion 连续窗 | 每 region 独立装配+基因型，再按 core 发表 | **已知缺陷**：连续链式多窗会丢后续 region 的 call，甚至等位翻错（chr20 20228/20429 簇） | WGS/WES 真实区间几乎都是多 region |
| 图中缺失的 indel/soft-clip 路径 | EventMap 只从装配 haplotype 来 | 默认 HC 已禁止 raw pileup 独立成等位；**图没装配出来的 indel 仍然没有** | 复杂 indel 敏感 |
| Linked de Bruijn / 人工 haplotype 恢复 | Mutect 实验路径 | 选项存在；完整 junction-tree 未当生产合同 | 体细胞重复区 |

证据：`GATK_ALIGNMENT.md` 装配行仍标 **partial**；`NEXT_PHASE_TASKS.md` A1 剩余「完整 SeqGraph 拓扑/pruning 1:1 仍未证」；A3 多 region 在执行日志里记为架构级缺陷。

chr17 69–70k 和 recover-all 双叉 fixture **不能**外推到 WGS：前者窗口小、后者是构造拓扑。

#### B. HaplotypeCaller genotyping / reference-confidence（gVCF）

GATK 源：`HaplotypeCallerGenotypingEngine`、`ReferenceConfidenceModel`、`AlleleFrequencyCalculator`。

| 子算法 | 缺口 |
|---|---|
| 任意倍性 / `--max-alternate-alleles` 全矩阵 | 有部分 polyploid/gVCF oracle；不是全部 calling mode |
| gVCF 参考置信度块边界、GQ band、homopolymer insertion | 有定向 oracle；真实 shard 组合仍在 Reblock 边界 |
| `--pileup-detection`（`PileupBasedAlleles`） | **缺席（by design）**。GATK 默认也关；打开 beta 时 native 没有对应算法 |
| 完整 annotation 插件集（AS_QD、FS、SOR、InbreedingCoeff…） | 常用 INFO 有一部分；VariantAnnotator 整工具未注册 |

A3（`NEXT_PHASE_TASKS`）明确未完成：字段级 exact corpus 扩到 ≥5 真实 fixture，且依赖 A1。

#### C. Mutect2 装配上下文与线粒体

公式层（`SomaticLikelihoodsEngine`）和 Java 转写 **机器级**：
全精度通道（`FASTGATK_DEBUG_TLOD=1`）下 6 biallelic 位点 max Δ ~5e-13（2026-09-11 重测：17:69631 处 4.69e-13，仍在 `<1e-8` 守门内；R24 记录的 1.4e-14 头条数字当时未复现，标 stale）。`fastgatk-mutect2-tlod-formula` CTest 双后端绿。  
注：`FASTGATK_TLOD_FULL` 只是 verifier 解析的日志 tag，不是 env 通道；真正的 env 通道是 `FASTGATK_DEBUG_TLOD=1`。  
生产差异来自 **送进公式的 haplotype/fragment 集合**：

| 子算法 | 缺口 | 影响 |
|---|---|---|
| 与 HC 共用的装配图 | 同上 A | chr17:69368 在 2026-09-11 重测：**native = GATK = 70.51**（全窗，`-L 17:69000-70000`）；单窗 GATK 自身漂到 67.43，native 也跟随漂到 67.43。R18/R30 时期记录的 Δ≈3.5（native 常数 74.01）在当前二进制下不复现——审计原文 "the previously recorded Δ3.5 TLOD residual at chr17:69368 is NOT reproducible now"。修复可归因于哪次提交未在历史中保留（工作树 clean）。 |
| 线粒体高深度重复（`--mitochondria-mode`） | ✅ **已闭环（2026-09-28）**：行集合全等 11/11，chrM:8372 byte-identical。根因 = somatic 活动证据缺 GATK `PileupQualBuffer` 的逐删除位点 INDEL 观测（C-stretch 7bp 缺失的 Q90×7 是播种位点 1-based 8373 的主导证据），修复见 `work/mutect2-priority/NEXT_STEPS.md` §1c。同会话闭环 OCM、htsjdk `formatVCFDouble` HALF_UP 渲染、FORMAT/AF 平坦先验。剩余共享行字段级差异（INFO/DP 插入簇、AF 小数位、9037/9070 证据集）pin 于 oracle `KNOWN_OPEN_DIVERGENCES`（双向失败）。 | byte-equality gate（含 pin 双向失败） |
| 默认 FilterMutectCalls 学习路径 | 默认仍用 native 近似 INFO 学习；Java-exact 要 `FASTGATK_FMC_JAVA_PASSES=1` | 未设 env 时阈值/模型不是 GATK |
| `SomaticClusteringModel` 默认化 | 完整 Beta+Binomial 峰分裂+BIC+5×EM 在 env 路径齐；默认路径未切换 | 生产 FILTER 阈值 |

A4/A5 的「引擎验收」已过；**任意肿瘤 WGS 上 TLOD/FILTER 1:1 仍被 A1 挡住**——注意 A4 的引擎 parity 与 A1 的装配输入是两条独立边界，前者已机器级封闭，后者仍是 chr17:69368 旧 74.01 残差对应的那条待办线。

---

### 3.2 P1 — 联合调用与 gVCF 收缩（cohort 必经）

GATK Best Practices germline 在 HC 之后是：

`GenomicsDBImport`（或 CombineGVCFs）→ `GenotypeGVCFs` → `VariantRecalibrator`/`ApplyVQSR` → `VariantFiltration`

| 算法 | GATK | Native | 缺口性质 |
|---|---|---|---|
| GenomicsDB 列式联合导入 | Java + GenomicsDB/TileDB 写 workspace | **adapter**：调 C/TileDB，不是从零重写；完整 batch/fd/restart 语义未宣称 | 大规模 cohort 的真实 import |
| GenotypeGVCFs spanning-deletion `*` | 默认 **0** 条 `*`（class-c C3） | **已修（M1.3 dense-materialize oracle, 12/12 strict-gated cases PASS）**：byte-for-byte 与 GATK 4.6.2.0 一致，包括 reverse-trim (`AAAA/AACA → AAA/AAC`)、spanning-locus 合成 (`*` rows at 3, 4 from a deletion record spanning 2-4)、dense-block-start 的 REF 处理、`-L` interval clip。Oracle `verify_genotype_gvcf_dense_materialize_gatk_oracle.py` 已注册 CTest `fastgatk-genotype-gvcf-dense-spanning-loci-gatk-oracle` 并加入 `verify_all.sh` serial gate（138s strict mode）。Wave-0 R36 报告中"exit 1 with 7 violations" 已过时——所有 gated cases 自 wave-0 commit `c971cf3`/`5a334d7` 等之后已 PASS | 闭合 |
| ReblockGVCF `--drop-low-quals` | 丢掉低 QUAL 位点 | **已修（wave-0 Case A）**：re-genotype 在 PL[0] 转换前发生（`ReblockGVCF.java:415-424`），`fastgatk-reblock-gvcf-droplowqual-gatk-oracle` PASS | 闭合 |
| ReblockGVCF 三倍体 diploidize | 收成 `0/0/0` + 4-item PL | **已修（M1.4 byte-equal）**：native `convert_to_ref_block` 在 `shouldBeReblocked` clause2（min-PL genotype no concrete ALT）触发时走 GATK's `changeCallToHomRefVersusNonRef` PL[0]==0 分支，用共享 Kokkos `remap_genotype_pl_kokkos` 子集 PL 到 REF+best-ALT（保留 source best ALT likelihood），并经 `order_reference_block_format` 与 `clear_info_except_end` 对齐 FORMAT 字段顺序与 INFO。`verify_reblock_gvcf_triploid_gatk_oracle.py` strict mode 0 violations（chr17:69000 real-data oracle + `verify_reblock_gvcf.py` 仍 PASS）。Trim/gap/NON_REF-AD 覆盖迁到独立的 trim-kept-variant fixture | 闭合 |
| CombineGVCFs 参考块 / interval | 工具已晋级 | 复杂 symbolic / 云端仍 fallback_boundaries | 中等 |

这是 **队列联合分型** 和 **单样本 gVCF 收缩** 的算法差，不是 CLI 差。对应 1:1 计划里的 P5a/P5b。

---

### 3.3 P1 — VQSR 残留与 CNV 概率模型

#### VQSR

VBEM/GMM 在 pinned 2-Gaussian fixture 上 VQSLOD **token 级相等**（B1）。仍缺：

| 算法 | 缺口 |
|---|---|
| `culprit` 注解 | **已修（wave-0 culprit-fix）**：native 写注解名（与 GATK `VariantRecalibratorEngine.calculateWorstPerformingAnnotation` 同语义），pinned oracle `verify_variant_recalibrator_culprit_gatk_oracle` 双后端 PASS（含 60-record rich-per-datum fixture + 6-record fixture）。仍缺完整多高斯 WGS 训练收敛与 AS 模型 |
| 完整多高斯 / 真实 WGS 训练收敛 | 小 fixture 位级 ≠ 全基因组 VQSR 合同 |
| allele-specific VQSR | 选项面有；完整 AS 模型未当 1:1 |
| ApplyVQSR 默认 cutoff 以外的 tranche 语义 | scatter→gather→apply 本地合同已通（B2）；复杂 AS + 真实 tranche 仍薄 |

`culprit` 不改 FILTER 也能过目视，但破坏与 GATK 报告/下游脚本的字段契约。

#### CNV（GATK `tools/copynumber`，19 个命令）

已注册的是计数/去噪/点估计分段：**CollectReadCounts、DenoiseReadCounts、CreateReadCountPanelOfNormals、ModelSegments、CallCopyRatioSegments、CollectAllelicCounts、Annotate/Filter/PreprocessIntervals**。

| 算法 | GATK | Native | 重要性 |
|---|---|---|---|
| ModelSegments **slice-sampler / latent-indicator MCMC / credible interval** | Java Gibbs | **仅 conditional-point**（B3 未做） | 体细胞 CNV 分段边界会变 |
| PoN HDF5 v7 **字节同一** + Spark metadata | 精确 schema | 未宣称 byte identity（B4） | 跨 GATK/native 不能混用 PoN |
| **GermlineCNVCaller** 整条（contig ploidy HMM + gCNV） | 19 命令中的核心 | **未注册** | germline CNV Best Practices 整段缺失 |
| PostprocessGermlineCNVCalls / DetermineGermlineContigPloidy / 绘图工具 | 配套 | 未注册 | 同上 |

没有 GermlineCNVCaller 的 HMM，就不能声称覆盖 GATK germline CNV 推荐。

---

### 3.4 P2 — 过滤、选择、JEXL、注释

| 算法 | 缺口 | 典型后果 |
|---|---|---|
| VariantFiltration **allele-specific FilterStatus** | **已修（wave-0）**：native 用 GATK 的 `SITE\|SITE` 编码（含 `--apply-allele-specific-filters` flag-only 路径），`fastgatk-variant-filtration-asfilterstatus-gatk-oracle` 与 `fastgatk-variant-filtration-flag-only-gatk-oracle` 双后端 PASS | 闭合 |
| SelectVariants **保留 ALT=`.` hom-ref** | **已修（wave-0 selectvariants-refonly）**：native 在 `--remove-unused-alternates` 下保留 hom-ref，`fastgatk-select-variants-refonly-gatk-oracle` PASS；REF-only 行的 ExcessHet/MQ/PGT/PID/PS 等 raw 修复已落地（git `5a334d7`、`554cc50`、`25c9558`） | 闭合 |
| **完整 JEXL** 引擎 | 简单谓词有；复杂 JEXL / concordance graph 仍 fallback_boundaries | `--filter-expression` 生产脚本 |
| **VariantAnnotator** + 80+ 注解插件 | 工具 **未注册** | 不能做 GATK 式独立注解 pass |
| **Funcotator**（基因/转录本/MAF） | 4 个命令整族 **未注册**；字符串+数据源查询 | somatic Best Practices 注释步缺失 |

Funcotator 在 GATK 里不是数值 kernel，但是 **体细胞推荐流水线的标准注释器**。缺它等于 Best Practices 在 FilterMutectCalls 之后断掉。

---

### 3.5 P2 — 整族缺席、但属于 GATK 的重要算法

这些不在 48 工具里，却是 GATK 4.6.2.0 源码里的独立算法族（见 `GATK_MODULE_INVENTORY.md`）。

| 算法族 | GATK 命令（代表） | 计算内容 | 在 Best Practices 中的位置 |
|---|---|---|---|
| CNN 替代 VQSR | `CNNScoreVariants`、`FilterVariantTranches` | 1D/2D CNN 对变异打分 | 可选 germline 过滤 |
| DragSTR / DRAGEN-GATK | `CalibrateDragstrModel`、`ComposeSTRTableFile` | STR 误差模型，喂给 HC | DRAGEN-GATK germline |
| RNA | `SplitNCigarReads`、`ASEReadCounter` | 剪接 CIGAR 拆分、等位表达 | RNAseq Best Practices |
| SV | `SVCluster`、`SVAnnotate`、`SVConcordance`、`SVStratify` + `CondenseDepthEvidence` | 证据聚类与分层 | SV 流水线 |
| 比对 / Spark | `BwaSpark`、`MarkDuplicatesSpark`、约 40 个 Spark 命令 | BWT 比对、分布式重复标记、PathSeq | 不复刻 Spark；native 也 **没有 BWA** |
| PathSeq | `PathSeqPipelineSpark` 等 | 微生物 k-mer 分类 | 感染组 |
| 线粒体专用配套 | `FilterAlignmentArtifacts` 等 | 比对伪影过滤 | 线粒体 calling 后处理 |
| RampedHaplotypeCaller | 同 HC 包第二条命令 | 未注册 | 少用 |

其中对「常用 GATK 流程」最关键的缺席是：**Funcotator、GermlineCNVCaller、CNNScoreVariants、SplitNCigarReads、SV 聚类、DragSTR。**  
Spark/BWA/PathSeq 按项目合同是 fallback-first，且当前 **不允许 Java fallback**，所以这些命令现在是 **直接不可用**，不是「调 jar」。

---

### 3.6 P3 — 支撑算法（不单独出 call，但改变数值）

| 算法 | 状态 |
|---|---|
| **BAQ**（Base Alignment Quality HMM） | 计划里 hybrid；BQSR/HC 路径未宣称与 GATK BAQ 一致。关 BAQ 的流水线不受影响；开了会改质量 |
| PairHMM **flow** 全 GKL corpus | 4 个 pinned case raw-bit；从零构造 flow read 的 Java harness 未建 |
| PairHMM/SW **GPU**（CUDA/HIP） | Track D 全空；现仅 OpenMP/Serial CPU |
| 云 I/O（GCS NIO） | 未做；HTSlib 本地/POSIX |
| 完整 Tribble 索引语义 | IndexFeatureFile 部分格式字节一致；ValidateVariants 完整字典/GVCF annotation 仍有边界 |

---

## 4. 映射到三条 GATK 推荐流水线

下面用「这一步的算法能不能替换 GATK」而不是「有没有同名工具」。

### 4.1 Germline short-read WGS/WES（最常用）

```
MarkDuplicates → SortSam → BQSR → HaplotypeCaller [-ERC GVCF]
  → GenomicsDBImport/CombineGVCFs → GenotypeGVCFs
  → VariantRecalibrator/ApplyVQSR → VariantFiltration
  → (Funcotator 可选)
```

| 步 | 替换判断 |
|---|---|
| MarkDuplicates / SortSam / BQSR | 可用（合同级） |
| HaplotypeCaller 小窗 SNP | chr17 级可用 |
| HaplotypeCaller **多 region / 复杂 indel / 完整图** | **不能当 1:1** |
| HaplotypeCaller gVCF | 定向 oracle 有；Reblock/联合 `*` 仍差 |
| GenomicsDBImport | adapter，非完整 native |
| GenotypeGVCFs | 常用位点可用；**`*` 行为反了** |
| VQSR 打分 | 小模型位级；真实 WGS + `culprit` 未齐 |
| VariantFiltration | 简单 FILTER 可用；AS 字段差 |
| Funcotator | **没有** |

### 4.2 Somatic short-read（Mutect2）

```
Mutect2 → GetPileupSummaries → CalculateContamination
  → LearnReadOrientationModel → FilterMutectCalls → Funcotator
```

| 步 | 替换判断 |
|---|---|
| Mutect2 工具链 | 在；DREAM/HCC1143 选定 fixture 严格过 |
| 任意肿瘤 + 线粒体 | 装配/重复区 **未齐**；TLOD 可差数个 log10 |
| FilterMutectCalls | env Java 路径 1:1；**默认路径不是 GATK 学习** |
| Funcotator | **没有** |

### 4.3 Germline CNV

```
PreprocessIntervals → CollectReadCounts → (PoN)
  → DetermineGermlineContigPloidy → GermlineCNVCaller → PostprocessGermlineCNVCalls
```

计数/区间预处理有；**ploidy HMM 和 gCNV caller 整段没有。**  
体细胞 CNV 的 ModelSegments 有点估计，**没有 GATK 的完整 Gibbs。**

### 4.4 RNAseq / SV / DRAGEN-GATK / CNN

对应算法族均 **未注册**。当前 native 不能接这些 Best Practices。

---

## 5. 未完成算法的性质分类

避免把「缺工具」和「工具在但算法错」混在一张表里。

### 5.1 工具在，算法未 1:1（最高优先级）

这些已经在 48 工具里，生产上最容易「看起来能跑、结果和 GATK 不同」。

1. ~~HC/Mutect2 **装配图 + 多 AssemblyRegion**~~ — **M1.5 closed**：native 的 graph window merge (`calling_pipeline.cpp:15891-15895`) 现在受 `max_assembly_region_size + 2*assembly_region_padding` 预算约束。`fastgatk-hc-chr20-real-contract` contiguous window 19901-10020710 byte-equal GATK 4.6.2.0
2. GenotypeGVCFs dense-materialize：emitted-ownership 已修（M1.1）+ 12/12 strict-gated cases PASS（M1.3）
3. ReblockGVCF **三倍体 diploidize** — **M1.4 closed**
4. Mutect2 **真实窗 TLOD（装配输入）** 与 **线粒体**
5. FilterMutectCalls **默认学习路径**
6. ModelSegments **MCMC**

**wave-0 已闭合（不再在本表）：**

- ~~VariantRecalibrator **`culprit`**~~ — wave-0 culprit-fix，`fastgatk-variant-recalibrator-culprit-gatk-oracle` 双后端 PASS（注解名语义对齐 `VariantRecalibratorEngine.calculateWorstPerformingAnnotation`）
- ~~VariantFiltration **AS_FilterStatus**~~ — wave-0，`SITE\|SITE` 编码与 flag-only 路径对齐，`fastgatk-variant-filtration-asfilterstatus-gatk-oracle` 与 `fastgatk-variant-filtration-flag-only-gatk-oracle` 双后端 PASS
- ~~SelectVariants **hom-ref**~~ — wave-0 selectvariants-refonly，`fastgatk-select-variants-refonly-gatk-oracle` PASS；REF-only 行的 ExcessHet/MQ/PGT/PID/PS 同步闭合
- ~~ReblockGVCF **drop-low-quals**~~ — wave-0 Case A，`fastgatk-reblock-gvcf-droplowqual-gatk-oracle` PASS

### 5.2 工具不在，但是 GATK 推荐流水线的算法

1. Funcotator（注释）
2. GermlineCNVCaller + contig ploidy HMM
3. CNNScoreVariants
4. SplitNCigarReads
5. DragSTR
6. SVCluster / SVAnnotate
7. VariantAnnotator

### 5.3 合同上故意不做或延后

- 210−48=**165** 个 GATK 命令：fail-closed，不实现、不调 jar。
- Spark 约 40 个命令：不复刻 Spark DAG。
- BWA 比对：GATK 自己也是调外部比对器；native 没有 BWT/FM-index。
- GPU PairHMM/SW：正确性 Track D 未开工。

---

## 6. 建议推进顺序（只按算法，不按「再加 CLI」）

若目标是 **GATK 范围内常用生产路径可替换**，而不是凑工具数：

1. **多 AssemblyRegion 独立装配+基因型**（修 chr20 连续窗丢 call / 等位翻转）。没有这一步，WGS 级 HC/Mutect2 1:1 不成立。
2. **SeqGraph + K-best 与 GATK 同图同路径**（至少再钉 ≥5 个真实窗，含 indel 密集区）。这是 TLOD Δ3.5 和复杂 indel 的根。
3. **GenotypeGVCFs `*` + ReblockGVCF C1/C2**。联合分型位点集才能和 GATK 比。
4. **FilterMutectCalls 默认路径切到 Java 学习模型**（现成 env 路径）。
5. **ModelSegments MCMC** 或直接补 **GermlineCNVCaller**（看 germ vs soma CNV 哪条要上生产）。
6. **Funcotator**（体细胞推荐的最后一公里，算法是数据源查询不是 kernel）。
7. 其余：JEXL、CNN、DragSTR、RNA、SV、GPU。

---

## 7. 和本仓库其它文档的关系

| 文档 | 关系 |
|---|---|
| `GATK_ALIGNMENT.md` | HC/Mutect2/GenotypeGVCFs **边界级** partial/aligned 账本 |
| `NEXT_PHASE_TASKS.md` | Track A/B/C/D 任务状态；A1 ◐，A3–A6 / B3–B5 / D 未完 |
| `A_B_EXECUTION.md` | Mutect2 TLOD / FMC 学习路径的实验结论 |
| `fastgatk-native/evidence/class-c-parity-bugs.md` | 工具在但结果钉反了的 6 条（C7 已关） |
| `GATK_MODULE_INVENTORY.md` | GATK 210 命令按包的建议状态 |
| 本文 | **按算法**回答：GATK 范围内还缺什么、缺了会怎样 |

本文不替代 A/B 任务清单，也不把 P3–P10 的 CLI 长尾算进「重要算法」。重要 = 会改 GATK Best Practices 的 call set、FILTER、gVCF 块或 CNV 分段。

---

## 8. 一页对照

```
GATK 短读变异核心
├─ 已有（可跑、部分 oracle）
│    ActivityProfile, 局部 k-mer 图, SW, PairHMM, 二倍体 PL/GQ,
│    BQSR, MarkDuplicates, VQSR VBEM, Mutect2 公式, FMC env 学习
├─ 未完成且会改结果          ← 当前真正的算法债
│    完整 SeqGraph/K-best, 多 region 装配,
│    gVCF * / Reblock, 线粒体图, FMC 默认学习,
│    ModelSegments MCMC, GermlineCNV HMM
└─ GATK 有、native 整族没有
     Funcotator, VariantAnnotator, CNN, DragSTR, RNA, SV, Spark/BWA
```
