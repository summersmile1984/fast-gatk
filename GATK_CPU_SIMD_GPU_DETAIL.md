# GATK 各模块 CPU / SIMD / GPU 计算细化与坑点

这份文档补充主报告，重点回答三个问题：

1. 每个模块到底计算什么、数据如何布局、并行粒度是什么；
2. CPU、SIMD、GPU 分别应该做哪一部分；
3. 哪些地方不能为了“上 GPU”而破坏 GATK 结果、排序、I/O 或资源契约。

## 0. 统一执行模型

### 三种并行层次

| 层次 | Kokkos 对应 | 适合的工作 | 不适合的工作 |
|---|---|---|---|
| 作业/样本/interval | Nextflow + SLURM | 样本 scatter、contig、interval、retry/resume | 在 kernel 内再启动 `sbatch` |
| batch/team | `TeamPolicy`、`TeamThreadRange` | 一个 assembly region、一批 read-haplotype pair、一个矩阵 tile | 依赖其他 team 的锁或等待 |
| lane/vector | `ThreadVectorRange`、`Kokkos::Experimental::simd` | 相同公式、固定长度、连续数组 | 可变长度字符串、随机 map/set、I/O |

Kokkos 的 team 可以共享 scratch，但团队之间没有可依赖的前进保证；不能用 spin-lock 等待另一 team。`TeamThreadRange`/`ThreadVectorRange` 的 collective 调用必须让所有 team 成员走相同分支。参考：[hierarchical parallelism](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/HierarchicalParallelism.html)、[TeamThreadRange](https://kokkos.org/kokkos-core-wiki/API/core/policies/TeamThreadRange.html)。

### CPU SIMD 的统一规则

- 对固定长度或分桶后的数据用 SoA；不要对 Java 风格 AoS/指针对象直接套 SIMD。
- lane 之间不能有写后读依赖；有依赖的 DP 用“跨 pair 向量化”或 anti-diagonal，而不是错误地向长度方向 vectorize。
- 尾部用 masked load/store 或 scalar cleanup，不读取 padding 之外的内存。
- `Kokkos::Experimental::simd<T>` 可为 AVX2/AVX-512/NEON 等提供统一源码，但 ABI/宽度是编译配置的一部分，需要按 `KOKKOS_ARCH` 生成不同 CPU 产物；它不是同一 fat binary 中任意 AVX2/AVX-512 runtime dispatch 的直接替代。启用 GPU backend 时默认 `simd<T>` 还可能成为单 lane，device kernel 应主要使用 Kokkos team/vector parallelism。参考：[Kokkos SIMD](https://kokkos.org/kokkos-core-wiki/API/simd/simd.html)。
- 浮点归约的 lane 顺序会改变最后几位；严格模式要固定树形归约，不能把 SIMD reduction 直接当作字节兼容。

### GPU 的统一规则

- 一个 team 对应一个独立 work unit；长度差异大的 work unit 先分桶。
- host→device 拷贝、kernel、device→host 拷贝使用同一 execution-space instance 或显式 event；不能在每条 read 后 `deep_copy`。
- 预分配 batch arena 和 scratch，禁止 device kernel 内部动态 `new`/`std::vector`。
- GPU 端尽量输出固定大小的数值结果；变长 CIGAR、变长 haplotype、VCF 字符串留在 host。
- 异步 `deep_copy` 在显式 fence 前不会完成；计时时必须 fence 正确的 execution space。参考：[deep_copy](https://kokkos.org/kokkos-core-wiki/API/core/view/deep_copy.html)、[fence](https://kokkos.org/kokkos-core-wiki/API/core/parallel-dispatch/fence.html)。

## 1. Engine：interval、shard、assembly region

源码：[AssemblyRegionWalker.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/engine/AssemblyRegionWalker.java>)、[Shard.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/engine/Shard.java>)。

### 计算方式

对每个 contig/interval：

1. 以 core interval 为输出范围，以 `assemblyRegionPadding` 扩展成 halo；
2. 按参考坐标读取 reads；
3. 每个位点计算 activity probability；
4. 用阈值、最大传播距离、最小/最大 region size 形成 active/inactive region；
5. 只把 core 内的 variant/block 写出。

`ActivityProfile` 要求状态连续，且 soft-clip 可能把概率传播到相邻位置；因此不是简单的逐位点独立 map。

### CPU / SIMD / GPU

- **CPU**：一个 contig 一个 streaming state machine；负责 read iterator、interval 边界、downsampling、active-region 状态转换。
- **SIMD**：在已经解码的 read batch 上并行做 base/quality/pileup 计数、mapping-quality/filter predicate、soft-clip 标记；状态转移仍按坐标顺序在 CPU 完成。
- **GPU**：只把固定窗口的 pileup/activity 统计放到 device；每个 window 返回概率数组和候选计数，CPU 再按顺序合并 profile。

### 坑点

- halo 不能当作输出区，否则 scatter/gather 会产生重复变异；
- 不同输入顺序不能改变 downsampling，deterministic 模式需要按 coordinate/read name 固定 tie-break；
- contig 首尾会截断 padding，不能越界读 reference；
- region 状态机有跨窗口依赖，不能让 GPU window 独立决定边界；
- 预取队列过大可能把 page cache 吃光，预算要包含 reads、reference、index 和 decoded buffer。

## 2. HTS I/O、read decode、filters、BAQ

源码：[ReadWalker.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/engine/ReadWalker.java>)、[GATKTool.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/engine/GATKTool.java>)。

### 计算方式

输入链路为：

```text
BGZF/CRAM block -> SAM record -> CIGAR/reference projection
  -> pre-filter transform -> read filters -> post-filter transform
  -> tool apply()
```

filter 组合是 AND/OR/NOT predicate；BAQ 根据 read/reference/CIGAR 计算 alignment uncertainty；很多工具需要按 coordinate 保序。

### CPU / SIMD / GPU

- **CPU**：HTSlib/BGZF/CRAM 解码、CIGAR 解析、reference range read、索引查询、过滤控制流。
- **SIMD**：批量比较 bases、质量阈值、mapping quality、CIGAR op 分类、BAQ 中的独立概率/质量变换；用 byte/uint16 SoA。
- **GPU**：只有在已经形成大块 fixed-layout reads 时做 base-quality transform、pileup histogram、简单 filter mask；不要把 BAM record 指针直接传给 GPU。

### 坑点

- CRAM 解码依赖 reference 和 codec 语义，不能只支持 BAM 就宣称兼容；
- BGZF block 边界不是 read 边界，不能按压缩块直接 scatter 计算；
- index seek 过于碎片化时网络 latency 远大于 kernel，应该合并相邻 ranges；
- GPU 传输一条 read 的成本通常超过 filter 本身，低覆盖/短区间应 CPU；
- output writer 的 header、program record、压缩等级和 index 生成必须保留 GATK/htslib 契约。

## 3. PairHMM

源码：[PairHMM.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/utils/pairhmm/PairHMM.java>)、[LoglessPairHMM.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/utils/pairhmm/LoglessPairHMM.java>)、[PairHMMModel.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/utils/pairhmm/PairHMMModel.java>)。

### 计算方式

对每个 read `r` 和 haplotype `h`，建立三状态 DP：match `M[i,j]`、insertion `I[i,j]`、deletion `D[i,j]`。logless 版本的核心为：

```text
M[i,j] = prior[i,j] * (M[i-1,j-1]*T_MM[i]
                       + I[i-1,j-1]*T_IM[i]
                       + D[i-1,j-1]*T_IM[i])
I[i,j] = M[i-1,j]*T_MI[i] + I[i-1,j]*T_II[i]
D[i,j] = M[i,j-1]*T_MD[i] + D[i,j-1]*T_DD[i]
```

最终只累加最后一个 read row 的 `M` 和 `I`，忽略以 deletion 结束的路径；`prior` 由 base quality、`N` 和三态修正决定。log10 版本用 log-sum-exp，避免下溢但更慢。

GATK 还按 haplotype 变化位置复用 read cache：只有 haplotype 前缀改变后的 cell 需要重新计算；若长度不同或差异位置倒退，必须重新 cache。

### CPU

- 一个 worker 处理一个或多个 read-haplotype pair；预分配 `M/I/D` 行列。
- 短 read（通常 100–300 bp）用两行滚动数组可减少内存；需要 haplotype prefix cache 时保留必要的旧行。
- transition/prior 用 0–93 Phred 的查表，避免每 cell 调 `pow/log`。
- NUMA 绑定：每个 worker 的 batch、matrix、quality table 放本地 NUMA node。

### SIMD

最安全的向量化方向是**跨独立 pair**：同一个 lane 处理不同 read-haplotype pair，保证每个 lane 内的 DP 依赖仍是标量顺序。另一种是 anti-diagonal vectorization：

```text
diagonal d: (i,j) where i+j=d
```

同一 diagonal 的 cell 依赖上一/上两条 diagonal，但要处理边界和不同长度 mask。建议先跨 pair，再考虑 diagonal。

`Kokkos::Experimental::simd<double>` 适合 CPU logless/transition 公式；quality lookup 可用 uint8 index + gather。不要假定所有平台都有相同 lane width。

当前可运行实现已经按这条路线完成：`Kokkos::View` + `RangePolicy` + `Kokkos::Experimental::simd<double>`，同源码构建 scalar/AVX2/AVX-512。M/I/D 利用从左到右依赖做三行原位更新：旧 `M/I/D[j]` 在覆盖前保存到 SIMD 寄存器，下一 cell 作为 diagonal 输入，因此只保留三条数组且不改变 `(A+B)+C` 顺序。512 records 的 8×8 matrix 共 4096 个结果在三个构建上都与固定 GATK strict oracle raw-bit 相同。

同机最终数据表明 AVX-512 matrix kernel 为 GKL 实际 `n*n` 吞吐的 1.09–1.28×；若把每次重新打包/分配/deep-copy 摊入十次计算，则为 1.05×、1.02×、1.01×、0.99×、0.94×（1/2/4/8/16 核）。因此 production gate 不只是 kernel：必须池化 workspace 和转换后的 read transitions，避免高核数下准备串行段主导。

### GPU

- 一个 team 处理一个 pair 或一小组相同长度 pair；team scratch 放 tile 化的 `M/I/D`。
- read/haplotype 以 2-bit/uint8 SoA 存储，offset 单独数组；每个 bucket 只包含接近的长度，减少 padding 浪费。
- DP 采用 anti-diagonal 或 tiled wavefront；tile 之间用 team barrier，不使用跨 team spin-lock。
- 输出 likelihood matrix 按 GATK 的 allele-major/read-major 索引顺序写入，host 端再恢复对象映射。
- GPU 只在 pair 数量大、长度分布窄、PCIe/NVLink 拷贝可摊平时启用。

### 坑点

- regular、flow、DRAGstr/non-symmetrical quality model 是不同数学契约，不能共用一个近似 kernel；
- logless 的 scaling/initial condition/tristate correction 必须逐项匹配，不能擅自改成 float；
- `hapStartIndex` cache 依赖 haplotype 顺序，GPU reorder 后必须保存 permutation 并恢复；
- `double` GPU 吞吐可能很低，fast mode 可评估 float，但必须先有误差/variant oracle；
- 每 pair 单独 kernel launch 会被 launch overhead 杀死；
- 原位三行 DP 能保持单 pair recurrence，但会放弃 GATK 的跨 haplotype common-prefix matrix cache；production matrix scheduler 仍应显式实现 prefix cache，并用 raw-bit oracle 验证；
- empty read、read/haplotype 长度约束、无效 log probability 要在 host 端先验证。

## 4. Smith-Waterman

源码：[SmithWatermanAligner.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/utils/smithwaterman/SmithWatermanAligner.java>)。

### 计算方式

局部对齐通常维护 `H`、`E`、`F`：

```text
E[i,j] = max(H[i-1,j] + gap_open, E[i-1,j] + gap_extend)
F[i,j] = max(H[i,j-1] + gap_open, F[i,j-1] + gap_extend)
H[i,j] = max(0, H[i-1,j-1] + match/mismatch, E[i,j], F[i,j])
```

需要记录最高分位置并 traceback 得到 CIGAR；GATK 还区分 SOFTCLIP、INDEL、LEADING_INDEL 等 overhang 策略。

### CPU / SIMD / GPU

- **CPU**：短序列直接滚动行 + traceback；CIGAR 生成和 tie-break 留 CPU。
- **SIMD**：跨独立 alignment request 做 lane-parallel；也可用 striped/wavefront layout，但要保持每 lane 的最大值位置。
- **GPU**：一个 team 处理一个 alignment；score matrix 用 anti-diagonal/tile；先输出 score、argmax、traceback checkpoints，host 再回溯 CIGAR。短序列可把多个 alignment 合并到一个 team。

### 坑点

- 相同分数的选择顺序会改变 CIGAR；必须复制 Java/GKL 的 tie-break，而不是任意 `max`；
- 只比较 score 不足以证明兼容；要比较 CIGAR、start/end、overhang；
- traceback 需要额外内存，不能为了省显存只保留 score 矩阵后再猜路径；
- 变长 CIGAR 写回会造成 GPU 原子分配和不确定顺序，建议固定上限 + host compaction；
- 参数是整数/负分，不能在 SIMD/GPU 中不小心转成 unsigned。

## 5. k-mer、read-threading assembly、De Bruijn graph

源码：[AssemblyBasedCallerUtils.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/AssemblyBasedCallerUtils.java>)、[ReadThreadingAssembler](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/readthreading/ReadThreadingAssembler.java>)。

### 计算方式

1. read error correction（可选）；
2. 取每条 read 的 k-mer；
3. 对 canonical k-mer（正向/反向互补取较小编码）计数；
4. 以 k-mer 为 node、相邻 k-mer 为 edge 建图；
5. prune low-support paths、找 reference-connected paths、生成 haplotypes；
6. 用 SW 对 haplotype/reference 对齐，再用 PairHMM 计算 read likelihood。

### CPU

- 2-bit rolling encoding，使用 per-thread hash/radix buffer；
- graph node/edge 用 flat arena + offset，不使用每 node 一个 heap object；
- 在 NUMA 节点内按 region 建 arena，region 完成后整体释放；
- 图搜索、路径剪枝、异常恢复和 debug dump 保留 host。

### SIMD

- 同一 read 的 rolling k-mer 用 shift/or/掩码；
- 多 read 以 SoA 并行做 base encode、reverse-complement、hash；
- radix sort 的 key compare/partition 用 SIMD；
- edge support/coverage 的整数累加使用 per-thread local bins，避免全局 atomic。

### GPU

- 两阶段：GPU 产生 `(canonical_kmer, read_id, position)`，GPU radix sort + run-length encode；host/设备都可做 unique/count；
- graph node id 由 sorted key 的 index 决定，edge 由相邻 key 查表；
- 只把固定长度的 edge arrays、support counts、candidate score 放 device；复杂 path traversal 在 host；
- GPU 适合很多 reads/region，单个小 active region 不值得传输。

### 坑点

- `N`、ambiguous base、低复杂度 k-mer 的处理必须跟 GATK 一致；
- canonical 方向改变时要保留 strand/orientation 信息，否则 CIGAR/allele orientation 错；
- hash collision 不能当作 key equality；
- graph prune 的 tie-break 必须稳定，否则 haplotype 顺序改变会连锁影响 PairHMM 和 genotype；
- region 内候选 haplotype 数可能爆炸，要在 host 端设置与 GATK 相同的上限和 deterministic ranking；
- 不能把 irregular graph pointer chasing 当作“GPU 适合的并行图”。

## 6. Haplotype likelihood、genotyping、annotation

源码：[HaplotypeCallerGenotypingEngine.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/HaplotypeCallerGenotypingEngine.java>)、[GenotypeLikelihoodCalculator.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/genotyper/GenotypeLikelihoodCalculator.java>)。

### 计算方式

输入是 allele/haplotype × read likelihood matrix。典型 diploid genotype likelihood 是对 read 的 allele likelihood 做乘积/对数求和，并按 genotype prior 计算 posterior/PL；多等位位点还要做 allele subsetting、AF calculation、MLE AC/AF、annotation。

GATK 的实际调用链包含：

```text
haplotype likelihoods
  -> marginalize haplotypes to alleles
  -> retain/filter evidence
  -> genotype likelihoods per sample
  -> allele-frequency calculation
  -> allele subsetting
  -> annotations + VariantContext writer
```

### CPU / SIMD / GPU

- **CPU**：allele mapping、read filtering、sample/allele list、genotype enumeration、VariantContext/annotation/plugin。
- **SIMD**：对 `likelihood[allele][read]` 做 log-add/log-sum-exp、per-read max、AD/DP/strand counts；对同一 ploidy/allele count 的 samples 分桶。
- **GPU**：把大量同形状的 likelihood matrices 做 batched log-sum-exp、genotype combination score、AF reduction；输出固定数组 `(sample, genotype, PL/GQ)`，host 负责 allele ordering 和 VCF fields。

### 坑点

- allele trimming/subsetting 会重排索引；GPU 结果必须携带 permutation；
- log10 likelihood 的归一化、underflow、`-Infinity`、no-call 规则不能用普通 `exp` 替代；
- 多等位 genotype 数组合增长，盲目 GPU 会耗尽显存；应先限制 max alternate alleles；
- annotation 插件有 side effect、字符串和 feature query，不能在 device 运行；
- 同一 likelihood 数值但不同 allele order 可能产生不同 VCF，必须在 host 统一排序。

## 7. Reference confidence、GVCF 与 GVCF writer

源码：[ReferenceConfidenceModel.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/ReferenceConfidenceModel.java>)、[ReferenceConfidenceVariantContextMerger](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/ReferenceConfidenceVariantContextMerger.java>)。

### 计算方式

- 对每个 reference locus 生成 non-ref likelihood、DP、GQ 等；
- 连续且兼容的 hom-ref loci 合并成 GVCF reference block；
- 遇到 variant、GQ band 边界、coverage/annotation 变化或 contig 结束时关闭 block；
- 多 sample GVCF merge 时按 locus 合并 allele/block，再 genotype。

### CPU / SIMD / GPU

- **CPU**：block state machine、locus ordering、GVCF merge、header、writer、index；
- **SIMD**：在 fixed locus batch 上做 base quality/likelihood、DP/GQ 计算；
- **GPU**：只做 per-locus likelihood/quality kernel，返回按 coordinate 排序的固定结果；block 合并永远 host 顺序执行。

### 坑点

- block 的 start/end、END 字段、GQ band 和 contig end 都是兼容性敏感点；
- halo 不能让同一个 block 被两个 shard 写两次；
- GPU batch 不能跨样本/contig 混用错误的 reference block state；
- writer 失败、index 缺失或 header 不一致必须使 tool 失败，不能只生成一个“看似可用”的 VCF。

## 8. Mutect2 somatic calling

源码：[Mutect2.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect/Mutect2.java>)、[Mutect2Engine.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect/Mutect2Engine.java>)。

### 计算方式

Mutect2 复用 active region/local assembly，但增加 tumor/normal、contamination、orientation artifact、germline resource、panel of normals、TLOD/NLOD、F1R2 和 filtering stats。每个 region 的主要数值热点仍然是：

```text
read/haplotype likelihoods
  -> tumor/normal allele evidence
  -> somatic likelihood/LOD
  -> artifact/germline/contamination filters
  -> VCF + .stats + optional F1R2/GVCF
```

### CPU / SIMD / GPU

- **CPU**：region orchestration、tumor-normal sample pairing、filters、F1R2 histogram、stats/VCF sidecars；
- **SIMD**：allele count、base/strand/orientation counts、likelihood vector arithmetic；
- **GPU**：复用 HC PairHMM/SW；对多个 tumor-normal region 做 batched somatic likelihood；
- **fallback**：任何 flow mode、特殊 artifact model 或未验证过滤参数回 Java。

### 坑点

- tumor/normal 样本 ID 和 read group 不能被 SoA packing 错配；
- TLOD/NLOD/contamination 的阈值和 `-Infinity/NaN` 行为要一致；
- `.stats`、F1R2 tar.gz、somatic GVCF 是输出契约，不是可选日志；
- 不能只验证 VCF；Mutect2 downstream FilterMutectCalls 依赖 stats 字段和 header。

## 9. BQSR / ApplyBQSR

源码：[BaseRecalibrationEngine.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/utils/recalibration/BaseRecalibrationEngine.java>)。

### 计算方式

每条 read：

1. consolidate CIGAR、clip adapter/soft clip、设置/恢复 base qualities；
2. 比较 reference 得到 `isSNP/isInsertion/isDeletion`；
3. 可选 BAQ；
4. 生成 read group、reported quality、cycle、context 等 covariates；
5. 对每个 offset/event type 更新 nested table：observations、errors、reported quality；
6. 全部 reads 完成后 collapse quality table 到 read-group table，并按 GATK 规则 round。

### CPU / SIMD / GPU

- **CPU**：reference/CIGAR/known-sites 判定、BAQ 状态机、hash/map 更新；
- **SIMD**：base/reference compare、quality→error probability、cycle/context 编码、error fraction；
- **GPU**：把 read batch 的 covariate matrix 和 `is*` arrays 生成放 device；每个 team 使用局部整数 histogram，再 deterministic merge；
- **ApplyBQSR**：读取、按 covariate 查表、改写 quality，通常 CPU streaming；只有读取已批量化且 GPU 队列有足够深度时才 offload。

### 坑点

- BQSR 表核心是整数 observation/error 计数，但 error fraction/quality correction 涉及 double；聚合顺序和最终 round 要固定；
- known sites、BAQ 开关、adapter clipping、original qualities 任何一个漏掉都会改变校准；
- covariate key 不是简单连续数组，必须保留 negative/invalid covariate 的 skip 语义；
- 每条 read 分配多个数组会造成 allocator/GC 压力，native 要使用 arena/reuse；
- Spark BQSR 的局部表合并顺序也要定义，否则与单机结果不一致。

## 10. MarkDuplicates、SortSam、外部内存排序

源码：[MarkDuplicatesSpark.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/spark/transforms/markduplicates/MarkDuplicatesSpark.java>)；非 Spark Picard 实现在外部依赖中。

### 计算方式

1. 以 read name 或 queryname group 收集 paired/secondary/supplementary/unmapped records；
2. 生成 duplicate key（library、unclipped 5' position、orientation 等）；
3. 在同一 key 中按 scoring strategy 选保留 read，其余打 duplicate flag；
4. 可选 optical duplicate tag、metrics；
5. coordinate sort 并写 BAM/CRAM/index。

### CPU / SIMD / GPU

- **CPU**：外部内存 run generation、多路 merge、read-name group 状态、duplicate semantics、metrics；
- **SIMD**：read name hash/compare、coordinate key extraction、CIGAR/reference position decode；
- **GPU**：只做一批固定 record keys 的 radix sort/partition；run 仍由 host 管理，最终 merge 和 stable order CPU；
- **I/O**：每个 run 有独立 BGZF writer，受 scratch byte/inode 预算控制。

### 坑点

- queryname-grouped 与 coordinate-sorted 输入性能和语义不同；
- stable sorting、secondary/supplementary/mate-unmapped 处理必须复现 Picard；
- optical duplicate metrics 使用 read name 中的 flowcell/tile/x/y 解析，不能当普通 hash；
- UMI 规则不等同于普通 duplicate marking；
- GPU sort 很快但磁盘 merge 仍是瓶颈，不能用 kernel benchmark 代替 end-to-end benchmark；
- spill 文件必须带 schema/version/CRC，作业中断后清理或恢复，不能留下错误的 partial BAM。

## 11. GenomicsDBImport / GenotypeGVCFs

源码：[GenomicsDBImport.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/genomicsdb/GenomicsDBImport.java>)、[GenotypeGVCFs.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/GenotypeGVCFs.java>)。

### 计算方式

- GenomicsDBImport：VCF/GVCF stream → interval partition → sample/locus transpose → sparse array fragments；
- GenotypeGVCFs：按 locus 读取所有 sample 的 reference block/variant → merge alleles → genotype likelihood/annotations → VCF。

### CPU / SIMD / GPU

- **CPU**：GenomicsDB/TileDB reader/writer、fragment、file descriptor、batch/partition、VCF header、locus iterator；
- **SIMD**：小规模 sample×allele likelihood、DP/AD/annotation numeric reductions；
- **GPU**：仅对足够大的同形状 locus batch 做 genotype math；GenomicsDB storage 不迁移到 GPU；
- **adapter**：保留 GenomicsDB backend，native 只负责 `batch-size`、VCF buffer、reader concurrency、tmp 和 memory telemetry。

### 坑点

- batch size 过大导致 reader、file descriptor、fragment memory 爆炸；过小会生成大量 fragment；
- GenomicsDB workspace 是持久状态，失败时不能随便覆盖或并发更新；
- sample order、contig partition、allele order 必须固定；
- GenotypeGVCFs 小矩阵通常受 branch/annotation/I/O 限制，GPU 可能比 CPU 慢；
- `gendb://` workspace 的 read-back 必须用原有 GenomicsDB/htsjdk/GATK 工具验证。

## 12. VQSR

源码：[VariantRecalibratorEngine.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/vqsr/VariantRecalibratorEngine.java>)、[GaussianMixtureModel.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/vqsr/GaussianMixtureModel.java>)。

### 计算方式

VQSR 是 VBEM/GMM：

1. 读取每 variant 的 annotation vector；
2. 初始化多个 Gaussian、mixture coefficients 和 covariance；
3. E-step：对每 datum、每 Gaussian 计算 log probability，log-normalize 得到 responsibility；
4. M-step：按 responsibility 累加均值、协方差、mixture weights；
5. 归一化并检查 mixture coefficient 收敛；
6. 对每 variant 计算 positive/negative model 的 log-odds，输出 VQSLOD/tranches。

### CPU / SIMD / GPU

- **CPU**：随机初始化、缺失 annotation 处理、协方差正定检查、收敛控制、tranche/VCF 输出；
- **SIMD**：variant×annotation 的中心化、平方/乘加、Gaussian log-density、worst annotation；
- **GPU**：按 Gaussian/annotation 维度做 batched matrix-vector 或小 GEMM；每轮输出 per-datum responsibility，host 做确定性归约；
- **内存**：annotation matrix 用 row-major/column-major 双视图，避免 GPU kernel 中 stride 访问。

### 坑点

- log10 normalization 必须用稳定 log-sum-exp，不能直接 `pow(10)` 后相加；
- covariance 接近奇异时的 regularization、逆矩阵和 NaN fallback 必须一致；
- GATK 的随机调用顺序影响模型，strict mode 要固定 RNG、数据顺序和 reduction；
- GPU 计算的协方差归约顺序可能改变 tranche 边界；
- 缺失 annotation 的随机 imputation、过滤和 contrastive evaluation 不能静默改语义。

## 13. CNV：read counts、SVD、segmentation、gCNV

源码：[CollectReadCounts.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/copynumber/CollectReadCounts.java>)、[DenoiseReadCounts.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/copynumber/DenoiseReadCounts.java>)、[ModelSegments.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/copynumber/ModelSegments.java>)、[SVDFactory.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/utils/svd/SVDFactory.java>)。

### 13.1 CollectReadCounts / CollectAllelicCounts

**计算**：读取经过 filter 的 read，取 read start，找到与 interval 的 overlap，给 interval count 加一；allelic counts 则在指定 heterozygous site 统计 ref/alt、base quality、strand 等。

- CPU：按 contig 排序的 interval cursor/overlap detector，避免每 read 扫全部 intervals；
- SIMD：batch read starts 与 interval start/end compare；base/quality/strand count；
- GPU：按 contig 分块，把 read starts sort/bin 后 histogram 到 intervals；不做全局 atomic 每 read 一次；
- 坑点：interval 必须不重叠或使用指定 merging rule；read start 与 read span 不是同一语义；HDF5/TSV 行顺序和 metadata 必须保持。

### 13.2 DenoiseReadCounts / SVD

**计算**：标准化 counts、GC bias correction，或对 panel-of-normals 做 SVD，保留指定 eigensamples，重构 denoised copy ratio。

- CPU：小 panel 使用 LAPACK/OjAlgo；保持 double 和稳定求解；
- SIMD：均值/方差/GC correction 的 row/column reductions；
- GPU：大 `samples × intervals` 矩阵用 Kokkos Kernels GEMM/SVD/QR；分块并把 U/V/singular values 计入 device budget；
- 坑点：SVD 的 singular vector 符号可翻转，不能用矩阵元素字节比较判断正确；奇异值截断 tolerance 要固定；HDF5 load/save 和 sample/interval order 不能变；小矩阵 GPU launch overhead 很高。

### 13.3 ModelSegments / gCNV

**计算**：kernel segmentation（copy ratio、allele fraction 或 combined kernel）、changepoint penalty、segment smoothing、MCMC posterior modeling；gCNV 另有 bias factors、copy-number states、posterior/ADVI 类模型。

- CPU：segment/changepoint 顺序控制、MCMC chain、state label/interval merge；
- SIMD：窗口 kernel 距离、log likelihood、posterior vector；
- GPU：大样本×interval 的 kernel matrix、batched state emission、posterior reductions；
- 坑点：segment boundary、contig/interval metadata、state label permutation、MCMC RNG 和 convergence；不能把每个 shard 的 segment 独立合并而忽略跨 shard boundary。

## 14. Funcotator、annotation、VCF/MAF/SEG

源码：[FuncotatorEngine.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/funcotator/FuncotatorEngine.java>)。

### 计算方式

按确定顺序遍历 data-source factories：variant/segment → reference context + feature context → transcript/gene/known-db lookup → funcotation map → VCF、MAF 或 SEG renderer。

### CPU / SIMD / GPU

- CPU：字符串、区间树/索引、transcript logic、factory dependency order、renderer；
- SIMD：只对固定宽度 numeric annotation（depth、quality、allele counts）做批量计算；
- GPU：默认不启用；只有数据源已预加载成 fixed-width table 且 lookup batch 很大时，才做 table scoring；
- I/O：memory-map/本地 SSD cache，LRU 受资源预算限制。

### 坑点

- factory 顺序是依赖关系，不可并行乱序；
- transcript/GENCODE 版本、reference build、B37↔HG19 contig 转换会改变结果；
- 缺失字段、默认值、annotation override 和 excluded fields 必须保持；
- 输出可能同时产生多个文件，不能只检查主 VCF。

## 15. SV、BWA、PathSeq、Spark 长尾

### SV

- CPU：evidence grouping、breakpoint clustering、contig graph、VCF assembly；
- SIMD：read evidence flag/count、k-mer/hash、depth/BAF numeric features；
- GPU：只加速大批量 evidence scoring 或 pairwise distance；graph clustering 先 host；
- 坑点：breakpoint evidence 的坐标、strand、mate orientation、split alignment 语义比 kernel 速度重要。

### BWA / read alignment

- CPU SIMD：BWA-MEM seed generation、FM-index/backtracking、Smith-Waterman extension；
- GPU：seed batch、DP extension、sorting/compaction；index 通常放 device/host pinned cache；
- 坑点：不同 BWA implementation 的 MAPQ、secondary/supplementary flag、CIGAR 和 tie-break 不一定等价；未完成前必须 fallback。

### PathSeq

- CPU/SIMD：k-mer encoding、taxonomic table lookup、quality filters；
- GPU：预加载 reference k-mer table 后做大批量 lookup；
- 坑点：hash table collision、taxonomy ordering、host/contaminant threshold 和输出排序。

### Spark tools

Spark 的 shuffle、partition、executor memory 和 spill 语义不是 Kokkos kernel 的职责。初期让 Nextflow/SLURM 进行 interval/sample 分片，Spark 工具走 Java fallback；只有 profile 确定某个 primitive（如 k-mer、PairHMM、sort key）值得移植时才替换 kernel。

### 15.1 其余 walkers 和顶层工具，不留“后面没方案”的空白

| 源码包/模块 | 计算方式 | CPU / SIMD | GPU 是否值得做 | 主要坑点 |
|---|---|---|---|---|
| `walkers/root`：CombineGVCFs、SelectVariants、ReblockGVCF、VariantFiltration、VariantAnnotator 等 | VCF/GVCF 顺序扫描、allele/block 合并、predicate/filter、header/annotation | CPU streaming；SIMD 做固定字段 compare/count | 默认不做；只有大批量 numeric annotation | allele order、header、END/block、过滤表达式和索引必须保持 |
| `consensus`：DownsampleByDuplicateSet | 按 duplicate set/coordinate 选择或下采样 reads | CPU grouping；SIMD 仅做 key compare | 不做 | duplicate set 定义、稳定选择和 read order |
| `variantutils`：LeftAlign/Trim、VariantsToTable、Split/merge | reference-aware normalization、allele trim、VCF field projection | CPU；SIMD 仅做 base compare/reverse-complement | 不做 | 左对齐的重复 indel、contig boundary、symbolic allele 不能近似 |
| `varianteval`：evaluators/stratifications | 按 variant/stratum 累计 TP/FP/FN、质量分布和统计量 | CPU map/reduce；SIMD counter/histogram | 大 cohort 可做 batch counter，但不优先 | stratification 顺序、missing field、浮点统计和报告格式 |
| `validation` | 交叉检查 VCF/BAM/reference、坐标、header、索引和统计 | CPU 可靠性优先；SIMD checksum/byte scan | 不做 | validation 不能因为 GPU 不支持就跳过；错误消息和退出码也属契约 |
| `coverage`：DepthOfCoverage、CallableLoci | locus/read coverage、callable threshold、区间统计 | CPU sorted sweep；SIMD depth/quality histogram | 高深度固定 bins 可做 GPU histogram | read overlap、interval merge、未覆盖区间、输出顺序 |
| `fasta`、`reference`、`GtfToBed` | FASTA 读写、reference base transform、GTF 字段解析/区间转换 | CPU I/O + SIMD reverse-complement/mask/ASCII classify | FASTA 很大且已分块时可做 base transform；通常 I/O-bound | line wrap、FAI/dict、1-based 坐标、GTF attribute 解析、reference build |
| `rnaseq`：ASEReadCounter、GeneExpressionEvaluation、SplitNCigarReads | N-cigar 分割、exon/allele count、strand/fragment 统计 | CPU CIGAR/interval；SIMD base/quality count | 默认不做 | splice junction、read pair/fragment、strand protocol、secondary reads |
| `featuremapping`：FlowFeatureMapper、FlowPairHMMAlignReadsToHaplotypes | flow-space base/flow 变换和 flow-aware PairHMM | CPU flow encoding；SIMD fixed flow vector | flow batch 足够大时 GPU，但必须独立 kernel | flow order、homopolymer collapse、flow quality 与 regular PairHMM 不兼容 |
| `groundtruth`：GroundTruthReadsBuilder/Scorer、AddFlowBaseQuality | truth read tagging、flow/base quality 派生、评分 | CPU；SIMD quality/base transform | 不优先 | truth tag、read name、原始 quality 和输出 header 不能丢 |
| `contamination`：Get/GatherPileupSummaries、CalculateContamination/Mixing | pileup allele counts、normal/tumor contamination likelihood | CPU pileup/interval；SIMD count/reduction | 大 cohort counts 可 batch GPU | site order、reference/alt orientation、zero-depth 和 prior |
| `readorientation`：CollectF1R2Counts、LearnReadOrientationModel | F1R2/F2R1 context histogram、orientation artifact model | CPU histogram + model control；SIMD context encode | 大 histogram 可 GPU reduce | context reverse-complement、sample/read-group key、tar/gzip sidecar |
| `realignmentfilter`、`filters` | alignment artifact/variant predicate | CPU branch-heavy；SIMD 只做无状态 mask | 不做 | filter 顺序和 short-circuit 会影响统计计数，不能任意重排 |
| `qc`：Pileup、CheckPileup、Flag/metrics、RSEM 后处理 | 读/位点统计、tag 转换、QC report | CPU streaming；SIMD numeric counters | 不优先 | metrics 字段、排序、read tags、外部脚本调用 |
| `variantrecalling`：HaplotypeBasedVariantRecaller | 依赖 haplotype/likelihood 的再召回和过滤 | CPU orchestration；复用 PairHMM/SW | 只复用已验证 kernel | region boundary、haplotype order、variant merge |
| `gnarlyGenotyper` | gVCF/variant 的特殊 genotyping 和合并 | CPU first；SIMD small matrices | 默认 fallback | 版本专用参数和输入格式，先建立独立 golden corpus |
| 顶层 `HtsgetReader`、`PrintReads`、`IndexFeatureFile`、`DumpTabixIndex`、`CRAMIssue8768Detector` 等 | range read、复制、索引生成/诊断 | CPU/I/O | 不做 | 这些命令的价值是兼容性和诊断，不能因“无 GPU”删掉 |
| `BwaMemIndexImageCreator`、`LocalAssembler`、`StructuralVariantDiscoverer` | index 构建、局部 assembly、SV evidence | CPU/SIMD 先 baseline | 只对独立 seed/k-mer/SW kernel | index byte layout、MAPQ/CIGAR、SV evidence coordinate |

这些工具多数不是算力瓶颈；正确策略是 CPU 高质量实现 + HTSlib/索引优化 + 明确 fallback，而不是强行 GPU 化。这样每个源码包都有可执行路径：`native CPU/SIMD`、`native GPU optional` 或 `Java fallback`。

## 16. 内存、I/O 和 kernel 选择器

每个模块都应提供同一个估算接口：

```text
estimate(work):
  host_bytes = decoded_records + metadata + output_buffer + cache
  device_bytes = packed_input + scratch + output + reduction_state
  scratch_bytes = spill_runs + index_ranges + restart_manifest
  arithmetic = estimated_cells / estimated_flops_per_second
  transfer = bytes_host_device / measured_link_bandwidth
```

选择规则：

```text
if unsupported or output contract not implemented:
    Java/Picard fallback
else if work < launch_threshold or transfer > compute:
    CPU scalar/SIMD
else if device_bytes > free_device * 0.6:
    smaller GPU batch or CPU
else:
    GPU team kernel
```

自适应 controller 要把以下全部计入预算：

- decoded read/haplotype/annotation arrays；
- Kokkos View、mirror view、pinned staging、device scratch；
- HTSlib index/cache、reference cache、GenomicsDB readers；
- output BGZF buffers、spill runs、sort merge heads；
- metrics、logs、restart manifest、allocator fragmentation。

## 17. 不能省略的验证矩阵

| 模块 | CPU oracle | SIMD/GPU 必测 | 输出级验证 |
|---|---|---|---|
| PairHMM | Java EXACT/ORIGINAL/GKL | likelihood matrix、empty/long read、flow/DRAGstr | HC/Mutect2 VCF/GVCF |
| SW | Java aligner/GKL | score、CIGAR、tie-break、overhang | assembly haplotype order |
| graph | GATK assembly result | k-mer count、path set、candidate order | HC calls/assembly debug |
| genotyping | Java GenotypingEngine | PL/GL、AF、allele permutation | genotype fields/annotations |
| BQSR | BaseRecalibrator tables | covariate matrix、table merge、round | recal table + recalibrated BAM |
| MarkDuplicates | Picard | duplicate flag/key/metrics | sorted BAM/index/metrics |
| GenomicsDB | official backend read-back | batch/fragment/tmp failure | workspace + GenotypeGVCFs |
| VQSR | GATK model/tranches | log density/responsibility/convergence | VQSLOD/tranches/VCF |
| CNV | Apache/OjAlgo/Spark path | SVD/posterior/segments | HDF5/TSV/VCF/segments |
| annotation | GATK Funcotator | numeric fields/order | VCF/MAF/SEG + sidecars |

没有通过这一矩阵时，native 实现只能标记为 `experimental`，不能被 `FAST_GATK_MODE=auto` 自动选中。

## 18. 实施优先级

```text
P0: I/O/engine contract + PairHMM + SW + oracle harness
P1: HC local assembly primitives + GVCF/reference confidence + Mutect2 reuse
P2: BQSR + genotype math + GenomicsDB adapter
P3: VQSR + CNV dense math + annotation cache
P4: external-memory MarkDuplicates/Sort + SV/BWA/PathSeq selective kernels
```

这里的 CPU/SIMD/GPU 不是三套不同算法，而是同一算法的不同执行策略；严格模式应优先保持数学和输出契约，只有 fast 模式才允许受控的浮点/归约差异。
