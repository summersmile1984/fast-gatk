# GATK 兼容的 Kokkos 计算库调研与实施方案

模块级 C++ Host、Kokkos kernel、依赖和验收拆分已整理到
[MODULE_IMPLEMENTATION_PLANS.md](</home/turing-agents/Documents/fast-gatk/MODULE_IMPLEMENTATION_PLANS.md>)。

## 结论先行

目标应定义为：

> 在不改 Nextflow 流程语义、不改 SLURM 资源管理方式、不改 GATK 文件/索引/参数契约的情况下，用一个 `gatk` 兼容入口替换 GATK 4.6.2.0；先由原生 Kokkos 实现接管高收益模块，未完成模块自动回退到 GATK。

不建议把 GATK 全部 Java 代码直接翻译成 Kokkos。GATK 的热点集中在少数数值内核；而读取、索引、VCF/BAM/CRAM 语义、图结构、字符串注释、外部内存排序和 GenomicsDB 存储属于不同类型的问题。最佳架构是：

```text
Nextflow process
      |
      v
SLURM allocation (cpus/memory/time/gpu/local scratch)
      |
      v
compatibility launcher: gatk
      |
      +--> native tool (C++/Kokkos/HTSlib/GenomicsDB adapter)
      |
      +--> Java GATK fallback (same arguments and same outputs)
```

第一阶段的“直接替换”应该是进程级替换，而不是跨 Nextflow 进程的隐式融合：每个 process 仍然产生 BAM/VCF/GVCF 及其索引，保证 `resume`、scatter/gather、缓存和失败重试不变。后续可增加显式的 fused 模式来减少中间文件，但不能把它作为兼容模式的隐含行为。

## 源码基线

本工作区已下载并固定 GATK 4.6.2.0 源码，当前提交为 `76edc75`（对应 4.6.2.0 tag）。源码目录：

- [gatk-source](</home/turing-agents/Documents/fast-gatk/gatk-source>)
- [GATK launcher](</home/turing-agents/Documents/fast-gatk/gatk-source/gatk>)
- [build.gradle](</home/turing-agents/Documents/fast-gatk/gatk-source/build.gradle>)
- [PairHMM.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/utils/pairhmm/PairHMM.java>)
- [AssemblyRegionWalker.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/engine/AssemblyRegionWalker.java>)
- [HaplotypeCaller assembly path](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/AssemblyBasedCallerUtils.java>)
- [GATK tests](</home/turing-agents/Documents/fast-gatk/gatk-source/src/test>)
- [完整模块清单与 native/fallback 状态](</home/turing-agents/Documents/fast-gatk/GATK_MODULE_INVENTORY.md>)
- [运行时和 Kokkos API 草案](</home/turing-agents/Documents/fast-gatk/FAST_GATK_RUNTIME_SPEC.md>)
- [CPU/SIMD/GPU 逐模块计算细节与坑点](</home/turing-agents/Documents/fast-gatk/GATK_CPU_SIMD_GPU_DETAIL.md>)
- [Rust-GATK（gatk-rs）源码调研与组合建议](</home/turing-agents/Documents/fast-gatk/RUST_GATK_RESEARCH.md>)
- [PairHMM native/SIMD 与 GATK 实测](</home/turing-agents/Documents/fast-gatk/PAIRHMM_SIMD_BENCHMARK.md>)

源码规模（Java 行数为粗略统计）：

| 区域 | 文件数 | 行数 | 说明 |
|---|---:|---:|---|
| `engine` | 131 | 18,375 | 遍历、分片、参考/reads/features 数据源 |
| `utils` | 401 | 78,005 | PairHMM、SW、BQSR、VCF/读取工具 |
| `tools/walkers` | 468 | 90,152 | 主体 variant walkers |
| `tools/spark` | 179 | 34,434 | 仅部分工具有 Spark 版本 |
| `tools/copynumber` | 122 | 18,260 | CNV 分析 |
| `tools/funcotator` | 62 | 19,008 | 注释与数据源查找 |
| `tools/genomicsdb` | 6 | 1,983 | GenomicsDB 导入/访问适配 |

源码中约有 211 个命令行程序注解类，但完整发行版还依赖外部 Picard、htsjdk、GenomicsDB、GKL、Spark/Hadoop 等。`build.gradle` 中的关键依赖包括 htsjdk 4.2.0、Picard 3.4.0、Spark 3.5.0、Hadoop 3.3.6、GenomicsDB 1.5.5 和 GKL 0.8.11。因此“全命令兼容”必须包含 Picard 命令的代理/回退，不应只移植 `src/main/java`。

GATK 4.6.2.0 目前是官方仓库标注的 latest release；官方仓库使用 Apache 2.0。实现时应保留 NOTICE、不要使用 Broad/GATK 商标暗示背书，并对复制的源文件做许可证审查。参考：[GATK releases](https://github.com/broadinstitute/gatk/releases)、[GATK repository](https://github.com/broadinstitute/gatk)、[GATK license](https://github.com/broadinstitute/gatk/blob/master/LICENSE.TXT)。

## 现有 GATK 数据流与可替换边界

### AssemblyRegionWalker 是天然的 work-unit 边界

[AssemblyRegionWalker.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/engine/AssemblyRegionWalker.java:38>) 明确把 reads 按 shard 读取，再划分 active/inactive assembly regions；`makeReadShards` 按 contig 分组并加 assembly padding（约 83–107 行），遍历阶段为每个 shard 设置过滤器、downsampler 和 transformer（约 157–174 行），最后对每个 assembly region 调用 `apply`（约 246–256 行）。

这给出 Kokkos 版本的基本契约：

1. Nextflow/SLURM 负责把样本和 interval list 分配到作业。
2. native engine 在作业内部按 contig/interval 建立 core interval + padded halo。
3. halo 只用于计算，输出只允许写 core 区域，避免 scatter 后重复记录。
4. 动态调度单位不能只按 bases 计数，应综合输入字节数、覆盖深度、active-region 密度、候选 haplotype 数和预计 PairHMM 工作量。

[Shard.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/engine/Shard.java:55>) 也明确规定 shard 的 core 区间和两侧 padding。native API 应直接保留这个概念，而不是把一条 contig 粗暴切成无 halo 的独立任务。

### GATK 已经有 I/O cache/prefetch，但不负责自适应预算

[GATKTool.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/engine/GATKTool.java:116>) 暴露 cloud prefetch buffer 和 BAM index caching；reads/features 初始化会把 prefetch buffer 传给数据源（约 458–507 行）。这是兼容层需要继承的参数语义，但 native 版还应加入：

- cgroup/SLURM 内存预算；
- local scratch 可用空间和 inode 预算；
- 输入格式、索引类型和远端对象存储延迟探测；
- 解码、计算、压缩、写出队列的背压；
- spill manifest、校验和、失败清理和 shard 级恢复。

### GATK launcher 行为必须保留

[gatk](</home/turing-agents/Documents/fast-gatk/gatk-source/gatk:1>) 不是简单的 `java -jar` 包装器：它处理 `--java-options`、`--gatk-config-file`、`--list`、`--dry-run`、Spark 参数分隔符 `--`，并为本地运行设置 htsjdk 压缩和异步 I/O 属性。官方 README 还明确提示直接绕过 launcher 会改变重要系统属性，包括压缩级别。因此兼容入口必须是同名 `gatk`，并保留这些选项。

## 模块源码分析与 Kokkos 方案

### 优先级矩阵

| 模块 | 源码事实/瓶颈 | Kokkos 实现建议 | 优先级 |
|---|---|---|---|
| PairHMM | 读 × haplotype 的动态规划；已有 Java、AVX、OpenMP AVX 多实现；可复用连续 haplotype 的缓存 | SoA ragged batch；按 read/haplotype 长度分桶；team/vector 并行；每个 team 使用 scratch；保留 logless 与精确参考模式 | P0 |
| Smith-Waterman | 单次 `align(ref, alt, params, overhang)`，需要保持 tie-breaking/overhang 语义 | 批量波前 DP；短序列走 CPU SIMD，较大批次走 CUDA/HIP/SYCL；结果回归逐字节/逐 CIGAR | P0 |
| HaplotypeCaller local assembly | De Bruijn/read-threading 图、k-mer 计数、错误纠正、图压缩、候选 haplotype 过滤；大量 Java map/set/JGraphT | 图拥有权放 host/NUMA arena；Kokkos 只做 k-mer 计数、排序/压缩、边评分、候选筛选；GPU 不作为第一版整图移植 | P0/P1 |
| Mutect2 | 继承 assembly-region；`Mutect2Engine.callRegion` 输出 VCF 和 stats，另有 somatic filters、F1R2、GVCF | 复用 HC 的 assembly/PairHMM/SW；肿瘤/正常逻辑和过滤先 CPU；保持 stats、F1R2、GVCF sidecar | P1 |
| GenotypeGVCFs | locus/variant 遍历、GVCF 合并、annotation、genotyping、VCF 输出；可访问 GenomicsDB | 先保留 GenomicsDB/TileDB backend；对每个 locus 的小矩阵做批量 CPU/GPU kernel；保留 allelic/annotation 顺序 | P1 |
| ReferenceConfidence/GVCF | reference block 合并、边界、非变异输出、header/索引严格 | CPU 流式 block combiner；PairHMM/likelihood 可加速，writer 不追求 GPU | P1 |
| BQSR | 每条 read 生成 covariate matrix，再更新 quality/read-group/额外 covariate 表；最后 collapse 和 round | SoA covariate + thread-local 整数计数表；Kokkos hierarchical reduction；ApplyBQSR 以流式 CPU/I/O 为主 | P1 |
| MarkDuplicates/Sort | 需要 read-name 分组、全局排序、内存和 spill；源码明确内存随 library complexity 增长 | 外部内存 sort：run generation、radix/key extraction、受控 spill；Kokkos 只做 key/比较/压缩；不能只做 GPU sort | P2 |
| GenomicsDBImport | sparse array；batch size 控制打开的 readers；大量 fragment 会放大文件句柄和内存；tmp 磁盘可能超过空间 | 先复用 GenomicsDB；自适应 batch/partition/VCF buffer/tmp；读取器并发受内存和 fd 限制 | P1 |
| VQSR | VBEM/GMM：`expectationStep`、`maximizationStep`、逐 annotation 评分 | Kokkos Kernels dense/batched GEMM、归约；模型控制/收敛/随机种子 CPU；大 cohort 按 feature block 分布 | P2 |
| CNV | CollectReadCounts 是 interval overlap + count；Denoise 是 SVD/矩阵；ModelSegments 是 kernel segmentation + MCMC | counts 用 SoA/并行 histogram；SVD/矩阵用 Kokkos Kernels；segmentation/MCMC CPU 首版，保留 HDF5/TSV 契约 | P2 |
| Funcotator | 多数据源按排序后的 factory 查找，字符串/区间/表格 I/O，输出 VCF/MAF/SEG | CPU cache、memory-map/索引批量查询；不应作为 GPU 优先项 | P3 |
| SV/PathSeq/BWA/Spark 长尾 | 算法和外部依赖多，部分命令只有 Spark 实现或依赖 Picard | 先 fallback；以后按 profiling 选择 BWA alignment、k-mer、SV evidence 等孤立 kernel | P3 |

### PairHMM：第一优先级

[PairHMM.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/utils/pairhmm/PairHMM.java:40>) 的 `Implementation` 枚举包含 EXACT、ORIGINAL、LOGLESS_CACHING、AVX、AVX+OpenMP 和 FASTEST_AVAILABLE。FASTEST_AVAILABLE 的回退顺序是 OpenMP AVX → AVX → Java logless，说明“按硬件选择实现”本来就是 GATK 的契约。

`computeLog10Likelihoods`（约 196–236 行）先找当前 batch 的最大 read/haplotype 长度，再按 `readCount * alleleCount` 写 likelihood matrix；`hapStartIndex` 用于连续 haplotype 的缓存。Kokkos 版不应把每个 read/haplotype 变成一次 kernel launch，而应：

```text
decode reads/haplotypes
  -> bucket by (read_length, hap_length, quality model)
  -> pack bases/qualities/indel penalties into SoA
  -> Kokkos team/vector DP with per-team scratch
  -> deterministic or fast reduction
  -> scatter likelihood matrix back to GATK-compatible order
```

必须独立实现 regular PairHMM 与 flow-based PairHMM；不能用普通四状态 DP 冒充 flow 模式。首版暴露 `--pair-hmm-implementation FASTEST_AVAILABLE` 的兼容参数，同时增加内部 `auto`：根据 Kokkos backend、序列长度分布、GPU 可用内存和 batch 规模选择 CPU/GPU。

### Smith-Waterman：第二个内核

Smith-Waterman 是接口干净、但结果细节敏感的内核。需要固定：

- SW 参数（match/mismatch/gap open/gap extend）；
- overhang 策略；
- 相同得分时的 tie-breaking；
- 输出 CIGAR 和 alignment start/end。

建议先建立 Java/GKL 对照 harness，再写 Kokkos CPU SIMD、CUDA、HIP、SYCL 三种 backend；通过同一批随机和真实 assembly-region 输入比较，而不是只比较总分。

### Local assembly/graph：host-first

[AssemblyBasedCallerUtils.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/AssemblyBasedCallerUtils.java:326>) 的 `assembleReads` 会构造 padded reference、可选 read error correction、SW 参数和 flow homopolymer collapsing，最终调用 `ReadThreadingAssembler.runLocalAssembly`（约 353–399 行）。后续 k-mer 计数和 haplotype 筛选大量使用 `HashMap`、`HashSet`、排序和字符串表示（约 430–486 行）。

因此第一版的 graph representation 应改成稳定的 flat arena/CSR-like edges，但所有权和复杂控制保留 host；Kokkos 只处理适合批量化的 primitive：k-mer 编码/计数、radix sort、unique/compact、edge score、候选排序。这样可以同时支持 CPU/OpenMP、CUDA、HIP、SYCL，而不把不规则图访问硬塞进 GPU。

### Mutect2 与 GenotypeGVCFs：复用而非重写

[Mutect2.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/mutect/Mutect2.java:220>) 继承 `AssemblyRegionWalker`，在 `apply` 中调用 `Mutect2Engine.callRegion` 并写 VCF（约 266–307 行），同时可能写 `.stats`、F1R2 和 somatic GVCF。最安全的 native 方案是复用 assembly/PairHMM/SW，先保留 Mutect2 engine 的过滤和 sidecar 语义。

[GenotypeGVCFs.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/GenotypeGVCFs.java:260>) 在启动时创建 annotation engine、GenomicsDB options 和 genotype engine，在每个 locus 调用 `gvcfEngine.callRegion`（约 288–330 行）。这是小矩阵、高分支、annotation 多的场景，应该先做 CPU batch/vector；GPU 仅在 profile 证明 genotyping math 占主导后加入。

### BQSR：非常适合 map/reduce

[BaseRecalibrationEngine.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/utils/recalibration/BaseRecalibrationEngine.java:115>) 对每条 read 计算 SNP/indel/BAQ/covariate，再在 `updateRecalTablesForRead` 中更新质量、read-group 和额外 covariate 表（约 255–288 行），最后 collapse 和 round（约 159–215 行）。源码注释还专门缓存 `EventType.values()`，说明对象分配在十亿级循环里是实际热点。

Kokkos 版用固定维度整数计数表、每线程/每 team 局部表和确定性 reduction；ApplyBQSR 维持顺序扫描和异步读写，只有在 batch 足够大时才把质量转换放到 device。分布式时先对局部表做整数合并，再按 GATK 的 round 规则输出，避免浮点归约顺序改变结果。

### MarkDuplicates/Sort：外部内存工程，不是单个 kernel

[MarkDuplicatesSpark.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/spark/transforms/markduplicates/MarkDuplicatesSpark.java:40>) 明确指出要把 read names 保存在内存；典型 30x WGS 至少建议 16 GB，内存随 library complexity 增长，过度 spill 会变慢（约 40–55 行）。源码先做 queryname sort，再只 shuffle 小的 read-name 元数据，避免 shuffle 大 read object（约 207–229、284–291 行）。

native 版应采用：

1. 受内存预算限制的 key/read run generation；
2. local SSD 上的压缩 spill runs；
3. 多路 merge，内存只保留 head buffers；
4. Kokkos 做 key extraction、radix sort 和批量比较；
5. 输出顺序、duplicate flag、metrics 和 sort order 严格匹配。

### GenomicsDBImport：保留后端，优化调度和资源控制

[GenomicsDBImport.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/genomicsdb/GenomicsDBImport.java:71>) 说明 GenomicsDB 是面向稀疏数组的样本转置存储；其参数文档指出 batch size 控制同时打开的 readers，数千 fragment 会造成约 20 倍文件打开和更高 bookkeeping 内存，tmp 磁盘需求可能超过默认目录（约 278–296 行及 177–183 行）。

因此不要在第一版从零替换 TileDB/GenomicsDB。native adapter 负责读取资源预算，动态决定 `batch-size`、VCF buffer、interval 并发、fragment consolidate 和 `--tmp-dir`，并把 SLURM 本地盘传给 GenomicsDB。

### VQSR、CNV、Funcotator

VQSR 的 [VariantRecalibratorEngine.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/vqsr/VariantRecalibratorEngine.java:34>) 是 VBEM/GMM：循环执行 expectation/maximization 直到收敛（约 103–128 行），适合 Kokkos Kernels 的批量 dense math，但模型控制和收敛判定先放 CPU。

CNV 的 [CollectReadCounts.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/copynumber/CollectReadCounts.java:176>) 是 read start 与 interval overlap 的计数，适合按 contig 的并行 histogram；[DenoiseReadCounts.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/copynumber/DenoiseReadCounts.java:184>) 使用 HDF5 读入和 SVD/矩阵运算；[ModelSegments.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/copynumber/ModelSegments.java:83>) 使用 kernel segmentation 和 MCMC，后两者先 CPU + Kokkos Kernels。

Funcotator 的 [FuncotatorEngine.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/main/java/org/broadinstitute/hellbender/tools/funcotator/FuncotatorEngine.java:227>) 按排序后的 data-source factory 逐个生成 funcotation，再由不同 renderer 输出 VCF/MAF/SEG（约 247–307 行）。它是字符串、区间查询和数据源 I/O 问题，不是 GPU 优先项。

## Kokkos 多硬件设计

Kokkos 提供 execution/memory abstraction，典型 backend 包括 OpenMP、CUDA、HIP、SYCL；Kokkos Kernels 提供 dense/batched/sparse/graph 数学。当前配置通常只能在一个构建中选择一个 device backend 和一个 host parallel backend，因此应交付多个构建产物，而不是假设一个二进制同时覆盖所有 GPU：

| 产物 | 目标 |
|---|---|
| `fast-gatk-cpu-x86` | x86 AVX2/AVX-512 + OpenMP，兼容无 GPU 节点 |
| `fast-gatk-cpu-arm` | ARM/SVE 或 OpenMP，避免 x86 指令假设 |
| `fast-gatk-cuda` | NVIDIA CUDA execution space |
| `fast-gatk-hip` | AMD ROCm |
| `fast-gatk-sycl` | Intel GPU/CPU SYCL |

运行时选择逻辑：

```text
FAST_GATK_MODE=auto|native|gatk
FAST_GATK_BACKEND=auto|openmp|cuda|hip|sycl

auto:
  if tool unsupported -> gatk fallback
  else if CUDA_VISIBLE_DEVICES/HIP/SYCL device valid -> native device path
  else -> native OpenMP path
  if correctness/format guard fails -> fail closed or explicit gatk fallback
```

Kokkos remote spaces或 MPI one-sided memory不应作为第一版的细粒度共享内存。跨节点的数据应以 interval/sample/shard 为单位批量传输；细粒度远端访问在高延迟网络上会抵消 kernel 收益。Nextflow+SLURM 已经提供了较好的任务级分布式边界。

## 自适应内存、设备内存和磁盘 I/O

### 预算模型

启动时读取：

- cgroup v2 `memory.max`/`memory.current`，回退到 `SLURM_MEM_PER_NODE`；
- `SLURM_CPUS_PER_TASK`；
- `CUDA_VISIBLE_DEVICES` 和 CUDA/HIP/SYCL 设备内存；
- `TMPDIR`/`SLURM_TMPDIR` 的 `statvfs` 可用字节、inode 和挂载类型；
- 输入 BAM/CRAM/VCF 大小、索引大小、是否远端 URI；
- reference 和 feature 数据源。

建议初始预算：

```text
job_mem = min(cgroup_limit, slurm_mem_limit_if_present)
reserve = max(1 GiB, 10% * job_mem)
compute_budget = job_mem - reserve
host_inflight_target = 0.60 * compute_budget
spill_limit = min(0.80 * scratch_free_bytes,
                  configured_spill_cap)
```

设备内存不直接等于 host budget；每个 backend 提供 `device_free_bytes()`，并以较小者决定 batch 大小。

### 三队列背压

把 I/O 和计算分开：

```text
read/decode threads -> bounded host queue -> Kokkos compute queue
       ^                                       |
       |                                       v
  reference/index cache <- compress/write threads <- result queue
```

Kokkos worker 只做计算，不阻塞网络/磁盘 I/O。队列元素必须带有 core interval、halo、输入偏移、输出序号和校验信息。

控制器每 250–1000 ms 采样 RSS、device free、scratch free、队列深度、read/write latency：

- RSS 或 device 使用率超过 80%，立即缩小 batch/并发；
- scratch 超过 75% 或 inode 低于阈值，停止扩大并发，必要时 spill/暂停 reader；
- 计算队列为空且 I/O 延迟低，按小步长增加 batch；
- I/O queue 满而 compute 空，增加预取/解码线程或压低压缩级别；
- compute queue 满而 I/O 空，减少 reader，避免 page cache 抢占；
- 任何分配失败都应退回较小 batch，而不是让作业被 OOM kill。

Batch 只能在安全边界（assembly region、read-name group、VCF block、GVCF locus block）缩放；不能在任意字节位置切断语义对象。

### 文件格式和 I/O 策略

- BAM/CRAM/VCF/BCF 继续用 htslib/HTS specs 兼容读写，不另造私有格式作为默认输出。
- BGZF 按 block 解压/压缩，配合 libdeflate；输入侧合并相邻 index ranges，输出侧多线程压缩。
- byte-compatible 模式固定 GATK/GKL 压缩等级、header 顺序和 writer 选项；fast 模式才允许根据 CPU/I/O 比例改变压缩等级。
- spill 使用版本化 run header、长度、CRC/校验和、输入 shard ID 和排序键；成功后按 manifest 清理，失败后可从最后一个完整 run 恢复。
- reference、FAI/dict、VCF index、annotation datasource 使用有上限的 LRU/cache；缓存大小计入同一 job budget。

## Nextflow + SLURM 兼容策略

Nextflow 继续负责 `cpus`、`memory`、`time`、container、queue、retry、resume、scatter/gather；SLURM 负责 allocation 和 GPU GRES。fast-gatk 默认绝不自行调用 `sbatch`，也不在一个 SLURM task 内偷偷启动第二套 MPI/Spark 集群。

典型配置只需增加 label/容器和 GPU 资源（GPU 版本）：

```groovy
process.executor = 'slurm'

process {
  withLabel: fast_gatk {
    cpus = 16
    memory = '64 GB'
    time = '24h'
    container = 'registry.example/fast-gatk:cuda'
    clusterOptions = '--gres=gpu:1'
  }
}
```

native binary 从 `SLURM_CPUS_PER_TASK` 读取线程数，从 cgroup/`SLURM_MEM_PER_NODE` 读取上限，从 `SLURM_TMPDIR` 选择 spill；GPU 只在 SLURM 已分配 GRES 时启用。CPU 节点无需修改命令行即可走 OpenMP/fallback。

需要明确两种兼容边界：

1. `gatk Tool ...`：同名入口可完全截获。
2. `java -jar gatk-package...`：绕过 launcher，无法由 PATH shim 截获；需要容器 entrypoint、模块文件或把流程命令改为 `gatk Tool ...`。

## `gatk` 兼容层设计

实现一个 Python/C++ 小型 dispatcher，安装到模块环境中：

```text
gatk Tool args
  -> parse launcher-only options
  -> preserve @args, --, --java-options, --gatk-config-file
  -> classify tool
  -> native registry lookup
  -> resource probe and backend selection
  -> native execution or Java fallback
  -> validate output/index/sidecar contract
```

必须支持：`gatk --list`、`Tool --help`、`Tool --version`、`--dry-run`、`--java-options`、`--gatk-config-file`、argument file `@file`、Spark 参数分隔符 `--`。兼容解析器不应把未知工具参数吞掉。

Java 选项映射：

- `-Xmx`：转为 native host budget 上限；不能超过 SLURM/cgroup；
- `-Djava.io.tmpdir`：作为 spill/tmp 优先目录；
- `-Dsamjdk.compression_level`：byte-compatible 模式固定 writer 参数；
- `-XX:*`：native 路径记录并忽略 JVM 专用项，fallback 原样传递；
- Spark-only arguments：native 路径默认拒绝误用并给出清晰提示，Spark 工具先 fallback。

输出契约至少包括：文件名、索引后缀、VCF/BAM header、contig/order、坐标排序、duplicate flags、GVCF block、stats/F1R2/metrics sidecar、退出码和错误类型。不要把“语义相同但 header 不同”误报为兼容完成。

## 验证与性能门槛

### 四层兼容性

1. CLI：参数、默认值、`--help`、`--dry-run`、错误码。
2. 格式：SAM/BAM/CRAM/VCF/BCF、索引和压缩可被 htslib/GATK 读取。
3. 语义：记录数、排序、header、过滤、GVCF block、metrics/sidecar。
4. 生物学结果：variant/haplotype/genotype 与 reference truth 的一致性。

字节级一致不是所有工具的合理目标；应定义 `deterministic` 模式（固定 reduction、排序、随机种子和压缩）与 `fast` 模式（允许浮点归约顺序变化，但结果在规定 tolerance 内）。

### 现有 GATK 测试可以直接变成 oracle

仓库包含 [PairHMMUnitTest.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/test/java/org/broadinstitute/hellbender/utils/pairhmm/PairHMMUnitTest.java>)、[VectorPairHMMUnitTest.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/test/java/org/broadinstitute/hellbender/utils/pairhmm/VectorPairHMMUnitTest.java>)、[HaplotypeCallerIntegrationTest.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/test/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/HaplotypeCallerIntegrationTest.java>) 和 [Mutect2IntegrationTest.java](</home/turing-agents/Documents/fast-gatk/gatk-source/src/test/java/org/broadinstitute/hellbender/tools/walkers/mutect/Mutect2IntegrationTest.java>)。建议：

- PairHMM：先逐矩阵比较 Java EXACT/ORIGINAL/LOGLESS/GKL，再比较 Kokkos；
- SW：随机序列 + 真实 assembly regions，比较 score、CIGAR、overhang 和 tie-breaking；
- HC/Mutect2：同一输入、同一 reference、同一 interval、同一压缩等级，比较 VCF/GVCF/stats；
- BQSR/MarkDuplicates/GenomicsDB：比较 table、metrics、workspace 可读性和重启行为；
- 全流程：用 hap.py/GA4GH/NIST GIAB truth/confident regions 做 haplotype-aware 评估，不能只看 VCF diff。

### 建议的门槛（目标，不是未经测试的承诺）

| 阶段 | 门槛 |
|---|---|
| kernel | PairHMM/SW 在各 backend 通过 oracle；CPU 至少不劣于 GKL 基线；GPU batch 达到明确吞吐目标 |
| single-node HC | native 结果与 GATK 语义一致；峰值 RSS < 分配内存的 70%；spill 可控 |
| cluster | interval scatter 后输出可 gather；在 I/O 不成为瓶颈前保持可观的强/弱扩展 |
| production | 失败重试不损坏输出；SLURM OOM/磁盘不足转为可读错误；fallback 可用 |

基线必须是 GATK 4.6.2.0 + GKL `FASTEST_AVAILABLE`，不能用纯 Java 慢实现来夸大加速比；固定 reference、input、interval、压缩等级、线程/内存和容器。性能报告至少拆成 decode、assembly、PairHMM、SW、genotyping、compression、write、spill 八项。

## 分阶段路线图

### P0：可替换但不改变结果

- `gatk` dispatcher 和 `FAST_GATK_MODE=auto|native|gatk`。
- native tool registry、资源探测、fallback、输出验证。
- HTSlib/BGZF/libdeflate 读写层，BAM/CRAM/VCF/BCF/index smoke tests。
- 容器/Environment Modules：CPU、CUDA、HIP、SYCL 版本。
- 用 GATK 测试资源建立 golden corpus。

### P1：单节点 HaplotypeCaller/Mutect2 热点

- Kokkos PairHMM regular mode。
- Kokkos Smith-Waterman。
- host-first local assembly + Kokkos k-mer primitives。
- AssemblyRegion core/halo、downsampling、reference/feature cache。
- 自适应 host/device batch、scratch、prefetch、spill controller。

### P2：完整 germline/somatic 直接替换

- HC `-ERC GVCF`、ReferenceConfidence、GVCF writer。
- Mutect2 VCF/stats/F1R2/somatic GVCF。
- GenotypeGVCFs 的 GenomicsDB adapter 和 batch genotype math。
- Nextflow+SLURM 端到端 scatter/gather、retry/resume。

### P3：BQSR、GenomicsDB、CNV/VQSR

- BQSR map/reduce table 与 ApplyBQSR streaming。
- GenomicsDBImport 自适应 batch/tmp/fragment 管理。
- VQSR dense math、CNV SVD/counts；Funcotator CPU cache。

### P4：外部内存和长尾

- MarkDuplicates/Sort external-memory engine。
- SV、BWA、PathSeq、Spark 工具的选择性 native kernels。
- fused mode（显式 pipeline profile），减少中间 BAM/VCF，但不影响 drop-in mode。

## 风险与决策

- **最容易成功的切入点**：PairHMM + SW + HTSlib I/O + dispatcher。
- **最容易失败的切入点**：从零重写 GenomicsDB、把整个 assembly graph 直接搬到 GPU、把 MarkDuplicates 当成普通 GPU sort、用一个 Kokkos binary 假设覆盖所有厂商 GPU。
- **兼容性的真正难点**：Picard 外部命令、header/索引/sidecar、GVCF block 边界、浮点 reduction、Java launcher 行为、Spark 工具的集群语义。
- **资源管理原则**：Nextflow/SLURM 管作业级资源；fast-gatk 只在 allocation 内部自适应，不自行 `sbatch`，不突破 cgroup，不把 I/O 阻塞放入 Kokkos worker。

最终建议把项目名称和对外命令定为“GATK-compatible”，而不是宣称官方 GATK 或替代所有 GATK 实现；在每个工具的兼容矩阵中明确 `native / fallback / unsupported` 三种状态。

## Rust-GATK 吸收边界

`gatk-rs` 只作为兼容性参考和验证资产吸收，不作为我们的运行时依赖。直接纳入：固定 reference/dependency/toolchain provenance、Java oracle/golden、t-wise 参数覆盖、差分 fuzzing、determinism gate、错误/边界语义、interval/pileup/read-filter 的测试向量，以及 byte-identical/bio-identical 分级。选择性借鉴其 Rust host 数据模型和 `@PG`/header/索引验证；不复制其未完成 caller、library-only 结构、无 CLI 形态或“Rust 天然更快”的假设。详见 [Rust-GATK 吸收清单](</home/turing-agents/Documents/fast-gatk/RUST_GATK_RESEARCH.md#11-只吸收有价值的部分>)。
