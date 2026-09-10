# fast-gatk 全模块执行方案（草案）

> 本文件是可执行的工程任务单，不替代总体架构规格。它把
> [GATK_MODULE_INVENTORY.md](</home/turing-agents/Documents/fast-gatk/GATK_MODULE_INVENTORY.md>)
> 和
> [MODULE_IMPLEMENTATION_PLANS.md](</home/turing-agents/Documents/fast-gatk/MODULE_IMPLEMENTATION_PLANS.md>)
> 展开为模块边界、输入输出、Host/Kokkos 分工、依赖、完成定义、oracle、性能门禁和风险。
> 现阶段文件名为 `DRAFT`，在 P0 设计评审通过后冻结为正式执行基线。

> **历史快照提示（2026-08-31）**：本草案保留了早期实施边界，不能直接作为当前
> 完成状态。例如 §1.1 中 HTSlib 的 CIGAR/index、HC 的 assembly/RCM 和集群能力
> 仍以“待完成”描述；这些条目后来已有不同程度的实现，但尚未因此达到 GATK
> 全算法/bit-identical 或生产替换门槛。当前可复核状态请以
> [FAST_GATK_EXECUTION_PLAN.md](</home/turing-agents/Documents/fast-gatk/FAST_GATK_EXECUTION_PLAN.md>)
> 顶部的最新增量和
> [PROGRESS_AUDIT_2026-08-31.md](</home/turing-agents/Documents/fast-gatk/PROGRESS_AUDIT_2026-08-31.md>)
> 为准；本文件仍作为“所有模块的执行任务分解”使用。

## 0. 目标、范围和不可妥协项

目标是保持 GATK + SLURM + Nextflow 的调用和产物契约，同时把适合批处理的数值阶段
改为 C++ Host + Kokkos API 的实现。`C++ Host` 负责对象语义、文件、状态机、排序、
异常和输出；`Kokkos kernel` 只接收平面化、定长或 offset 描述的 batch。所有后端共享
同一份 C++/Kokkos 源码，不能为 AVX、CUDA、HIP、SYCL 分叉算法语义。

本计划覆盖清单中的所有源码区域、utils 子包和工具包。`native` 不等于“所有代码都
搬到设备”：字符串、索引、复杂图遍历、压缩、GVCF block 状态和 writer 明确留在 Host。
不满足 oracle 或输出契约的工具只能登记为 `fallback`，不能宣称可替换。

必须满足的硬约束：

- Strict 模式固定输入顺序、排序 tie-break、归约树、随机种子和数学运算顺序；可达到
  的阶段要求 raw-bit 或逐 CIGAR 一致，不能用 fp32/fp64 宽松比较掩盖实现差异。
- Fast 模式允许 backend 优化，但要有工具级 tolerance 和生物学结果门槛；manifest
  必须记录模式，不能把 Fast 输出标成 Strict。
- `HostBatch -> DeviceBatch -> KernelPlan` 是统一生命周期边界。prepare（packing、
  allocation、deep-copy）和 execute 必须独立计时，Plan 在 region/batch 生命周期内
  复用。
- 所有写出工具都必须先写临时文件，校验 header、排序、index、sidecar、CRC 和
  `OutputManifest`，再原子 rename；任何缺件都返回非零错误码。
- 设备 kernel 不调用 HTSlib、文件系统、网络、SLURM、Java，也不保存 Host 对象图指针。
- 输入、输出和错误分类必须可序列化，支持 SLURM cgroup 限制、scratch 背压、失败
  重试和 Nextflow resume。

## 1. 现状证据与基线

### 1.1 已完成且可复用的资产

| 资产 | 当前路径 | 当前证据 | 仍缺的生产能力 |
|---|---|---|---|
| 统一 batch/plan | `fastgatk-core/include/fastgatk/core/{batch,plan}.hpp` | PairHMM 和 SW/BQSR/assembly/region 原型已链接同一 library | persistent workspace pool、schema/version、allocator 和多阶段 plan |
| PairHMM production kernel | `fastgatk-kernels/src/pairhmm_kokkos.cpp` | scalar/AVX2/AVX512 4096 结果 bit-different=0；使用 `Kokkos::Experimental::simd<double>`；demo 仅为客户端/oracle | ragged length、regular/flow 分离、GATK full matrix 已接入；设备 wavefront 与完整 GKL flow corpus 仍待完成 |
| 四类 Kokkos production primitives | `fastgatk-kernels/src/{smith_waterman_kokkos,kmer_graph,bqsr,activity_profile}.cpp` | SW、BQSR、k-mer、region 已纳入共享 kernel target，并通过跨 backend signature/API smoke | 更完整的真实 read/CIGAR/graph corpus、稳定内存和 tool-level composition |
| HTSlib reader | `fastgatk-core/src/hts_reader.cpp` | BAM 全量、BAM region、CRAM+FASTA 可读 | CIGAR-aware projection、index seek、VCF/BCF writer、云 range |
| HC smoke | `fastgatk-native/src/hc_smoke.cpp` | batch reader + Kokkos reductions，manifest/telemetry | 完整 assembly、SW、PairHMM、genotype、VCF 等价 |
| HC call prototype | `fastgatk-native/src/calling_pipeline.cpp`、`hc_call.cpp` | Kokkos pileup/CIGAR indel、bounded local haplotype combinations、marginalized PairHMM likelihood、genotyping；多 ALT VCF/GVCF、完整 indel REF-span block 排除、跨 batch/thread signature 稳定 | GATK-equivalent pruning/traversal、full prior/reference-confidence semantics、GATK oracle bit identity |
| 兼容/集群 smoke | `fastgatk-native/scripts/gatk_compat.sh`、`workflow/*.sh`、`*.nf` | 参数别名、Nextflow process 等价路径、SLURM 本地资源模拟 | 真实 SLURM allocation、registry、fallback 和多 shard gather |

### 1.2 必须保留的性能基线

- PairHMM matrix8 的 Kokkos steady-state 为 GKL 实际矩阵吞吐的约 1.09--1.28 倍；将
  一次 prepare 摊入十次迭代后，16 核约为 0.94 倍。所有后续优化报告必须同时给出
  kernel-only、prepare、amortized 和端到端，不得只报告热 kernel。
- SW demo 从约 355M 到 4.94G cell updates/s（1 到 16 核）是规则 DP 的 CPU 参考；
  新实现不得因 kernel 化引入超过 10% 的单节点回归。
- BQSR 全局 bin atomic 和 k-mer 图高核数退化已被观察；默认采用 per-team/local
  table、稳定 sort 和 Host graph ownership。
- 小 region 计算受 I/O 和启动成本限制；batch 以 bytes 背压，低于阈值走 Host fast path。

### 1.3 现阶段不应误标的能力

当前 `fastgatk-hc-call` 已支持有界多 ALT 位点、局部 haplotype combination 和
多等位 VCF/GVCF AD/PL 输出；其 call-set 与 GATK 在 chr17 fixture 上只证明存在
可比较 sentinel/交集，不能称为 HaplotypeCaller 等价。P0/P1 的完成标准仍需
补齐 GATK pruning/traversal、全量 prior/reference-confidence 和 bit-level oracle。

## 2. 统一实现契约

### 2.1 每个模块必须交付的五个对象

```cpp
struct HostBatch;                         // std::vector/arena 所有权
struct DeviceBatch<ExecSpace>;           // typed Kokkos::View 和 byte accounting
template<class ExecSpace> class KernelPlan;
PrepareResult prepare(KernelPlan&, const HostBatch&, ExecSpace);
ComputeResult execute(KernelPlan&, const DeviceBatch&, ExecSpace, DeterminismMode);
CollectResult collect(const ComputeResult&, HostBatch&, DeterminismMode);
```

实际工程可把 `HostBatch`/`DeviceBatch` 按模块扩展为强类型结构，但不能绕过 core
plan 的生命周期和 telemetry。`prepare` 禁止在每一轮 execute 中重新分配大 View。

### 2.2 统一数据 schema

| schema | Host 输入布局 | Device 布局 | Host 输出 |
|---|---|---|---|
| `ReadBatch/v1` | `offsets:uint32[n+1]`、`bases:uint8[]`、`qualities:uint8[]`、`tid:int32[n]`、`pos:int32[n]`、`mapq:uint8[n]`、CIGAR/flags/read-group side arrays | SoA Views；变长字段只以 offset/length 传递 | read order、filter reason、统计和可重放 chunk |
| `AlignmentBatch/v1` | reference/read offsets、scoring params、stable id | bases/score/trace checkpoint Views | score、start/end、CIGAR opcode/length，CIGAR 由 Host compaction |
| `LocusBlock/v1` | coordinate key、observation offsets、base/quality/mapQ | compact observation、histogram、activity/coverage Views | stable coordinate order、pileup、candidate loci |
| `KmerGraph/v1` | canonical key、node/edge arena、read path offset | encoded key、sort/unique、support、edge CSR primitive | graph ownership、prune、path、haplotype strings |
| `PairHmmBatch/v1` | pair ids、read/haplotype offsets、quality/transition tables | bucketed/ragged Views、3-row or tiled workspace | likelihood matrix，恢复原 pair 顺序 |
| `VariantBlock/v1` | allele/sample order、likelihoods、annotations、GVCF state | fixed small likelihood/genotype matrices | VCF/GVCF records、header/index、sidecar |
| `ModelBlock/v1` | dense matrix、feature schema、seed/model metadata | block matrix / Kokkos Kernels views | model parameters、report、convergence/segments |

每个 schema 增加 `schema_version`、`input_ordinal`、`output_ordinal`、`contig_id`、
`core_interval`、`halo_interval`、`determinism_mode` 和 byte budget。未知版本必须 fail
closed，不能静默按旧布局解码。

### 2.3 统一结果/错误/遥测

错误枚举固定为 `BAD_INPUT`、`RESOURCE_EXHAUSTED`、`BACKEND_UNAVAILABLE`、
`NUMERICAL_CONTRACT_FAILURE`、`OUTPUT_CONTRACT_FAILURE`、`FALLBACK_EXECUTED`。
每一阶段记录 `prepare_bytes/device_bytes/prepare_seconds/execute_seconds/bytes_read/`
`bytes_written/spill_bytes/queue_depth/backend`。manifest 的 `status=complete` 只有
primary output 和所有 sidecar 都可回读时才允许。

## 3. 执行顺序、任务 ID 和依赖图

任务状态只能取 `planned -> implementing -> verified -> registered`。`verified` 需要
对应的 oracle、格式、资源和性能测试；`registered` 才能被 launcher 自动选择。

```text
FG-001 contract/registry
  -> FG-010 runtime/resource + FG-020 HTSlib/model + FG-030 output
  -> FG-040 traversal/queue
  -> FG-100 PairHMM + FG-110 SW + FG-120 k-mer/assembly primitives
  -> FG-130 pileup/activity + FG-140 genotype/reference confidence
  -> FG-200 HaplotypeCaller
  -> FG-210 Mutect2 / FG-220 BQSR
  -> FG-230 GenomicsDB/GenotypeGVCFs
  -> FG-240 MarkDuplicates/Sort + FG-250 VQSR/CNV
  -> FG-300 annotations/utility + FG-310 SV/BWA/PathSeq long tail
  -> FG-400 GPU backends + FG-410 real SLURM/Nextflow + registry promotion
```

### 3.1 P0：共用底座和第一个可组合工具

| ID | 可执行任务 | 依赖 | 交付物/完成定义 |
|---|---|---|---|
| FG-001 | 把当前 registry 结构化为 JSON schema；记录参数白名单、输出/sidecar、backend、strict/fast、golden corpus、资源估算和 fallback 命令 | 无 | unknown 参数不静默忽略；`--dry-run` 输出选择理由；单元测试覆盖 native/fallback/fail-closed |
| FG-010 | 实现 `ResourceProbe`/`ResourceBudget`/`AdaptiveController`；读取 cgroup v2、SLURM CPU/memory/GRES、scratch/inode/fd | 001 | 在受限 cgroup 内不会按物理机内存分配；压力时降低 batch 而非 OOM 无限重试 |
| FG-020 | 扩展 HTSlib reader 为 index seek、CIGAR decode、flags/read group、reference cache、VCF/BCF/BGZF/CRAM writer；保留 `ReadBatch/v1` | 001 | BAM/CRAM/VCF/BCF/samtools read-back、header/index/PG 一致；压缩失败不产生假成功 |
| FG-030 | 建立 `OutputManifest`、atomic writer、sidecar schema、partial cleanup/checkpoint | 001 | 删除/中断/磁盘满场景可恢复；所有生产工具返回前 manifest complete；runtime 已提供可复用 bundle validator/publisher，覆盖缺 index、临时残留和失败回滚 |
| FG-040 | 实现 core/halo interval、deterministic shard order、byte-bounded queues、reader/compute/writer CPU pool | 010/020/030 | scatter/gather 无重复漏计；队列按 bytes 背压；单 shard 可 replay |
| FG-100 | 提炼 `PairHmmPlan`：length bucket、persistent workspace、quality/transition cache、regular/flow API、collect permutation | 020/040 | GATK/GKL matrix strict oracle；prepare amortized 不低于当前基线；长短/空/异常输入全覆盖 |
| FG-110 | 提炼 SW score + Host traceback；输出 checkpoint、CIGAR、overhang/tie-break | 020/040 | score、CIGAR、start/end 和边界与 GATK/GKL 一致；不允许只比较 score |
| FG-120 | k-mer encode/sort/unique/count/edge compact primitive；Host graph ownership/traversal | 020/040 | GATK ReadThreadingGraph/SeqGraph corpus 的 edge/path/haplotype 顺序一致 |
| FG-200 | HC fixed-parameter native path：traversal → activity → assembly → SW → PairHMM → genotype → VCF | 030/040/100/110/120 | HC integration golden、VCF/header/index、GVCF smoke、fallback unknown options；才可标 `native-candidate` |

### 3.2 P1：生产 HC、Mutect2、BQSR

| ID | 可执行任务 | 依赖 | 完成定义 |
|---|---|---|---|
| FG-130 | activity profile、locus iterator、downsampling、BAQ/pileup；Host 保持 filter 顺序和 seed，Kokkos 处理窗口 histogram/quality primitive | 020/040 | 固定 seed、window/halo、filter/count 与 GATK 一致；低覆盖窗口 Host fast path |
| FG-140 | allele/read likelihood transform、log-sum-exp、genotype likelihood、reference confidence、GVCF block state | 100/110/130 | AD/DP/GQ/PL、block 边界、underflow/NaN、header/index 对照通过 |
| FG-200A | 完整 CIGAR-aware assembly：read threading、indel haplotype、prune、reference-connected path | 120/130 | indel corpus、homopolymer、complex graph、haplotype stable order |
| FG-200B | HC orchestration、AssemblyRegion/halo/downsample、region-level PipelinePlan | 200A/140 | GATK 4.6.2 HC test resources/real truth；同样输入可 resume/gather |
| FG-210 | Mutect2 复用 HC；Host 维护 tumor/normal、orientation、contamination、filters、F1R2/stats/somatic GVCF | 200B/140 | Mutect2 integration 与所有 sidecar；native/fallback 可解释 |
| FG-220 | BQSR covariates/table/merge；ApplyBQSR 首版 Host stream，批量 transform 可 Kokkos | 020/040/130 | report/table、read qualities、order/header、分片 merge/restart 和 GATK 对照 |

### 3.3 P2：联合分型、外部内存和 dense math

| ID | 可执行任务 | 依赖 | 完成定义 |
|---|---|---|---|
| FG-230 | GenomicsDB/TileDB adapter：batch、reader、VCF buffer、fragment/consolidation、tmp/fd budget；GenotypeGVCFs locus block math | 030/140 | 原 GATK 可读取 workspace；sample/allele order、数千 sample 资源压测 |
| FG-240 | MarkDuplicates/SortSam external-memory runs、spill manifest、stable multiway merge；Kokkos 只抽取/排序 key | 020/030/040 | duplicate flags/metrics、sort order、index、memory cap、spill/restart |
| FG-250 | VQSR VBEM/GMM、CNV counts/SVD/GEMM；HDF5/TSV/model state Host | 140/230 | VQSR tranche/report、CNV HDF5/segments、SVD tolerance、模型 restart |
| FG-260 | root walkers（SelectVariants、CombineGVCFs、ReblockGVCF、GenotypeGVCFs 等）接入统一 VariantBlock/writer | 030/140/230 | VCF/GVCF fields/header/index/ordering 全部 read-back；不支持模式 fallback |

### 3.4 P3：注释、SV、Spark/Picard 长尾和 GPU

| ID | 可执行任务 | 依赖 | 完成定义 |
|---|---|---|---|
| FG-300 | VariantUtils、VariantEval、annotation、Funcotator；Host 负责字符串/data-source/plugin，只有纯数值字段批处理 | 030/260 | VCF/MAF/SEG/header/missing value/plugin 输出一致；profile 证明收益才 native |
| FG-310 | SV evidence count/cluster/stratify、BWA seed/extend、PathSeq k-mer classify；一项一项注册 | 120/030 | 每个 primitive 有独立 oracle/benchmark；没有 oracle 继续 fallback |
| FG-320 | Spark/Hadoop、Picard、unsupported、report/text/config/python/R 保持 fallback/adaptor；只提供 native launcher、日志和 resource mapping | 001/010/410 | 原命令可执行、exit code/日志/输出契约透明；不在 runtime 复刻 Spark DAG |
| FG-400 | Kokkos CUDA/HIP/SYCL backend：先 correctness scalar-equivalent，再 team/vector/tile 优化 PairHMM/SW/k-mer/genotype | 100/110/120/140 | CPU/GPU strict corpus；显式 scratch/device OOM；跨 backend signature/format oracle |
| FG-410 | 真实 SLURM allocation、Nextflow scatter/gather/resume/retry/cache、容器和多节点 I/O | 200/210/220/230/400 | 真实 scheduler job 中资源、telemetry、sidecar、failure/retry 通过；才允许默认 native |

## 4. 平台和基础模块实施明细

下表对 inventory 中非 tool 的所有区域给出具体边界。`D0` 指必须先通过 P0 公共门禁；
`D1` 指依赖对应数值 kernel；`D2` 指可在长尾阶段执行。

| 模块/任务 | Host-only 边界 | Kokkos kernel 与 batch | 输入/输出 | 依赖、oracle、性能门禁和风险 |
|---|---|---|---|---|
| `cmdline/compat` (`FG-001`) | `gatk` dispatcher、`@args`、`--java-options`、registry、fallback、错误码、原始 argv | 无；kernel 不解析参数 | argv/registry → selection/manifest/exit code | D0；GATK help/version/error/Nextflow 命令 oracle；启动额外开销 < 1% wall；风险是误吞参数，unknown 必须 fail closed |
| `engine/traversal` (`FG-040`) | shard/core/halo、contig order、AssemblyRegion state、downsample seed、read order | fixed window activity/filter/count | interval/read chunks → ordered `LocusBlock` | D0；GATK traversal + scatter/gather oracle；无重复/漏计；风险是 halo double output |
| `exceptions` | 分类、上下文、cause chain、可读日志、fallback reason | 无 | internal error → stable code/manifest | D0；故障注入 oracle；错误分类不混淆；风险是异常跨 backend 泄漏 |
| `metrics` (`FG-030`) | telemetry aggregator、RSS/device/scratch/I/O、JSON schema | optional device counters/reduction | stage events → telemetry sidecar | D0；重复运行 schema 可解析；计时误差 < 3%；风险是只记录 kernel 造成错误结论 |
| `runtime/io/nio/gcs` (`FG-010/020/040`) | HTSlib/BGZF/CRAM、reference/cloud range、prefetch、spill、bounded queues、cgroup/scratch/fd | decode 后 filter/transform，不能传 file pointer | BAM/CRAM/VCF/BCF/cloud → schema batches/files | D0；samtools/GATK read-back、fault/retry；I/O bandwidth 不低于 baseline 90%；风险是缓存超 budget |
| `read/codecs/variant/fasta/reference` (`FG-020`) | CIGAR、header/PG、VCF semantics、reference cache、string/codec | base/quality transform、CIGAR op classification、fixed annotation primitive | flat records ↔ HTS files | D0；HTSlib/GATK byte/field oracle；禁止无索引全文件扫描替代 region seek |
| `transformers` (`FG-130`) | filter order、read group、soft clip、duplicate/adapter semantics | stateless mask/quality transform | ReadBatch → filtered/transformed batch | D1；read count/base quality oracle；单 batch kernel setup 不得超过计算 20%；风险是改变过滤顺序 |
| `collections/iterators/functional/param` | arena、stable vector、offset iterator、schema descriptors | 无 | Host views/descriptors | D0；ASan/UBSan/property tests；对象图不得进 device |
| `report/tsv/text/help/logging/config/python/R` | 全部保留 CPU/fallback；统一日志和版本 | 无 | text/config/report ↔ files | D2；golden text/header/locale；CPU 时间占比 profile 后再考虑；风险是非核心移植消耗吞吐 |
| `spark` | 由 Nextflow+SLURM 编排；保持原 Spark/Picard fallback | 无 | Spark command → fallback output | D2；process exit/output oracle；不重建 Spark API |
| `bwa/dragstr/illumina/fragments` | seed/read special semantics、模型 table、STR context | 可批量的 encode/score/count primitive | read/features → scores/tables | D2/FG-310；只有 profile 后注册；风险是专项模型的 hidden semantics |

## 5. 数值和核心算法模块明细

### 5.1 PairHMM、SW、alignment/BAQ/clipping

| 模块 | Host 任务 | Kokkos 任务 | 依赖 | 完成定义/性能门禁/风险 |
|---|---|---|---|---|
| `utils/pairhmm` (`FG-100`) | regular/flow dispatch、length bucket、quality table、pair permutation、strict log10、workspace pool | CPU `RangePolicy` + `Kokkos::Experimental::simd<double>` 跨独立 pair；GPU `TeamPolicy` tiled wavefront；三行 M/I/D 或 scratch | FG-020/040 | GATK EXACT/ORIGINAL/LOGLESS/GKL matrix raw-bit 或逐元素 oracle；empty/long/ragged/quality boundary；amortized ≥ 当前 0.94× 16-core baseline 后再优化；风险是 recurrence/reduction 改变一位 bit |
| `utils/smithwaterman` (`FG-110`) | scoring/overhang/tie-break、traceback、CIGAR compaction、变长输出 | affine score `H/E/F`，SIMD across alignments，GPU tile/anti-diagonal；输出 checkpoint 不动态分配 | FG-020/040 | score+CIGAR/start/end 与 GATK/GKL；≥ 90% 当前 355M→4.94G cell/s 强扩展；风险是 checkpoint 不足导致 traceback 错 |
| `utils/alignment` | CIGAR semantics、reference projection、clip/flag、read coordinate | CIGAR op mask、base mapping、短 read transform | FG-020/110 | CIGAR/reference projection oracle；只在 batch ≥ threshold 设备化；风险是 indel 坐标错位 |
| `utils/baq/clipping` | BAQ HMM 边界、clip/filter 顺序 | quality/mask primitive，fixed window reduction | FG-130/110 | GATK BAQ/clipping test corpus；quality histogram exact；小窗口 Host fast path；风险是 BAQ 浮点非确定 |

### 5.2 activity profile、locus iterator、downsampling、pileup

| 模块 | Host 任务 | Kokkos 任务 | 依赖 | 完成定义/性能门禁/风险 |
|---|---|---|---|---|
| `activityprofile` (`FG-130`) | active/inactive state、probability threshold、window boundary | per-locus base/mismatch/quality/mapQ histogram、activity score | FG-020/040 | region boundary/active mask 与 GATK；高覆盖吞吐 ≥ Java 2× 或 baseline；风险是跨 window 状态重排 |
| `locusiterator` | coordinate merge、halo、read order | fixed locus compaction/count | FG-040 | no duplicate/missing loci；I/O+compute pipeline queue depth 可观测；风险是乱序导致 Strict 失败 |
| `downsampling` | reservoir/seed、per-sample/read-group semantics | candidate mask/score only | FG-130 | same seed exact read subset；随机数不可由 backend 改变；风险是 per-team reduction 改变样本选择 |
| `pileup` | filter predicate、CIGAR-aware projection、annotation order | base/quality/mapQ histogram、coverage, mismatch | FG-020/130 | pileup counts/DP/AD exact；小 region 启动成本不超过 Host；风险是 soft clip/indel projection |

### 5.3 BQSR、genotyper、haplotype、dense numerical

| 模块 | Host 任务 | Kokkos 任务 | 依赖 | 完成定义/性能门禁/风险 |
|---|---|---|---|---|
| `recalibration`/BaseRecalibrator (`FG-220`) | read-group/quality/cycle/context descriptor，report/table schema，固定顺序 merge | per-team local integer table，deterministic merge；不默认 global atomic | FG-020/130 | report/table 和 integer counts 与 GATK；比当前 Java wall 至少 2×，内存按 batch 上限；风险是 covariate context/indel event 语义 |
| `ApplyBQSR` | recal table lookup、quality cap、read/output order | quality lookup/transform on large batches | FG-220 | qualities/header/read order exact；短 batch Host；风险是 round/collapse 不一致 |
| `genotyper`/`haplotype` (`FG-140`) | allele/sample order、prior、annotation、special allele、NaN/underflow | log-sum-exp、likelihood transform、genotype math、reference confidence | FG-100/130 | PL/GQ/AD/DP/GT/gVCF block exact or declared tolerance；fixed reduction tree；风险是 vendor BLAS/fast math |
| `pileup` utility | allele evidence schema、sample mapping | counts/likelihood primitive | FG-130/140 | same as above；不重复实现 reader |
| `mcmc` | chain state、seed、accept/reject、checkpoint | batched log probability/proposal score | FG-250 | same seed/restart chain oracle；吞吐只在长链 profile 后门禁；风险是 random order |
| `svd/clustering` | metadata、model state、segmentation/cluster control | block GEMM/SVD/GMM dense math via Kokkos Kernels | FG-250 | tolerance、singular/NaN、segment boundaries；设备内存 block budget；风险是 library reduction order |

## 6. 工具包逐项实施路线

下表覆盖 inventory 中所有主要 `tools` 包。工具均必须先完成对应 kernel 和 output
contract，再在 registry 变为 native；没有列出的命令默认走 GATK fallback，并记录
`FALLBACK_EXECUTED`。

| 工具包/任务 | Host 边界 | Kokkos/adapter | 输入/输出契约 | 依赖 | oracle/性能门禁/风险 |
|---|---|---|---|---|---|
| `walkers/haplotypecaller` (`FG-200`) | AssemblyRegion、filters、downsample、graph/path、SW traceback、VCF/GVCF writer | activity/pileup、k-mer primitives、SW score、PairHMM、genotype | BAM/CRAM+FASTA+interval → VCF/VCF.GZ/GVCF+index+manifest | 100/110/120/130/140/020/030 | GATK HC integration/GIAB/HC test resources；端到端包含 decode/write，不能只算 kernel；风险最高：indel、CIGAR、gVCF |
| `walkers/mutect` (`FG-210`) | tumor/normal mapping、contamination、orientation、filters、stats/F1R2/somatic GVCF | 复用 HC assembly/SW/PairHMM/likelihood；per-team histogram | tumor/normal BAM/CRAM → VCF+stats+F1R2+GVCF+manifest | 200/140 | Mutect2/FilterMutectCalls integration、truth、sidecars 全部存在；风险是 somatic filter hidden state |
| `walkers/bqsr` (`FG-220`) | report/table schema、merge、Apply stream | covariate/table reductions、quality transform | BAM/CRAM → recal report/table、BQSR BAM/CRAM | 130/020 | GATK report/quality/read order；高覆盖 throughput、memory cap；风险见 §5.3 |
| `walkers/root` (`FG-260`) | SelectVariants/CombineGVCFs/Reblock/annotations、VCF writer、allele/block order | locus block likelihood/genotype/field reductions | VCF/GVCF/GenomicsDB → VCF/GVCF | 140/230 | GATK root walkers golden、header/index/block boundaries；风险是 cross-sample order |
| `walkers/vqsr` (`FG-250`) | annotation schema、training resources、seed/convergence/tranche/report | VBEM/GMM expectation/maximization、batch score | VCF+resources → recal/table/tranches/VCF | 250/030 | GATK VQSR report/tranche/VCF；model tolerance + deterministic seed；风险是 BLAS/reduction差异 |
| `walkers/annotator` | plugin lifecycle、data source interval lookup、strings/header | fixed numeric depth/AC/quality summary only | VCF/BAM/resources → annotated VCF | 300/030 | plugin output/header/field order；必须 profile 证明 cache/lookup 是瓶颈；风险是插件生态不可静态化 |
| `walkers/variantutils` | normalize/left-align、VCF conversion/table、header/missing values | optional fixed numeric field count | VCF/BCF → VCF/TSV/normalized VCF | 020/030 | byte/field exact；CPU/cache-first；风险是 representation/left-alignment |
| `walkers/varianteval` | evaluator/stratification/plugin/report | optional aggregate counts | truth+eval VCF → report/metrics | 030/300 | GATK metrics/report；不能为报告逻辑牺牲 native correctness |
| `walkers/sv` | evidence graph/cluster/stratify/output | evidence count、key sort、cluster primitive | BAM/VCF/SV evidence → SV VCF/metrics | 310/030 | SV test corpus/VCF sidecars；先 fallback，单 primitive 独立 oracle；风险是复杂结构事件 |
| `walkers/coverage` | traversal/filter/metrics writer | depth/coverage histogram | BAM/CRAM → coverage metrics/TSV | 130/020 | DepthOfCoverage metrics/header；收益低于 I/O 时保持 Host |
| `walkers/fasta` | FASTA index/header/reference checks | checksum/parallel fixed-base transform | FASTA/dict → FASTA/dict/checksum | 020 | byte/checksum exact；I/O-bound 不强行设备化 |
| `walkers/filters` | predicate semantics、error reporting | bitmask/filter count | reads/variants → filtered records/counts | 130 | read/variant acceptance exact；风险是 predicate order |
| `walkers/validation` | ValidateVariants/Reads header/index/schema | checksum/count primitive | files → validation report/exit code | 020/030 | malformed corpus 和 exit code；不把验证失败吞掉 |
| `tools/spark` | Spark/Hadoop API、partition/reducer、PathSeq orchestration | 无默认 kernel；可调用 FG-310 primitive | Spark command → 原格式输出 | 001/410 | 原 GATK Spark fallback，Nextflow/SLURM process 通过；不重建 DAG |
| `tools/copynumber` (`FG-250`) | HDF5/TSV/sample/interval metadata、MCMC/segmentation/model | CollectReadCounts histogram、Denoise GEMM/SVD、批量 likelihood | BAM/HDF5/TSV → HDF5/segments/plots/metrics | 250/020 | GATK CNV test resources、HDF5 readable、segment boundaries；风险是模型状态/数值 tolerance |
| `tools/funcotator` | data source cache/index、gene/transcript lookup、string rendering | optional numeric consequence summary | VCF/segments+data sources → VCF/MAF/SEG | 300/030 | GATK Funcotator fields/header/MAF；CPU/cache-first；风险是数据源版本漂移 |
| `tools/genomicsdb` (`FG-230`) | GenomicsDB/TileDB API、workspace、fragment/consolidation、fd/tmp | locus block likelihood/genotype only | VCF/GVCF samples → workspace/manifest | 230/030 | original GATK can reopen workspace；fd/memory/tmp stress；不伪造 sparse kernel |
| `tools/reference` | reference dictionary/compatibility、FASTA I/O、checksum | parallel checksum/compare fixed blocks | FASTA/dict → compatibility report/exit | 020 | CheckReferenceCompatibility/CompareReferences exact；风险是 contig alias |
| `tools/dragstr` | model calibration/compose table、seed/table schema | optional batch score/count | reads/STR tables → DragSTR model/table | 310 | GATK DragSTR table/model oracle；profile 后决定 native，否则 fallback |
| `tools/report/unsupported` | 全部 CPU/fallback，统一 manifest 和日志 | 无 | report/text/scripts → original files | 001/030 | exit/output transparent；风险是被误注册为 native |

## 7. 三条主路径的实施顺序

### 7.1 HaplotypeCaller

```text
FG-020 HTSlib + reference
  -> FG-040 traversal/core+halo/filter/downsample
  -> FG-130 activity/pileup/BAQ
  -> FG-120 k-mer encode/sort + Host graph/path
  -> FG-110 SW score + Host traceback
  -> FG-100 PairHMM regular/flow
  -> FG-140 genotype/reference confidence/GVCF state
  -> FG-030 VCF/GVCF/index/manifest
  -> FG-200 registry + GATK oracle
```

每个 region 只能由一个 `PipelinePlan` 持有这些 batch/workspace，不能在 HC Host
中重复 deep-copy。halo 计算但不输出；core 结果按 `(contig, position, stable id)`
排序。初版只登记支持的参数子集，`-ERC GVCF` 在没有 reference-confidence oracle
前继续 fallback。

### 7.2 Mutect2

```text
HC traversal/assembly/SW/PairHMM
  -> tumor-normal likelihood
  -> orientation/contamination/somatic filters (Host)
  -> VCF + .stats + F1R2 + somatic GVCF + manifest
```

sidecar 必须和主 VCF 一起校验；只有 VCF 而没有 stats/F1R2 的结果不能返回成功。

### 7.3 BQSR

```text
HTSlib read stream
  -> Host covariate descriptor
  -> Kokkos local integer tables
  -> deterministic merge/report
  -> Host ApplyBQSR stream
```

table merge 在 shard 结束处可序列化；中断后从最后一个 input ordinal 继续，不重复计数。

## 8. 测试、oracle 和性能门禁

### 8.1 测试层级

| 层 | 任务 | 必须执行的证据 |
|---|---|---|
| Unit/property | schema offset、CIGAR、base encoding、quality table、error enum、stable sort | ASan/UBSan、空/最大/恶意输入、property/fuzz |
| Kernel strict | PairHMM/SW/k-mer/pileup/BQSR/genotype 小向量和边界 | GATK/GKL/reference raw-bit、逐 CIGAR、固定 signature |
| Format | BAM/CRAM/VCF/BCF/BGZF/index/header/PG/sidecar | HTSlib/samtools/GATK read-back，gzip/CRAM checksum |
| Tool integration | HC/Mutect2/BQSR/root/VQSR/CNV corpus | GATK 4.6.2.0 golden、真实 truth、失败和 fallback |
| Resource | cgroup、SLURM CPU/memory/GPU、scratch/inode/fd、spill/retry | 不 OOM 越界、不产生 partial success、manifest 可解释 |
| Pipeline | Nextflow scatter/gather/resume/cache、SLURM allocation/retry、容器 | 原命令格式、process exit、output ordinal、日志/telemetry |
| Cross backend | scalar/AVX2/AVX512/SVE/CUDA/HIP/SYCL | strict signature 或明确 Fast tolerance；所有 backend 能读同一输出 |

### 8.2 Golden corpus 分层

1. `unit/`：合成 reads、空序列、单 base、质量边界、indel/homopolymer、NaN/underflow。
2. `gatk-small/`：`gatk-source/src/test/resources` 中 HC/BQSR/Mutect2/GATK 工具 fixture。
3. `real/`：chr17 69k--70k BAM、CRAM+FASTA、真实 interval/shard。
4. `truth/`：GIAB/confident regions、tumor-normal、cohort gVCF、CNV/VQSR resources。
5. `stress/`：高覆盖、长 read、数千 sample、有限内存/scratch/fd、断点恢复。

golden 记录输入 SHA256、reference/dictionary SHA256、GATK 版本、命令 argv、backend、
determinism mode、输出 hashes 和可接受差异。oracle 脚本不应只比较一个 sentinel。

### 8.3 性能门禁

- 单 kernel：相对当前 C++/Kokkos baseline 不得回归超过 10%；报告 cells/observations/
  pairs/s、bytes/s、occupancy（设备时）、prepare/execute。
- 端到端：相对 GATK Java/GKL 报告 wall、CPU-hours、RSS、scratch、I/O bytes、网络
  bytes；P0 先达到可替换正确性，P1 才要求工具级 speedup，不能用一项 kernel 速度
  代表整个 tool。
- 扩展性：1/2/4/8/16/32 threads；强扩展效率、弱扩展、batch size sweep、region
  size sweep；超过 2 倍核数仍下降时必须记录原因而不是隐藏点。
- 设备：host→device copy 和初始化单列；若 prepare 占 wall > 30%，必须先做 plan/cache
  优化；device OOM 必须回退或 fail-closed，不能在同一预算无限重试。

## 9. 兼容、集群和发布执行

### 9.1 Registry promotion gates

工具从 `fallback` 到 `native-candidate` 需要：CLI 参数白名单、输出/sidecar schema、
Strict oracle、格式 read-back、资源压力、至少一轮多线程签名、benchmark 报告。
从 `native-candidate` 到默认 `native` 还需要真实 Nextflow/SLURM E2E 和 GATK truth。
任何阶段失败则保留 fallback 并在 manifest 写明原因。

### 9.2 Nextflow/SLURM

- Nextflow process 仍执行 `gatk HaplotypeCaller ...`，launcher 读取 `SLURM_CPUS_PER_TASK`、
  `SLURM_MEM_PER_NODE`、`SLURM_TMPDIR` 和 GPU GRES，不在作业内自行 `sbatch`。
- scatter 以 sample/interval/shard 为边界；gather 只接收 manifest complete 的产物。
- `resume` 使用 input/reference/argv/backend/schema 的 hash；retry 不能复用不完整 output。
- 多节点只做作业级并行，Kokkos 不做细粒度跨节点同步；远端 I/O 计入 queue/backpressure。

### 9.3 产物矩阵

| 产物 | Kokkos backend | 优先模块 |
|---|---|---|
| `fast-gatk-cpu-x86` | OpenMP + scalar/AVX2/AVX512 | 全部 P0/P1 |
| `fast-gatk-cpu-arm` | OpenMP + scalar/SVE-safe | P0 correctness，再做 SIMD |
| `fast-gatk-cuda` | CUDA + host OpenMP | PairHMM/SW/k-mer/genotype |
| `fast-gatk-hip` | HIP + host OpenMP | 同上，需独立 oracle |
| `fast-gatk-sycl` | SYCL + host backend | correctness first |

## 10. 每个模块的 PR/交付清单

提交任何模块前必须附：

1. schema 和 Host/Device ownership 说明；
2. `prepare/execute/collect` API 和 plan 生命周期测试；
3. strict oracle 命令、golden 输入 hash、差异报告；
4. output/sidecar/manifest read-back 测试；
5. 资源预算、队列背压、spill/retry 行为；
6. 1/2/4/8/16/32 thread 和 batch-size benchmark，分列 kernel/prepare/E2E；
7. fallback 条件、错误码、兼容参数和 registry 更新；
8. 未实现范围和后续任务，不得在 README 中使用“兼容 GATK”而没有门禁证据。

## 11. 第一轮实际执行清单（可直接开工）

按以下顺序创建 issue/PR，每项完成后运行对应脚本并把 JSON 放到 results：

1. `FG-001/010/030`：registry、ResourceBudget、OutputManifest/atomic writer。
2. `FG-020`：扩展 HTSlib `ReadBatch` 为 CIGAR-aware/index seek + VCF/BCF writer。
3. `FG-040`：实现 core/halo、byte-bounded queues、deterministic shard replay。
4. `FG-100`：从 demo 提炼 PairHMM regular/flow Plan、persistent workspace、GATK table。
5. `FG-110`：SW checkpoint/CIGAR Host collector，随机和真实 assembly oracle。
6. `FG-120`：k-mer sort/unique/edge compact + Host graph/path，并接入 assembly region。
7. `FG-130/140`：activity/pileup/BAQ/downsample 和 genotype/reference confidence。
8. `FG-200`：用真实 HC 资源替换当前 bounded local-caller prototype，先支持有限参数子集。
9. `FG-220/210`：BQSR、Mutect2 和 sidecar；分别注册，不共用未验证的 output contract。
10. `FG-230/250/240`：GenomicsDB/GenotypeGVCFs、CNV/VQSR、external-memory Sort/MD。
11. `FG-300/310/320`：annotation/Funcotator/variantutils/varianteval/SV/BWA/PathSeq 长尾。
12. `FG-400/410`：GPU correctness/optimization、真实 SLURM/Nextflow E2E、registry promotion。

P0 的停止条件不是“demo 能跑”，而是 `FG-200` 在指定 GATK golden、VCF/GVCF/index、
资源压力和真实 process 中通过；否则 HC 继续 fallback，后续模块仍可并行开发 kernel，
但不能宣称整个 fast-gatk 已经可直接替换 GATK。

## 12. 风险台账和应对

| 风险 | 触发信号 | 应对/决策 |
|---|---|---|
| FP 运算顺序改变 | raw-bit mismatch、GQ/PL 漂移 | Strict 固定函数/归约树；Fast 单独标识；不能用 epsilon 掩盖 schema 错 |
| CIGAR/indel 语义错误 | SW score 对但 CIGAR/VCF 错 | CIGAR-aware Host oracle 先于 GPU；完整 haplotype corpus |
| prepare/IO 吞掉 kernel 收益 | kernel 快但 E2E 慢，prepare >30% | persistent Plan、packing/cache、队列和批量阈值；按 bytes 调节 |
| device memory/cgroup 超限 | OOM、swap、scratch 满 | ResourceBudget、chunk/spill、device→CPU fallback 或 fail-closed |
| global atomic/graph contention | 高核数 throughput 下降 | local table、ownership、sort/compact；保留 Host traversal |
| registry 误选 native | 不支持参数仍生成结果 | 白名单 + capability check + fail-closed + fallback reason |
| sidecar/索引缺失 | VCF 存在但 downstream 失败 | OutputManifest 原子状态；gather 只接收 complete |
| 外部依赖漂移 | HTSlib/GATK/GenomicsDB 版本改变结果 | pin 版本和 hash；manifest 记录 provenance；CI 定期重跑 golden |
| 长尾移植消耗过大 | profile 无明显 kernel 热点 | 保持 fallback；只按独立 primitive/收益进入 registry |
| 真实集群差异 | 本地 smoke 通过、SLURM/Nextflow 失败 | 真实 allocation、cgroup、retry、remote I/O E2E gate |

## 13. 完成判定

全项目只有在以下条件全部满足时才可称为“按计划完成”：

- inventory 中每个模块都有 registry 状态：`native`、`adapter` 或带原因的 `fallback`；
- P0/P1 工具的 Strict oracle、格式/sidecar、资源、pipeline gate 全部通过；
- 不同 Kokkos backend 使用同一源代码，结果差异符合工具定义；
- Nextflow/SLURM 原命令可执行，scatter/gather/resume/retry 不丢失结果；
- 所有 benchmark 同时提供 kernel、prepare、I/O、E2E、CPU-hours 和内存/scratch/网络；
- 尚未实现的 biological semantics（indel、gVCF、somatic filters、模型训练等）不能
  用 smoke 结果替代，必须继续 fallback 或单独标注 prototype。
