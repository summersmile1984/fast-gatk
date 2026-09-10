# fast-gatk 运行时和 Kokkos API 草案

按模块落地的 C++ Host + Kokkos kernel 实施蓝图见
[MODULE_IMPLEMENTATION_PLANS.md](</home/turing-agents/Documents/fast-gatk/MODULE_IMPLEMENTATION_PLANS.md>)。

这份草案把调研结论收敛为可以开始编码的边界。它不绑定具体 Kokkos backend；同一套 tool contract 在 OpenMP、CUDA、HIP、SYCL 上复用。

## 组件分层

```text
compat/       gatk launcher, argument files, Java fallback, output contract
tools/        HaplotypeCaller, Mutect2, BQSR, GenotypeGVCFs adapters
traversal/    interval, shard, core/halo, deterministic output ordering
io/           HTSlib, BGZF, reference/index, cloud range reads, spill runs
kernels/      PairHMM, Smith-Waterman, k-mer, reductions, dense math
runtime/      resource probe, adaptive controller, queues, telemetry
backends/     Kokkos execution/memory spaces and optional vendor fast paths
```

依赖关系必须单向：`kernels` 不依赖 `io` 或 SLURM；`io` 不调用 Kokkos worker；`compat` 可以选择 native 或 fallback，但 native tool 不应再次调用 launcher。

## 资源探测接口

```cpp
namespace fastgatk::runtime {

struct ResourceBudget {
  uint64_t host_hard_bytes;       // cgroup memory.max or SLURM memory
  uint64_t host_target_bytes;     // hard minus reserve
  uint64_t device_hard_bytes;     // backend-reported free memory cap
  uint64_t scratch_hard_bytes;    // configured fraction of free local scratch
  uint64_t scratch_inode_limit;
  uint32_t cpu_threads;           // SLURM_CPUS_PER_TASK, validated against affinity
  bool local_ssd;
  bool remote_input;
};

class ResourceProbe {
 public:
  ResourceBudget snapshot() const;
  uint64_t resident_bytes() const;
  uint64_t scratch_free_bytes() const;
  uint64_t scratch_free_inodes() const;
  double io_read_latency_ewma_ms() const;
  double io_write_latency_ewma_ms() const;
};

}  // namespace fastgatk::runtime
```

探测优先级：cgroup v2 → SLURM 环境 → host total memory。不能把物理机总内存当作 job 可用内存。`TMPDIR` → `SLURM_TMPDIR` → 配置文件；不存在可写 local scratch 时，spill 应显式降级到远端或失败，不得静默写到共享根目录。

## 自适应控制器

```cpp
struct BatchLimits {
  uint32_t max_reads;
  uint32_t max_haplotypes;
  uint64_t max_host_bytes;
  uint64_t max_device_bytes;
  uint64_t max_inflight_bytes;
};

enum class Pressure { Normal, HostMemory, DeviceMemory, Scratch, FileDescriptors, IoBound };

class AdaptiveController {
 public:
  BatchLimits initial(const ResourceBudget&, const WorkEstimate&) const;
  BatchLimits next(const Telemetry&, const BatchLimits&, Pressure) const;
  bool should_spill(const Telemetry&) const;
  bool should_pause_reader(const Telemetry&) const;
};
```

控制规则：

- 初始目标为 host target 的 60%、device free 的 60%、scratch free 的 80% 上限。
- 连续两个采样窗口超过 80% 使用率时减半；低于 50% 且计算队列空时每次增加 10–20%。
- scratch 或 inode 达阈值时停止扩大 batch，优先完成现有 run 并合并。
- 分配失败只允许缩小 batch 后重试一次；再次失败返回 `RESOURCE_EXHAUSTED`，让 Nextflow/SLURM retry 策略处理。
- controller 只调整安全边界：assembly region、read-name group、VCF block、GVCF locus block，不能在对象中间截断。

## 三队列 I/O 设计

```cpp
template<class Record>
class BoundedQueue;  // byte capacity, not only item count

BoundedQueue<ReadChunk> decoded_reads;
BoundedQueue<KernelBatch> compute_batches;
BoundedQueue<ResultChunk> encoded_results;

// reader: HTSlib/BGZF/reference range reads
// compute: Kokkos execution space
// writer: deterministic reorder + BGZF/libdeflate + index writer
```

队列以字节容量为主，避免一条大 read/chunk 造成 item-count 失真。每项携带：`job_id`、`shard_id`、`core_interval`、`halo_interval`、`input_offset`、`output_ordinal`、`schema_version`、CRC。

Kokkos execution space 只消费已解码的 `KernelBatch`；reader/writer 运行在独立 CPU pool。设备 staging 使用 pinned host memory，但 pinned bytes 必须计入 `host_target_bytes`。

## PairHMM API

```cpp
struct PairHMMBatch {
  Kokkos::View<const uint8_t*> read_bases;
  Kokkos::View<const uint8_t*> read_qualities;
  Kokkos::View<const uint32_t*> read_offsets;
  Kokkos::View<const uint8_t*> hap_bases;
  Kokkos::View<const uint32_t*> hap_offsets;
  Kokkos::View<const uint32_t*> pair_read_ids;
  Kokkos::View<const uint32_t*> pair_hap_ids;
  ErrorModel model;                 // regular or flow, never implicit
};

template<class ExecSpace>
PairHMMResult compute_pair_hmm(const PairHMMBatch&, ExecSpace,
                               DeterminismMode);
```

实施规则：

- 先按 read/haplotype 长度和 error model 分桶。
- CPU 基线用 `RangePolicy` 每个 work item 处理一个 SIMD pair group；GPU 性能版再按 tile/wavefront 使用 team policy。DP workspace 必须是 Kokkos 管理的 memory-space allocation。
- 连续 haplotypes 的 common prefix 用显式 `hap_start_index`/prefix cache 表示。
- `DeterminismMode::Strict` 固定归约、排序和输出矩阵顺序；`Fast` 允许 backend 浮点顺序差异，但必须通过 tolerance oracle。
- regular PairHMM 和 flow PairHMM 使用不同 kernel tag、测试向量和资源估算器。

## Smith-Waterman API

```cpp
struct AlignmentRequest {
  Span<const uint8_t> reference;
  Span<const uint8_t> alternate;
  SWParameters parameters;
  SWOverhangStrategy overhang;
  uint64_t stable_id;
};

struct AlignmentResult {
  int32_t score;
  int32_t ref_start, ref_end;
  int32_t alt_start, alt_end;
  Cigar cigar;
};

template<class ExecSpace>
void align_batch(Span<const AlignmentRequest>, Span<AlignmentResult>,
                 ExecSpace, DeterminismMode);
```

首版必须以 GATK Java/GKL 输出的 score、CIGAR、overhang 和 tie-break 为 oracle；不要只验证 score。

## Tool contract

```cpp
struct ToolContext {
  ResourceProbe resources;
  AdaptiveController controller;
  IoContext io;
  TraversalPlan traversal;
  DeterminismMode determinism;
  Backend backend;
};

struct ToolResult {
  ExitCode code;
  OutputManifest outputs;       // primary, index, stats, metrics, sidecars
  TelemetrySummary telemetry;
};

ToolResult run_haplotype_caller(const HaplotypeCallerArgs&, ToolContext&);
ToolResult run_mutect2(const Mutect2Args&, ToolContext&);
ToolResult run_base_recalibrator(const BqsrArgs&, ToolContext&);
```

`OutputManifest` 在关闭 writer 前必须验证：文件存在、索引类型正确、header/contig dictionary 一致、core interval 不重复、排序顺序合法、sidecar 数量符合工具契约。验证失败时不应把 partial output 标记为成功。

当前读文件 QC 工具的 contract gate 已落地：`CountReads` 和 `FlagStat` 使用同一套
HTSlib bounded decode 与 C++ Host/Kokkos compute 边界；前者以 `RangePolicy` 做记录归约，
后者以 12-counter atomic kernel 计算 FlagStat 计数。两者的 registry、binary
OutputManifest/telemetry 和 `--help` 均标记 `contract-compatible`，并以 pinned chr17
fixture 对重复 `-I`、`-L/-XL`、interval padding/merging、已实现 read filters、空/坏输入
边界和 GATK 4.6.2.0 oracle 做验证。该等级只覆盖已声明的 CLI、文本输出、遍历、过滤、
确定性和 manifest 契约；remote/cloud reader 及未声明的 GATK filter/launcher 行为必须
继续走显式 fallback。

## CMake/Kokkos 构建矩阵

```text
fast-gatk-cpu-x86  Kokkos_ENABLE_OPENMP=ON, native SIMD baseline
fast-gatk-cpu-arm  Kokkos_ENABLE_OPENMP=ON, no x86-only assumptions
fast-gatk-cuda     Kokkos_ENABLE_CUDA=ON
fast-gatk-hip      Kokkos_ENABLE_HIP=ON
fast-gatk-sycl     Kokkos_ENABLE_SYCL=ON
```

The native build entry point selects these artifacts with
`FASTGATK_KOKKOS_BACKEND=cpu|serial|cuda|hip|sycl` (or an explicit
`*_openmp` variant) and keeps the same C++/Kokkos source and tool registry.
The launcher explicitly resets non-selected CMake backends when reusing a build
directory, so switching CPU/OpenMP to serial cannot silently retain OpenMP.
Missing device compilers/runtimes fail at build or backend initialization and
are never silently reported as a successful accelerator run.

每个产物发布同一版本的 `tool_registry.json`、schema 和 golden corpus。launcher 根据 `FAST_GATK_BACKEND`、设备可见性和输入规模选择产物；找不到可用 native backend 时回退 Java，不把 CUDA 初始化失败当作输入错误。

## 测试层次

```text
kernel oracle       PairHMM/SW exact/reference/GKL vectors
format oracle       samtools/htslib validate, index/read-back
tool oracle         GATK integration tests and expected VCF/GVCF/stats
pipeline oracle     Nextflow + SLURM scatter/gather/resume/retry
resource tests      cgroup cap, tiny scratch, device OOM, interrupted spill
```

每个工具在 registry 中登记支持参数、输出文件、索引、sidecar、determinism、backend 和测试集；未登记工具只能走 fallback。这样“兼容”是可验证的状态，而不是 launcher 通过后就算完成。

## PairHMM 首个可运行 kernel

`pairhmm-demo/` 是这份 runtime spec 的最小实现样板，当前主路径已经统一为 Kokkos 5.2.0 API：

```text
PairHMMInput batch
      -> host validation / quality-table selection
      -> Kokkos Views/deep_copy + RangePolicy
      -> Kokkos SIMD scalar/AVX2/AVX-512 independent-pair kernel
      -> likelihoods + checksum + telemetry
```

对应的后端合同：

```cpp
enum class PairHMMMode { Strict, Fast };

struct PairHMMContext {
  int workers;                 // must equal the SLURM CPU allocation used by the process
  PairHMMMode mode;
  DeterminismMode determinism;
  const uint64_t* gatk_table_hash;
};

PairHMMResult compute_pair_hmm(const PairBatch&, const PairHMMContext&);
```

`Strict` 使用与 GATK 运行时导出的 IEEE-754 quality/transition table、OpenJDK fdlibm log10 和固定运算顺序，要求 raw bit comparison。SIMD 以独立 pair 为 lane，不改变单 pair DP 的依赖和最终累加顺序；当前 4096-value matrix oracle 在 Kokkos scalar/AVX2/AVX-512 上均为 `bit_different=0`。Java oracle 固定 `-XX:DisableIntrinsic=_dlog10`，因为默认 HotSpot x86 Math intrinsic 本身会产生 3/4096 个 1-ULP 平台差异。实测与 GATK/GKL 对照见 [PAIRHMM_SIMD_BENCHMARK.md](</home/turing-agents/Documents/fast-gatk/PAIRHMM_SIMD_BENCHMARK.md>)。

### Host 语言边界决策

Kokkos core 和 host/runtime 首选 C++：避免在 kernel、每条 read 或每个 event 上跨 Rust/C++ ABI，并直接复用 Kokkos 的 execution/memory space。若未来需要 Rust，限制为批次级 C ABI（region/shard 一次提交连续 buffer）或独立编排进程；不得把 FFI 放入 DP cell/read 热循环。C ABI 仍可作为稳定的外部兼容边界，但不是内部调度层。

### CPU backend 分派

CPU backend 当前从同一 `pairhmm_kokkos.cpp` 构建 generic、ZEN3/AVX2、ZEN4/AVX-512 三个产物；算法不含 raw intrinsics，launcher 后续根据 CPU feature 选择已构建的 Kokkos architecture。测试入口是 `python3 pairhmm-demo/scripts/verify_kokkos.py`。AVX-512 不能只按向量宽度推导收益，但最终同核 benchmark 显示它在本机已优于 AVX2：independent kernel 是 GATK Java 的 6.89–7.81×，matrix8 kernel 是 GKL 实际矩阵吞吐的 1.09–1.28×。

准备阶段仍是明确瓶颈：当前一次 input packing、allocation 和 deep-copy 摊到十次 matrix iteration 后，在 16 核为 GKL 的 0.94×，而 steady-state kernel 为 1.28×。runtime 必须持久化 Kokkos 实例、池化三行 workspace、缓存 read quality transitions，并以 region/batch 级调用摊销准备成本。GPU 后端下 `simd<T>` 是 scalar correctness baseline；优化版需要 Kokkos team/vector wavefront，不能伪称当前 CPU SIMD kernel 已获得 GPU 性能可移植性。

## 四类代表模块验证（Kokkos 5.2.0）

为确认上述模式不只适用于 PairHMM，仓库新增
[kokkos-modules-demo](</home/turing-agents/Documents/fast-gatk/kokkos-modules-demo/README.md>)：

- Smith-Waterman：affine-gap 规则 DP、滚动行 workspace、`RangePolicy`（当前验证最大分数，CIGAR/tie-breaking 尚未接入）；
- BQSR 风格统计：质量/周期 binning、`parallel_reduce` 和整数 `atomic_add`；
- local assembly graph：k-mer edge list、竞争度数更新和不规则图归约；
- region stream：输入适配、region 分块、pileup 原子更新，并拆分 I/O 与 compute 计时。

三个 CPU 变体（generic、ZEN3/AVX2、ZEN4/AVX-512）均使用同一份 kernel
源码，验证脚本逐项对照主机参考实现并比较跨变体结果签名：

```bash
bash kokkos-modules-demo/scripts/build_variants.sh
python3 kokkos-modules-demo/scripts/verify_modules.py --records=512 --threads=4
```

当前验证结果记录在
[verify-modules.json](</home/turing-agents/Documents/fast-gatk/kokkos-modules-demo/results/verify-modules.json>)；
三种变体四个模块全部通过。AVX-512 代表性 benchmark 记录在
[benchmark-avx512.json](</home/turing-agents/Documents/fast-gatk/kokkos-modules-demo/results/benchmark-avx512.json>)。

当前已进一步提炼出共享的 [fastgatk-core](</home/turing-agents/Documents/fast-gatk/fastgatk-core/README.md>) library target；
PairHMM 和四个原型都使用 `HostBatch/DeviceBatch/KernelPlan` 记录 prepare、
device bytes、execute calls 和 kernel 时间。真实 HTSlib/BAM/CRAM 文件边界由
[fastgatk-hc-smoke](</home/turing-agents/Documents/fast-gatk/fastgatk-native/README.md>)
验证；它现在覆盖 GATK `-I/-R/-L/-O` 参数别名、VCF/VCF.GZ 文件边界、
OutputManifest/telemetry、Nextflow 本地 process 等价路径和 SLURM 本地资源模拟，
但仍不是完整 HaplotypeCaller。Nextflow 通过
`fastgatk-native/workflow/nextflow.config` 默认使用 local profile，生产集群
显式使用 `-profile slurm`；这样由 Nextflow 管理 sbatch/squeue allocation，
`slurm_smoke.sh` 只负责参数/资源边界适配，不在已有 allocation 内二次提交。
下一阶段的 `fastgatk-hc-call` 已接通单 SNP assembly → likelihood → genotyping，
支持 VCF/VCF.GZ、GATK 参数别名和 manifest；`verify_gatk_oracle.py` 对同一 GATK
fixture 做结构性 oracle，但当前明确不宣称 bit-identical 或完整 HC 生物学等价。
结果显示规则 Smith-Waterman 随 CPU 核数扩展较好；BQSR 的细粒度原子
binning、图度数更新和很小的 region batch 在高核数会受同步/调度开销限制，
因此生产实现需要分桶、局部归约、批量化和持久化 workspace，不能简单复制
PairHMM 的并行粒度。


### Nextflow/SLURM failure boundary（2026-08-31）

当前 scheduler smoke runner 的可审计实现是：每个 shard/gather attempt 使用 outdir
下的私有目录，失败后清理，成功后以 primary/index → manifest 的 manifest-last
rename 顺序发布；旧 bundle 在提交期间进入同一 attempt 的备份目录，提交失败时回滚。
--retries N 是有界的本地/fake-SLURM 模拟，成功 manifest 的 workflow signature
包含 input/reference/binary/gather binary/threads 和 shard manifest SHA-256；只有
signature 与 complete manifest 同时匹配才允许 resume。
fastgatk-scheduler-retry-contract 验证 retry、gather retry、无执行 resume 和 retry
exhaustion fail-closed。它不宣称真实 SLURM 多节点、preemption 或 cgroup failure
已验证；生产环境仍由 Nextflow 的 errorStrategy retry/maxRetries 管理 process 重试，
slurm_smoke.sh 只做资源边界适配且不会在已分配 job 内二次 sbatch。
