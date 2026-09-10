# FastGATK 硬件与性能执行计划（Draft）

日期：2026-08-21

本文把 GATK_CPU_SIMD_GPU_DETAIL.md、PAIRHMM_SIMD_BENCHMARK.md、Kokkos module
validation 和当前 benchmark 脚本转成可执行路线。目标不是只优化一个 kernel，而是
让数据准备、I/O、内存、网络、kernel、writer、SLURM 资源在同一个性能门禁下可测量、
可回归、可解释。

## 0. 当前基线与必须保留的事实

| 项目 | 当前证据 | 对路线的含义 |
|---|---|---|
| CPU | AMD Ryzen 9 7945HX，1 socket，16 physical / 32 logical，1 NUMA node | 可做物理核 1/2/4/8/16 的可重复基准；不能把本机结果外推到双路 NUMA |
| ISA | avx2、avx512f/dq/bw/vl、avx512_vnni/bf16 等 | AVX-512 可用，但需同时测频率、功耗和 AVX2；向量宽度不等于端到端收益 |
| Kokkos | 5.2.0；OpenMP + Serial 已启用 | 当前只有 CPU 实测；CUDA/HIP/SYCL 仍需设备/编译器 CI |
| CPU 变体 | generic/scalar、Kokkos_ARCH_ZEN3、Kokkos_ARCH_ZEN4 | Kokkos SIMD 通过不同构建获得 lane width 1/4/8；不是一个可任意 ISA runtime dispatch 的 fat binary |
| 热路径 | PairHMM 使用 View、RangePolicy、Kokkos SIMD；module 使用 HostBatch/DeviceBatch/KernelPlan | 继续只使用 Kokkos API；raw SIMD intrinsics 不进入生产路径 |
| strict 正确性 | scalar/AVX2/AVX-512 对 512 records、8×8 matrix、4096 likelihood 已 bit_different=0 | strict 保持 FP64、固定表、固定递推顺序和浮点环境；GPU 先做 tolerance/bio-identical |
| PairHMM steady-state | Kokkos AVX-512 相对 GKL 实际 matrix work 在 1/2/4/8/16 核为 1.09–1.28× | kernel 有潜力，但不能把它当作完整 HaplotypeCaller 加速 |
| PairHMM amortized | packing/allocation/deep-copy 摊入十次计算后为 GKL 的 1.05/1.02/1.01/0.99/0.94× | 第一优先级是持久 workspace、转换缓存和 batch 复用 |
| module scaling | SW 规则 DP 接近线性；BQSR 原子 binning、图度数和小 region 在高核数受同步/调度限制 | 不同模块需要不同粒度和归约策略 |
| I/O | HTSlib BAM/CRAM reader 已接通；VCF writer/manifest 已有 smoke/prototype 边界 | 必须分离 decode/prepare/kernel/encode/index，不能只报告 kernel 秒数 |

当前基准脚本已经有 warmup、固定 affinity、交替 baseline/native 和 median，但还需
统一硬件快照、I/O/内存/功耗计数器、p95 以及自动 gate。当前数字应标为 CPU
kernel/API slice，不代表完整 HaplotypeCaller、Nextflow 或远程对象存储端到端性能。

## 1. 总体执行模型

~~~text
Nextflow process
  -> SLURM allocation/cgroup (CPU, memory, GPU, local scratch)
  -> GATK-compatible launcher / backend registry
  -> resource probe + bounded controller
  -> HTSlib/reference reader (CPU)
  -> flat HostBatch / persistent arena
  -> Kokkos DeviceBatch + KernelPlan
  -> CPU SIMD or Kokkos team/vector GPU kernel
  -> deterministic host reorder / VCF-BAM writer / index
  -> OutputManifest + telemetry
~~~

职责边界：

- Nextflow 负责样本、contig、interval 的 DAG、cache、resume 和 retry；不在 kernel
  内启动 sbatch。
- SLURM 负责 allocation、cpuset、memory cgroup、GPU、local scratch 和失败重试；
  native 进程只读取 SLURM_CPUS_PER_TASK、SLURM_MEM_PER_NODE、SLURM_TMPDIR、
  CUDA_VISIBLE_DEVICES 等信息。
- reader/writer/VCF header/index、CIGAR、字符串、排序和状态机留在 Host；只有固定
  布局数组和数值结果进入 device。
- kernel 不持有 Java/Rust/C ABI 对象，不访问 HTSlib record 指针，不执行 device-side
  动态 new/vector，不依赖另一个 team 的 spin-lock。
- KernelPlan 必须报告 prepare bytes/time、device bytes、execute time/calls；benchmark
  不能把 deep_copy 隐藏在 kernel throughput 里。

### 两种数值契约

| 模式 | 用途 | 规则 | 验收 |
|---|---|---|---|
| strict | GATK pure Java oracle、回归、稳定输出 | FP64；固定 quality/transition table；固定 M/I/D 递推和归约；禁 fast-math/FMA；固定 MXCSR；必要时 scalar fallback | PairHMM raw IEEE bits；VCF record/header/index 与 golden 一致 |
| fast | GKL/硬件替换、吞吐优先 | 允许 backend 归约顺序、FTZ/DAZ、SIMD/tile 变化，需定义误差和生物学输出合同 | likelihood tolerance、variant/GT/DP/AD/PL 字段 oracle |

Kokkos::Experimental::simd 只承担 CPU lane 抽象。GPU 默认 SIMD ABI 可能是 scalar
lane；GPU 性能版使用 TeamPolicy、TeamThreadRange、ThreadVectorRange、scratch 和
tile/wavefront，不能把 CPU SIMD 源码直接当成 GPU 优化完成。

## 2. CPU 路线：generic、AVX2、AVX-512、ARM

### 2.1 构建和分派

每个 CPU 产物使用同一源代码和 API，编译配置不同：

~~~text
cpu-generic  -> all supported CPUs, strict scalar fallback
cpu-avx2     -> requires AVX2, Kokkos SIMD width 4 on x86 FP64
cpu-avx512   -> requires AVX-512 contract, width 8 where supported
cpu-arm-sve  -> ARM SVE target, width is hardware/compiler dependent
~~~

对应 configure 变体是：generic（OpenMP+Serial）、ZEN3（AVX2 candidate）、ZEN4
（AVX-512 candidate），ARM runner 再加 ARMV84_SVE/实际 compiler target。launcher 根据
FAST_GATK_BACKEND、CPU feature、cpuset 和输入规模选择 artifact；缺少 ISA/binary 时
回落 generic/Java，不能执行后收到 SIGILL 才 fallback。artifact 记录 compiler、
Kokkos、KOKKOS_ARCH、strict/fast、Git SHA 和依赖 hash。

### 2.2 PairHMM 优化顺序

1. 持久 arena：三行 M/I/D、quality/transition table、pair index、host staging 从
   pool 获取；一次 prepare 后重复 matrix/batch。
2. 按 read/haplotype 长度、error model、sample/region 分桶，减少 SIMD tail、padding
   和 branch；增加 common-prefix/hapStartIndex cache 及 permutation。
3. 以独立 read-haplotype pair 为 SIMD lane，保留每 lane DP 依赖及 (A+B)+C 顺序；
   不沿单矩阵错误的依赖方向 vectorize。
4. 测 AVX-512 的频率下降、功耗和 bandwidth；若 AVX2 端到端更高，auto 选择 AVX2，
   显式 AVX-512 仍可复测。
5. common-prefix cache、pool 和 batch 调度完成后，才评估 anti-diagonal/tile。

### 2.3 CPU 资源与 NUMA

- 默认 physical cores，不把 SMT sibling 当作额外线性容量；由 SLURM_CPUS_PER_TASK
  设置 Kokkos threads，并校验 affinity。
- OpenMP 记录 OMP_PROC_BIND=close/spread 和 OMP_PLACES=cores；PairHMM 先用 close，
  I/O+compute pipeline 再测 spread。
- 当前只有一个 NUMA node；双路机器须 first-touch read batch、workspace、table、writer
  buffer，并分别报告 local/remote memory。
- batch 根据 bytes 和工作量调节；不能只按 records。大 batch 不能超过 LLC/NUMA budget，
  小 batch 不能让 launch/reader 占据大部分时间。

## 3. GPU 路线：Kokkos CUDA/HIP/SYCL

当前 cache 为 Kokkos_ENABLE_CUDA=OFF、HIP=OFF、SYCL=OFF，因此 GPU 目前只是路线设计。
每种设备单独 configure/build，不能用 CPU OpenMP binary 伪充 GPU：

~~~text
CUDA: Kokkos_ENABLE_CUDA=ON, Kokkos_ENABLE_CUDA_LAMBDA=ON, actual Kokkos_ARCH
HIP:  Kokkos_ENABLE_HIP=ON, actual gfx target
SYCL: Kokkos_ENABLE_SYCL=ON, actual oneAPI/device target
~~~

runner 必须附 nvidia-smi/rocminfo/SYCL device report、driver、toolkit、GPU memory 和
architecture；没有设备时只做 compile smoke，不发布 GPU throughput。

### 3.1 投放顺序

| 阶段 | 投放内容 | Host 保留内容 | 启用条件 |
|---|---|---|---|
| G0 | View/deep_copy/fence、memory telemetry | 全部 I/O | build、device init、checksum 通过 |
| G1 | PairHMM fixed-length buckets，team/tile wavefront | bucket、permutation、compaction | pair 数大、长度窄、copy 可摊薄 |
| G2 | SW score tile + argmax checkpoint | traceback、CIGAR、overhang、tie-break | alignment batch 足够大 |
| G3 | k-mer encode/sort/RLE、BQSR local histogram | graph traversal、deterministic merge、report | flat arrays 足够大 |
| G4 | likelihood/genotype numeric batch | allele order、annotation、VCF/GVCF、state machine | 同 shape matrix、显存和顺序可控 |

GPU 规则：host 先将 BAM/CRAM/CIGAR/字符串变成 SoA；使用 pinned staging arena、异步
copy 和 execution-space instance；device 预分配 workspace；输出固定宽度数值数组；
变长 CIGAR/haplotype/VCF 在 host compaction；OOM、初始化失败、输入过小、长度过宽时
回 CPU并写 fallback_reason。

### 3.2 CPU/GPU 阈值

~~~text
T_gpu = copy_in + kernel_gpu + copy_out
T_cpu = prepare_cpu + kernel_cpu
启用 GPU 当 T_gpu < T_cpu * (1 - safety_margin)
且 pair_count >= min_pairs、padding <= max_padding、
device_free >= estimated_batch_bytes + reserve
~~~

初始 safety_margin=0.20、reserve=15% device memory、max_padding=1.5；由 calibration
数据拟合并写入版本化 registry。小 region、低覆盖、随机 seek 和单条 read 不上 GPU。

## 4. 内存路线

### 4.1 预算模型

~~~text
host_hard = cgroup v2 memory.max or SLURM memory allocation
host_target = 0.75 * host_hard
host_live = decoded reads + CIGAR/reference + batch arenas
           + pinned staging + kernel workspace + writer/index + page-cache allowance
device_target = device_total - driver/reserve - staging reserve
device_live = input Views + workspace + scratch + output
~~~

不能把物理机总内存当作 job 可用内存，也不能排除 pinned bytes、page cache、HTSlib
decoder、writer buffer。缺少 cgroup/SLURM limit 时只标记 local-development。

### 4.2 数据结构和池化

- flat SoA + offsets：bases、qualities、cigar_ops、positions、tids、mapq；不使用
  每条 read 的 heap object。
- PairHMM 按长度 bucket 复用 transitions、three-row workspace、pair index、output。
  DeviceBatch 绑定 typed Kokkos Views，KernelPlan 记录 prepare/device bytes。
- graph/sort 使用 per-thread bins 和 flat edge/key arrays；合并时按稳定 key 排序。
- output reorder queue 以 byte capacity 限制；partial writer 有 schema/version/CRC，
  manifest 在 writer/index/header 验证后才 complete=true。
- 分配失败只允许在安全边界缩小 batch 后重试一次；仍失败返回 RESOURCE_EXHAUSTED。

### 4.3 内存门禁

每个 batch 记录 host live/peak、device live/peak、arena reuse、allocation count、
spill bytes、queue depth、page-cache hint、prepare bytes/time。

- max RSS <= 0.80 × SLURM/cgroup hard limit；
- device peak <= 0.85 × usable device memory；
- 无 unbounded queue、device OOM、partial output；
- workspace reuse >=95%，medium PairHMM prepare 占总 compute <=20%；
- spill 时 manifest 列出 spill bytes、run count、scratch path、CRC 状态。

## 5. I/O、存储和网络

### 5.1 本地文件

HTSlib/BGZF/CRAM decode、FAI/reference seek、CIGAR projection 留 CPU；按 coordinate/
interval 合并相邻 range，形成 batch 后进入 Kokkos。BGZF block 不是 read 边界，不能
按压缩块直接 scatter。CRAM 必须验证 reference、codec、index、dictionary。

writer 在 host 顺序完成 header、program record、sort、BGZF/BCF、tbi/csi；GPU 只返回
固定数值。关闭 writer 前验证文件、header、contig、排序、index、sidecar，最后写
OutputManifest complete。decode、prepare、compute、encode/index 分开计时。

### 5.2 远程对象和网络

先记录 request count、range size、first-byte latency、retries、compressed/decompressed
bytes、cache hit ratio；合并相邻 ranges，限制 in-flight requests/prefetch bytes；在
SLURM_TMPDIR 做有 byte/inode/TTL 的 cache；reference/index 使用锁和校验。reader 与
compute queue 解耦，backpressure 按 bytes 触发，结果按 shard/output ordinal 重排。

网络门禁：medium fixture 的 p50/p95 latency、cache hit、throughput 完整；warm-cache
decode >= local 0.80×；in-flight bytes <= host target 10%；transient retry、checksum/
index mismatch 产生稳定错误码，不生成 partial-success manifest。

## 6. SLURM / Nextflow

### 6.1 allocation 合同

~~~text
srun --cpus-per-task=task.cpus --mem=task.memory --tmp=task.disk fastgatk ...
srun --cpus-per-task=task.cpus --gpus=task.gpus --mem=task.memory fastgatk ...
~~~

native 校验 workers == SLURM_CPUS_PER_TASK（或显式 threads 不得超过 allocation）；
OMP/Kokkos/cpuset 一致；GPU visibility 与 registry 一致；memory.max 和 SLURM_TMPDIR
可写且预算可见；禁止工具内部再 sbatch/srun。

scatter 以 contig/interval/core 区间划分，halo 只用于计算；shard 携带 core、halo、
output_ordinal、schema/version。gather 按 contig/coordinate/ordinal deterministic merge，
不能按完成顺序拼接。OOM、scratch、transient network 可重试；输入/header/index/
strict mismatch 不应盲目重试。Nextflow cache key 包含 backend/ISA、compiler/Kokkos、
input/reference/index hash、parameters、determinism mode。

当前本地 smoke：

~~~text
bash fastgatk-native/scripts/verify_pipeline_local.sh
~~~

它验证 aliases、VCF/manifest、Nextflow process 等价路径和 SLURM_CPUS_PER_TASK/
memory/tmp 模拟；安装 Nextflow 后额外执行真实 DSL2 local process。真实 SLURM 用
FASTGATK_USE_SBATCH=1 验证 workflow/slurm_smoke.sh。

## 7. Benchmark 设计

### 7.1 JSON schema

每次 benchmark 至少保存 schema_version、tool、backend、strict/fast、hardware（CPU、
ISA、GPU、NUMA）、software（compiler、Kokkos、Git SHA）、resources（cpus/memory/GPU）、
input fixture hash、command/affinity、warmup/repeats、timing（prepare/H2D/kernel/D2H/
I/O/total）、memory（host/device peak/spill）、throughput、checksum/oracle、p50/p95。

### 7.2 分层

| 层次 | 输入/对照 | 指标 | 目的 |
|---|---|---|---|
| K0 correctness | fixed synthetic；GATK Java strict | raw bits、NaN/Inf、排序 | 浮点/算法错误 |
| K1 kernel | PairHMM、SW、BQSR、k-mer SoA；scalar Kokkos | cells/s、obs/s、active lanes、occupancy | SIMD/team kernel |
| K2 preparation | real batches；pool on/off | packing、allocation、H2D/D2H | amortization |
| K3 file boundary | BAM/CRAM/reference fixture | decode、seek、cache、VCF/index、RSS | I/O 合同 |
| K4 tool | HC/Mutect2/BQSR | wall、VCF/GVCF/stats、memory | 真实 tool |
| K5 pipeline | scatter/gather/resume/retry + SLURM | makespan、utilization、network/storage | 替换是否成立 |

数据集至少 small（512 synthetic）、medium（69k BAM、CRAM+reference）、large（真实
高/低 coverage、多 contig、长 read/indel、remote/cache-cold）。每个 fixture 记录
SHA256、reference/index、interval、filters。

公平性：相同 physical CPU set、输入、线程、warmup、迭代；GPU 显式计 H2D/D2H；Java
pure/GKL/native strict/fast/GPU 分开；GKL matrix 按 n*n 实际工作量校正；同时报告
kernel-only、amortized、file-to-file total、median/p95、temperature/power。

## 8. 性能与正确性门禁

### 8.1 阻断门禁

1. generic/AVX2/AVX-512 在 x86 编译；ARM/SVE、CUDA/HIP/SYCL 在对应 runner compile
   smoke；artifact metadata 正确。
2. PairHMM strict raw-bit=0；SW 比较 score/start/end/CIGAR/tie-break；BQSR/graph/
   region 比较 integer bins、排序、checksum。不通过不能进入 fast benchmark。
3. BAM/CRAM 可读回；VCF/GVCF header、contig、record、index、stats/sidecar 完整；
   manifest 只在 writer/index/validation 成功后 complete。
4. 无 SIGILL、OOM、device illegal access、data race、unbounded queue、partial-success；
   超预算返回稳定错误码或受控 fallback。
5. local Nextflow/SLURM smoke 和真实 scatter/gather/resume/retry 保持 output ordinal/
   manifest 一致。

### 8.2 性能 gate

| 指标 | PR/no-regression | 发布目标 |
|---|---:|---:|
| strict kernel | 当前 scalar/Kokkos baseline 的 0.95× | AVX2/AVX-512 在支持输入上 >= scalar 的 1.5×/2.0× |
| PairHMM vs GKL actual matrix | >=0.95× | steady-state >=1.10×，逐个报告 1/2/4/8/16 核 |
| PairHMM amortized | >=0.90× | workspace pool 后 medium >=1.00× GKL |
| regular DP scaling | 8 核效率 >=0.60 | SW/PairHMM 8 核 >=0.70，16 核 >=0.55 |
| irregular/atomic module | pinned baseline 的 >=0.90× | local reduction/bucket 后 8 核 >=0.50 |
| GPU offload | 无 GPU 不设数字 | 含 copy 的 medium/large total >= CPU 的 1.5×，否则 CPU |
| host/device memory | 无 OOM | host <=0.80 allocation；device <=0.85 usable memory |
| I/O | 正确输出 | warm-cache decode/encode >= local 0.80×且 telemetry 完整 |

这些是阶段性 gate，不是把小 fixture 数字外推到所有机器。硬件、compiler、Kokkos 或
数据分布变化都要产生新 baseline，gate 输入必须带 baseline hash 和环境。

### 8.3 fast 数值门禁

PairHMM 先定义 per-likelihood abs/rel tolerance，再比较 genotype、allele、DP/AD/PL/
GQ 和 VCF record。GKL 相对纯 Java 已有约 2.61e-8 级差异，因此 fast/GPU 不要求所有
路径 byte-identical；但 tolerance、FTZ/DAZ 和 output fields 必须写入 manifest。
VCF/GVCF record order、sample/allele order、END/block、header、index 必须 exact。

## 9. 分阶段交付

### P0：基线冻结（当前→1 周）

交付 hardware/software probe、统一 benchmark JSON、现有 PairHMM/module/HTS 结果归档、
strict oracle CI、perf gate checker。退出：generic/AVX2/AVX-512 JSON 可复现，oracle、
BAM/CRAM smoke、affinity、memory limit、Kokkos/hash 齐全。

### P1：CPU production path（1–3 周）

交付 persistent KernelPlan/arena、length bucket、PairHMM prefix cache、CPU registry、
NUMA/SLURM probe、prepare/H2D/D2H/total telemetry。退出：medium amortized >=1.0× GKL、
RSS gate、strict 不退化、AVX2/AVX-512 frequency/power report。

### P2：算法 kernel breadth（3–6 周）

交付 SW score+CIGAR/tie-break、BQSR local histogram/deterministic merge、k-mer
encode/sort/RLE、HC assembly→likelihood→genotype batch API。退出：K0–K3 通过、每个
模块有 file-boundary/memory/IO profile、未支持语义明确回退 Java。

### P3：GPU proof-of-value（6–10 周）

交付一个 CUDA 和一个 HIP 或 SYCL runner；PairHMM team/wavefront、SW tile/checkpoint、
async pipeline、device arena、threshold controller、CPU fallback。退出：真实 device
test，含 H2D/D2H 的 medium/large 达到目标才默认启用，否则 registry 保持 experimental。

### P4：端到端 GATK/SLURM/Nextflow（10–16 周）

交付 HC VCF/GVCF/index/manifest、GATK oracle、Nextflow scatter/gather/resume/retry、
SLURM CPU/GPU profiles、remote/cache-cold I/O、性能 dashboard。退出：K4/K5 通过，
真实样本 variant/VCF/sidecar 契约通过；此时才能讨论直接替换 GATK。

## 10. 下一步命令

~~~text
bash pairhmm-demo/scripts/build_kokkos_variants.sh
python3 pairhmm-demo/scripts/verify_kokkos.py --records=512 --workload=matrix8 --variants=scalar,avx2,avx512
bash kokkos-modules-demo/scripts/build_variants.sh
python3 kokkos-modules-demo/scripts/verify_modules.py --records=512 --threads=4
python3 pairhmm-demo/scripts/benchmark_kokkos.py --variants=scalar,avx2,avx512 --cores=1,2,4,8,16 --repeats=5 --iterations=10 --output=pairhmm-demo/results/benchmark-perf-baseline.json
python3 pairhmm-demo/scripts/benchmark_same_cores.py --native-backend=avx2 --cores=1,2,4,8,16 --repeats=5 --iterations=10 --output=pairhmm-demo/results/benchmark-same-cores-avx2.json
bash fastgatk-native/scripts/build_native.sh
python3 fastgatk-native/scripts/verify_native.py
bash fastgatk-native/scripts/verify_pipeline_local.sh
~~~

下一段实现应优先做 P0/P1 的 persistent arena、统一 telemetry 和 gate checker；在这
之前继续优化单个 SIMD kernel，会把真正的 prepare/I/O/内存瓶颈隐藏起来。
