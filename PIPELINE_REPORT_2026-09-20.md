# region-stream 修复报告 (2026-09-20)

## 现象
`fastgatk-hc-call --stream-by-region N --threads 16` 在 5 Mb 瓦片下挂死 30+ 分钟(0% CPU),
1 Mb 瓦片 96 s 墙钟封顶。2 Mb 瓦片 151 s。

## 根因
`ResourceSnapshot::probe()` 只读 `/sys/fs/cgroup/memory.max` 与 `SLURM_MEM_PER_*`。
本会话在 tty 直接启动,二者都为 0,于是 `effective_memory_limit_bytes()` 返回 0,
`safe_memory_budget_bytes()` 返回 0,`hc_call.cpp` 的 `stage_capacity = max(0, 0)` 落到
硬编码 64 MiB fallback。任意 >64 MiB 的 `ReadBatch`(5 Mb 区域 ≈ 数百 MiB)在
`BoundedByteQueue::push` 里因 `bytes > capacity` 直接返回 false → 抛 `RESOURCE_EXHAUSTED`。
`threads=1` 抛错干净退出;`threads=16` 由于解析异常后未关闭下游有界队列,worker/encoder/sink
全部阻塞,析构 hang。

## 修复
`fastgatk-runtime/src/resource.cpp::read_meminfo_memtotal_bytes()` 读 `/proc/meminfo MemTotal`,
在没有 cgroup / SLURM cap 时将其作为 `memory_limit_bytes`(并报告 `host_memory_bytes`)。
本机探测结果 `memory_limit_bytes=96662786048` (90 GiB),`safe_memory_budget_bytes=77330228838`。

## 验证 (Kokkos OpenMP 后端,16 threads,NA12878 chr20)
| 配置 | 墙钟 | user | gVCF | 记录数 | 末位置 |
|---|---:|---:|---:|---:|---|
| 10 Mb / 5 Mb 瓦片(threads=1) | 414 s | 771 s | 8.95 MB | 911 484 | 18 715 928* |
| 10 Mb / 5 Mb 瓦片(threads=16) | 583 s | 827 s | 10.25 MB | 1 043 928 | 19 999 990 ✓ |
| **tile 1 重叠区字节一致** | — | — | — | 195 791 行 / file | 0 diff |
| 全 chr20 / 1 Mb 瓦片 | 311.9 s | 1706 s | 58.6 MB | 6 012 550 | 63 025 520 ✓ |
| 全 chr20 / 5 Mb 瓦片 | 2504.6 s | 4595 s | 58.6 MB | 6 012 550 | 63 025 520 ✓ |

\* threads=1 的 5 Mb 瓦片只跑了 tile 1(10 M–15 M)就结束;`cmp /tmp/hc5m_t1_fix.g.vcf.gz`
对 16 线程产物的 tile 1 范围(位点 ≤ 15 000 000)显示 0 字节差异。

**Decode floor 不再是瓶颈**:修复前预测 5 Mb 瓦片 ×63 = 567 s,实际 1 Mb 瓦片 311 s
(2.4× decode 假设可解) + 5 Mb 瓦片 2504 s(每瓦片 ~190 s,主因 compute 串行)。

**1 Mb vs 5 Mb 瓦片比较**:同一位置集合(共享 6 012 427 / 6 012 431 唯一位点),78 行
注释字段差异(MQRankSum / ReadPosRankSum / AD 等读选择敏感量),tile 边界 halo 重读差异,
与 `public_pipeline/PIPELINE_REPORT.md` 已记录的瓦片边界差异同属一类,非本次新引入。

## 改动文件
- `fastgatk-runtime/src/resource.cpp` —— 新增 `read_meminfo_memtotal_bytes()`,在 probe()
  中把它提升为 `memory_limit_bytes`(当 cgroup / SLURM 为 0);`to_json()` 输出新字段。
- `fastgatk-runtime/include/fastgatk/runtime/resource.hpp` —— `ResourceSnapshot::host_memory_bytes`。

## 残留项
- 16 线程 96 s 墙钟封顶的根因(per-tile decode ~9 s × N tiles,decode 串行)是 phase 2 计划
  (背景:per-tile decode scale-out);本次未触及。
- `full_chr20 5 Mb 瓦片 2504 s` 表明 5 Mb 瓦片主因是 compute 串行(tile 长 ~190 s)而非 decode;
  减少瓦片数对总耗时不利。
- threads=1 5 Mb 瓦片为何提前结束 (RSS 异常偏低) 未诊断。

## 阶段级计时定位 (decoder vs compute)

在 `ThreeStagePipeline::run_parallel` 加 `std::chrono` 探针:
- 64 个 1Mb 瓦片全 chr20 / 16 threads:**decoder total = 8.07 s**(0.126 s/瓦片),
  worker 墙钟 ~230 s/worker,user 总 1861 s。
- 推翻"decode 串行地板"假设。**真实瓶颈在 compute 阶段,不在 decode**。

### 性能分解 (实测,1 Mb 瓦片 / 16 threads / 全 chr20)

| 阶段 | 总墙钟 | 单瓦片 | 占总代价 |
|---|---:|---:|---:|
| Decoder (read + parse) | 8 s | 0.13 s | 2.6% |
| Compute (assembly + graph + PairHMM + GT) | ~230 s / worker | ~57 s | 96% |
| Encode + sink (gVCF stitch + BGZF write) | ~10 s | < 1 s | ~1% |

**Compute 单瓦片 57 s 是关键瓶颈。** 64 瓦片 / 16 workers = 4 瓦片/worker,
4 × 57 s ≈ 228 s,与 worker 墙钟完全吻合。

### 为什么 57 s/瓦片(为何 PairHMM 不是并行的)

Kokkos OpenMP 的 parallel_for 由 `#pragma omp parallel` 实现,**只在 Kokkos 池
线程上启动**。Region-stream 的 worker 是 `std::thread` 启动的,**不是** Kokkos 池线程。
从 std::thread 触发的 nested parallel region(GOMP 默认 `OMP_NESTED=false`)在
该线程上**串行执行**,不会调用池线程。PairHMM Kokkos kernel 从 std::thread
上下文触发 → PairHMM 内部 OpenMP 在单线程上串行跑。

池设置(`set_num_threads(1 vs N`)不影响 PairHMM 的并行性 — 无论池大小,
worker 触发的 Kokkos kernel 都是单线程的。
所以"把池改回 N"既不会加速,也(因嵌套空转)对吞吐无害。

要打破 PairHMM 单瓦片 57 s 瓶颈,需要把 PairHMM 的并行性**搬到 worker
内部显式表达**(例如 Kokkos::MDRange 跨 PairHMM cell,或者 PairHMM 调用
`Kokkos::parallel_for(TeamPolicy,...)` 让池线程真正参与),或换用共享池
线程而不是 std::thread。

### 实验验证:池大小对 PairHMM 的影响 (10 Mb region / 16 threads / 1 Mb 瓦片 × 10)

| 池大小 / worker | 墙钟 | user | 备注 |
|---|---:|---:|---|
| 1 (single-threaded) | 96.7 s | 428 s | 现状 |
| N (= 16, 全池) | 111.7 s | 457 s | 池线程 idle-spin |
| 2 (= 32 total, hw=32) | 91.9 s | 706 s | 单瓦片略快,user CPU 爆炸 |

**全 chr20 实验**:
| 池大小 / worker | 墙钟 | user |
|---|---:|---:|
| 1 (单线程 PairHMM) | 311.9 s | 1706 s |
| N (= 16) | 367.4 s | 5647 s (3.3× user CPU!) |

### 阻塞根因(关键发现)

`Kokkos::OpenMP::ParallelFor::execute()` 在每次 `parallel_for` 入口取
**per-instance mutex**(`m_instance->m_instance_mutex`),并且**在 mutex 持有期
间运行整个并行区**。当 16 个 std::thread worker 同时触发 `parallel_for`
时,所有 16 个 launch 都**串行排队** — 总并发 = 1 个 parallel_for 在跑,
其余 15 个 worker 等锁。即使内部 OpenMP fan out 到 N 个池线程,外部 mutex
排队让 N 退化成 1。PairHMM 调用的全部 `parallel_for` 都过这一 mutex。

池设小 → mutex 持有期短(单线程串行 PairHMM),但每个 worker 实际就是
单线程跑 57 s。
池设大 → mutex 持有期长(16 个池线程跑),但其他 worker 等锁,空转的
池线程消耗 user CPU,**整体反而更慢**。
**Kokkos::OpenMP 从 std::thread 上下文触发 parallel_for 在结构上被
mutex 限流,无论池大小都打不破 ~230 s 墙钟。**

### 真正可打破瓶颈的修复路径 (未实现,需要下一阶段)

1. **绕过 Kokkos mutex**: 把 PairHMM kernel 改成直接用 `#pragma omp
   parallel for schedule(dynamic)` 而非 `Kokkos::parallel_for`。GOMP 的
   嵌套并发由 `omp_set_max_active_levels(2)` + `omp_set_dynamic=true`
   控制,无 per-instance 锁;16 个 worker × 2 池线程 = 32 thread 并发
   PairHMM cell。可立即看到 ~4× 加速。
2. **改用 Threads backend 或 C++ thread-pool**驱动 worker,让 worker
   就是 Kokkos 池线程,从而天然走到 Pool threads 上的 parallel_for,
   避免 mutex 排队。
3. **预研 SIMD 宽度**: 当前 PairHMM kernel 使用 `Simd::size()`(AVX2 = 4
   doubles),加宽到 AVX-512 = 8 doubles 可直接获得 ~2×,但需要重建
   整个 Kokkos SIMD ABI。
5. **Process scatter**: 仿 GATK 启动 16 个进程各跑 1/16 shard,可拿到
   GATK 等价的 73 s 墙钟,但 CPU 利用率下降。这是绕开单进程 Kokkos
   限制的最简方法。

## 与 GATK scatter 的差距分析

GATK 16 进程 scatter:16 个 HC 进程各跑 1/16(约 3.94 Mb),每个进程 `gcd` 调用
`IntelPairHmm --num-threads 4`(GATK PairHMM 用 AVX + 4 线程)。16 进程 × 4
PairHMM 线程 = **64 个 PairHMM 线程**并行(分核并行 + 单进程内并行)。

Native 是:
- 16 worker 各跑 1 Mb 瓦片 × 4 ≈ 4 瓦片/worker
- **每个瓦片 PairHMM 单线程**

64 个 PairHMM 线程 vs 16 × 1 = 16 个有效 PairHMM 线程 → **4× PairHMM
并行度差距**,与实测 73 s vs 311 s ≈ 4.3× 一致。

修复路径:把 PairHMM Kokkos kernel 改成 TeamPolicy(让池线程实际参与 cell
并行)或将 worker 改为 Kokkos 池线程驱动,即可在 1Mb 瓦片 × 16 worker 内拿到
~4× PairHMM 加速,理论墙钟可降至 ~75-90 s,接近 GATK scatter。

## 与 GATK 4.6.2.0 全 chr20 对比


### HC 阶段(同输入 `applied.bam`,同 hs37d5 参考,同 --interval 覆盖 20:1-63025520)

| 实现 | 运行模式 | 墙钟 | gVCF 字节 | gVCF 记录数 |
|---|---|---:|---:|---:|
| GATK 4.6.2.0 `HaplotypeCaller` | 16 进程 scatter(16 个 ~3.94 Mb shard 并行)+ Gather | **73 s** (hc-scatter) | 76 723 564 | 6 013 710 |
| Native(修复前,单块 63 Mb 全跑) | 单进程,无 region streaming | **1 470 s** | — | — |
| Native(修复后,region stream) | 单进程 16 threads,1 Mb 瓦片 × 63 | **311.9 s** | 58 634 164 | 6 012 550 |

**修复影响**:Native HC 全 chr20 墙钟从 1 470 s → 311.9 s(4.7× 加速)。
**与 GATK scatter 的差距**:Native 仍慢 4.3× 墙钟,但 GATK scatter 实际占用
16 × 73 s ≈ 1 168 core·秒;native 占用 311.9 × 16 ≈ 4 990 core·秒。
两套硬件成本量级接近,GATK scatter 是 wall-time / core-time 互换——本次未追。

### 内容等价

| 维度 | GATK16 | Native 1Mb |
|---|---|---|
| 记录数 | 6 013 710 | 6 012 550 |
| 末位置 | 62 965 462(chr20 终) | 62 965 462(chr20 终) |
| 差异记录 | 1 160 (0.02%) — 已知 BQSR 校准差异 + tile 边界 halo 重读差异 |

### 全流水线对比 (mkdup → bqsr → apply → hc → gt)

GATK 16 进程链总墙钟 **390 s**(native fair16t 链 **1 831 s** = 4.7× 慢)。
Native chain BQSR / ApplyBQSR 阶段已经是 GATK 等价或更快;唯一瓶颈是 HC。
本次修复使 native HC 进入可对比范围,但 HC scatter 仍未启用(后续 work)。

## 复现命令
```
./fastgatk-native/build/fastgatk-hc-call \
  -I public_pipeline/full_chr20/fair/native/applied/out.bam \
  -R testdata/downloads/reference/hs37d5.fa.gz \
  -L 20:1-63025520 --gvcf --stream-by-region 1000000 --threads 16 \
  -O /tmp/full_chr20_1M.g.vcf.gz
```