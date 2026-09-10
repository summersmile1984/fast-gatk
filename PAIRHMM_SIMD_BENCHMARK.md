# PairHMM native/SIMD 示例：GATK 兼容性与性能实测

日期：2026-08-21  
源码：`/home/turing-agents/Documents/fast-gatk/fastgatk-kernels/src/pairhmm_kokkos.cpp`（demo 仅保留客户端、oracle 和 benchmark）
参考：GATK 4.6.2.0 官方 package，纯 Java `LoglessPairHMM`；另测试 Intel GKL AVX 路径。

## 当前结论：统一 Kokkos API 实现（取代下文 raw-intrinsics 原型）

当前主实现是 `fastgatk-kernels/src/pairhmm_kokkos.cpp`。算法、内存和调度只使用 Kokkos 5.2.0 API：

- `Kokkos::View` 管理输入、pair ID、workspace 和结果，`deep_copy` 负责 memory-space 传输；
- `Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>` 分发 SIMD group；
- `Kokkos::Experimental::simd<double>` 负责 scalar/AVX2/AVX-512 lane；源码中没有 `_mm*`；
- 同一源码以 generic、`Kokkos_ARCH_ZEN3`、`Kokkos_ARCH_ZEN4` 构建，实测 SIMD width 为 1/4/8；
- read/haplotype ID 若连续则用 Kokkos SIMD contiguous load，矩阵组内相同 read 用 broadcast，任意 ID 退回 portable generator gather；
- M/I/D 从六行滚动 workspace 优化成三行原位更新：被覆盖的上一行 cell 暂存在 SIMD 寄存器，recurrence 和加法顺序不变；
- OpenJDK fdlibm `log/log10` 运算顺序内置为 device-callable strict 路径，避免结果依赖 host `libm`。

固定版本：Kokkos 5.2.0、CMake 4.3.4；构建入口 `pairhmm-demo/scripts/build_kokkos_variants.sh`。

### 严格一致性

`python3 pairhmm-demo/scripts/verify_kokkos.py` 对真实 GATK 4.6.2.0 Java `LoglessPairHMM` 的 512 records、8×8 blocks、4096 个交叉 likelihood 做 raw IEEE-754 比较：

| Kokkos build | SIMD width | bit_different |
|---|---:|---:|
| generic scalar | 1 | 0 / 4096 |
| ZEN3 / AVX2 | 4 | 0 / 4096 |
| ZEN4 / AVX-512 | 8 | 0 / 4096 |

这里的 oracle 固定 `-XX:DisableIntrinsic=_dlog10`，使 GATK Java 使用 StrictMath/fdlibm。原因不是 FP32：全部计算都是 FP64，且 4096 个 DP scaled sum 在默认 HotSpot 下也全部逐位相同；默认 HotSpot x86 会把 `Math.log10` 换成 Intel LIBM intrinsic，最终有 3/4096 个 likelihood 相差 1 ULP。Java `Math` 本身允许这种平台实现差异，所以“任意默认 JVM 上都相同”不是稳定合同。项目的 strict 合同必须固定 GATK/JDK/VM flags；fast/GKL 合同继续用 tolerance 和最终 VCF oracle。

### 最终 benchmark 方法

- 主机：AMD Ryzen 9 7945HX，16 physical / 32 logical cores，支持 AVX2 和 AVX-512；
- 两边绑定同一物理 CPU 集 `0..N-1`，不使用 SMT sibling；
- 512 records，read 150、haplotype 160、seed 42；
- independent 30 iterations，matrix8 10 iterations，完整 warmup 1 次；
- 重复 3 次取 median，并交替 baseline/Kokkos 顺序；
- kernel rate 不含 input packing/deep copy；amortized rate 把一次准备成本计入整次 invocation；
- GKL 按 8×8 block 实际 `n*n` 工作量校正，不再拿 diagonal record count 比较。

结果文件：[benchmark-kokkos-5.2.0-range-final.json](/home/turing-agents/Documents/fast-gatk/pairhmm-demo/results/benchmark-kokkos-5.2.0-range-final.json)。

### Independent pairs：Kokkos AVX-512 vs GATK pure Java

| 物理核 | GATK Java pair/s | Kokkos kernel pair/s | kernel/Java | Kokkos amortized pair/s | amortized/Java |
|---:|---:|---:|---:|---:|---:|
| 1 | 11,081 | 76,297 | 6.89× | 70,639 | 6.37× |
| 2 | 20,993 | 154,489 | 7.36× | 133,252 | 6.35× |
| 4 | 40,736 | 304,627 | 7.48× | 236,664 | 5.81× |
| 8 | 74,153 | 567,531 | 7.65× | 374,997 | 5.06× |
| 16 | 120,928 | 944,452 | 7.81× | 495,517 | 4.10× |

### Matrix8：Kokkos AVX-512 vs GKL AVX-OMP 的实际矩阵工作量

| 物理核 | GKL actual pair/s | Kokkos kernel pair/s | kernel/GKL | Kokkos amortized pair/s | amortized/GKL |
|---:|---:|---:|---:|---:|---:|
| 1 | 70,432 | 76,874 | 1.09× | 73,880 | 1.05× |
| 2 | 137,800 | 152,374 | 1.11× | 141,109 | 1.02× |
| 4 | 260,639 | 296,841 | 1.14× | 263,757 | 1.01× |
| 8 | 463,772 | 553,893 | 1.19× | 461,023 | 0.99× |
| 16 | 724,019 | 924,297 | 1.28× | 681,994 | 0.94× |

因此可以下两个不同层次的结论：steady-state Kokkos AVX-512 kernel 已在本机 1–16 核超过 GKL 9%–28%；若每十次矩阵计算都重新分配、打包和 deep-copy，一次性端到端优势在高核数会被准备阶段吃掉。生产 runtime 必须池化 workspace、复用已转换 read qualities，并把多个 region/batch 放进持久 Kokkos 实例；当前结果仍只是 PairHMM kernel/API slice，不代表整个 HaplotypeCaller 已超过 GATK。

下文 round-5 raw-intrinsics 数字保留为优化历史，不再代表当前后端选择。

## 实现边界

本例实现 GATK `LoglessPairHMM` 的以下合同：

- 矩阵为 `(readLength + 1) × (haplotypeLength + 1)`；
- 第 0 行 deletion 初始化为 `2^1020 / haplotypeLength`；
- match/insertion/deletion 三个 recurrence 的计算顺序与 Java 源码相同；
- final likelihood 只累加最后一行的 match + insertion，忽略 deletion；
- 默认 tristate correction 为 3；
- 使用 Java GATK 运行时导出的质量概率和 match-to-match 表，避免 host `libm` 的最后几位差异。

实现分为两条路径：

| 路径 | 用途 | 一致性 |
|---|---|---|
| native scalar | strict compatibility reference | 与 GATK pure Java raw IEEE-754 bit-exact |
| native AVX2 | 性能路径，4 个独立 pair/批次 | 每 lane 保持 scalar recurrence；当前等长回归集与 strict raw bit-exact |
| native AVX-512 | 性能路径，8 个独立 pair/批次 | 每 lane 保持 scalar recurrence；当前等长回归集与 strict raw bit-exact |

SIMD 沿独立 read/haplotype pair 展开，而不是重排单个动态规划矩阵；这样不会改变 PairHMM 的依赖拓扑和单 lane 累加顺序。`-ffp-contract=off -fno-fast-math` 禁止编译器用 FMA/fast-math 改写浮点合同。

## 本机和固定条件

- CPU：AMD Ryzen 9 7945HX
- 16 physical cores / 32 logical CPUs，单 NUMA node
- CPU flags 含 AVX2、AVX512F；Kokkos 是同源码的 scalar、ZEN3/AVX2、ZEN4/AVX-512 构建，launcher 按 host flags 和现有 binary 安全选择，不在不支持的节点启动高 ISA binary
- 两边都用 `taskset` 绑定到相同物理核：`0`、`0-1`、`0-3`、`0-7`、`0-15`
- 输入：512 个 pair，read 150 bp，haplotype 160 bp，seed 42
- 每次测量 10 iterations，重复 3 次取 median
- Java 侧是真实 GATK 4.6.2.0 `LoglessPairHMM`，不是自写 Java 参考
- native 侧加载 Java 导出的 `gatk-tables.hex`，保证 strict 常量一致

本次实验资产 SHA-256：GATK package `32d2f90bf13fcb3a8ac765bb2cb8ec1fc9a6cc447055d0156bd1db2092d4e3e8`；JDK 17 `be7668bc030d578b83d6d5ef9221d6d6729bbbca8cf94a7d52e16ac68b5a5a35`；Java 导出的表 `84e1ce7a3e9dcfc0b6a21ade7f4e292067389297a0f955954bd74134a6e2010e`。

运行命令（分别固定 backend，避免把两种实现混为一个数字）：

```bash
python3 pairhmm-demo/scripts/benchmark_same_cores.py \
  --native-backend=avx2 \
  --repeats=5 --iterations=10 \
  --output=pairhmm-demo/results/benchmark-round5-avx2.json

python3 pairhmm-demo/scripts/benchmark_same_cores.py \
  --native-backend=avx512 \
  --repeats=5 --iterations=10 \
  --output=pairhmm-demo/results/benchmark-round5-avx512.json
```

Java 和 native 都先完整 warmup 一轮，首次 JIT、线程启动和 workspace 分配不计入 kernel 时间。早期未 warmup 的多线程数字不再作为性能结论。

## 结果：GATK pure Java vs native AVX2（round 5，steady-state）

| 相同物理核心数 | GATK Java pairs/s | native AVX2 pairs/s | native/GATK |
|---:|---:|---:|---:|
| 1 | 10,974 | 14,498 | 1.32× |
| 2 | 20,621 | 28,026 | 1.36× |
| 4 | 37,484 | 44,860 | 1.20× |
| 8 | 76,515 | 89,805 | 1.17× |
| 16 | 110,386 | 150,099 | 1.36× |

结果文件：[benchmark-round5-avx2.json](/home/turing-agents/Documents/fast-gatk/pairhmm-demo/results/benchmark-round5-avx2.json)。本次所有 native checksum 为 `-22478.5810851`。

## 结果：GATK pure Java vs native AVX-512（round 5，steady-state）

| 相同物理核心数 | GATK Java pairs/s | native AVX-512 pairs/s | native/GATK |
|---:|---:|---:|---:|
| 1 | 10,865 | 10,092 | 0.93× |
| 2 | 18,508 | 19,792 | 1.07× |
| 4 | 40,881 | 31,863 | 0.78× |
| 8 | 75,087 | 61,267 | 0.82× |
| 16 | 103,804 | 93,229 | 0.90× |

结果文件：[benchmark-round5-avx512.json](/home/turing-agents/Documents/fast-gatk/pairhmm-demo/results/benchmark-round5-avx512.json)。`auto` 不根据历史单次吞吐硬编码选择 ISA，而是按 host flags 和 binary 可用性选择 AVX-512 → AVX2 → scalar；需要比较不同构建时显式指定 backend。

### 本轮优化的实际作用

* SIMD workspace 从三张完整矩阵改为两行滚动 DP。150×160 时，AVX-512 每线程约从 4.7 MB 降到 62 KB，显著降低 8/16 核共享缓存压力。
* `2^1020 / hapLength` 从每个 haplotype cell 重复除法改为每 lane 一次。
* read quality 对应的 match/mismatch prior 从内层 haplotype 循环提到每个 read row 一次，同时保持逐位相同的运算顺序。
* Java/native 都加入一轮不计时 warmup，消除早期 benchmark 对 Java 多线程的不公平惩罚。

最后一项也使数字看起来没有早期版本“惊艳”，但它更接近生产 steady-state。当前可信结论是 AVX2 kernel 比 pure-Java `LoglessPairHMM` 快约 17%–36%，而不是多倍数量级。round 5 还交替 Java/native 的运行顺序，降低 boost/温度造成的系统性偏差。

### 兼容性证据

512 个 pair 的 strict native 与 GATK Java raw bits 比较：

```text
bit_different = 0
max_abs        = 0.0
all_within_1e-12 = true
```

对应命令：

```bash
taskset -c 0 ./pairhmm-demo/pairhmm-demo --mode=scalar \
  --pairs=512 --read-len=150 --hap-len=160 --threads=1 --iterations=3 \
  --seed=42 --values-out=pairhmm-demo/results/native-medium-strict.hex

taskset -c 0 third_party/jdk17/bin/java \
  -cp pairhmm-demo/java-classes:third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
  org.broadinstitute.hellbender.utils.pairhmm.GatkPairHmmBenchmark \
  --mode=java --pairs=512 --read-len=150 --hap-len=160 --threads=1 --iterations=3 \
  --seed=42 --values-out=pairhmm-demo/results/gatk-medium-strict.hex

python3 pairhmm-demo/scripts/compare_pairhmm.py \
  pairhmm-demo/results/native-medium-strict.hex \
  pairhmm-demo/results/gatk-medium-strict.hex
```

当前 `fastgatk-kernels`/Kokkos 构建会默认加载同一份 pinned Java 表；只有
package 将表移到非 source-tree 位置时，才需要设置 `FAST_GATK_PAIRHMM_TABLES`
或在 CMake 中指定 `FASTGATK_PAIRHMM_TABLES_FILE`。

本轮修正了 SIMD 中 `A+B+C` 的加法顺序，显式保持 Java/C++ 标量的 `(A+B)+C`。在 512 个 pair、150/160、seed 42 的当前回归集上：

```text
scalar-strict vs GATK Java: bit_different = 0, max_abs = 0.0
AVX2          vs GATK Java: bit_different = 0, max_abs = 0.0
AVX-512       vs GATK Java: bit_different = 0, max_abs = 0.0
```

另外用 `(64,37,41,seed=1)`、`(73,99,113,seed=7)`、`(128,150,160,seed=123)` 三组不同尺寸/seed 重复验证，scalar、AVX2、AVX-512 均为 `bit_different=0`。

这证明“跨独立 pair 的 SIMD”可以在受控的 FP 环境和固定运算顺序下做到逐位一致；但仍不能把单一 x86/toolchain 回归结果外推成所有 CPU/libm 的永久保证。因此 runtime 仍提供：

```text
FAST_GATK_MODE=deterministic  -> native scalar / strict golden
FAST_GATK_MODE=fast           -> runtime AVX2/AVX-512，声明 tolerance 并持续做 raw-bit 回归
```

严格模式固定 MXCSR 为 round-to-nearest、关闭 FTZ/DAZ；GKL 的 fast 路径则主动启用 FTZ，二者数值契约必须分开。

## GATK Intel GKL AVX 的单独结果

官方 GKL 使用批量矩阵 API。为了给每个 pair 保持不同 haplotype，本例按 8×8 block 调用 GKL，并只统计对角线结果；GKL 实际计算了额外的交叉 pair，因此这个数字不能直接与独立 pair AVX2 数字做严格吞吐排名。

本轮 warmup 后的 effective diagonal pairs/s 为：

| 物理核 | GKL AVX-OMP diagonal pairs/s |
|---:|---:|
| 1 | 8,827 |
| 2 | 17,110 |
| 4 | 32,687 |
| 8 | 59,164 |
| 16 | 85,294 |

这里每个 8×8 block 实际计算 64 个 likelihood，只取 8 个对角输出；因此实际 DP pair 工作量约是表中数字的 8 倍。在真实 HaplotypeCaller 中 read×haplotype 交叉积正是所需工作，当前“独立 pair”native demo 还不能据此宣称超过 GKL。下一步公平 gate 必须实现相同的矩阵 API、复用 read/haplotype 数据，并以总 DP cell 数比较。

GKL 与纯 Java 的输出最大差约 `2.61e-8`，这是 GATK 自身硬件路径相对于纯 Java `LoglessPairHMM` 的数值差异；它应标为 bio-identical，而不是 byte-identical。下一步如果要公平比较 GKL，需要让 native 也采用同样的 read×haplotype block API，并以计算的总 cell 数作为吞吐单位。

## 这对总体框架的意义

1. **兼容层先行。** `native scalar` 是 PairHMM 的 reference backend，能作为 GATK golden 的严格门槛。
2. **SIMD 跨 pair。** 不改变单 pair DP 顺序，比在 haplotype 方向强行 vectorize 更容易维持 GATK 语义。
3. **资源合同显式化。** `threads` 与 `taskset` 同时控制；SLURM 中应将 `--cpus-per-task` 映射到 native worker 数，并避免额外线程池。
4. **性能结果分层。** pure Java、GATK GKL AVX、native strict、native fast 必须分开报告，不能混为一个“GATK 基线”。
5. **生产实现还需要。** 预分配 workspace、批次化 reads/haplotypes、NUMA/缓存布局、AVX-512 频率/功耗测量、真实 HaplotypeCaller 输入分布、跨 ARM/AVX512/GPU determinism 和端到端 I/O 测试。

当前结论只针对 PairHMM kernel，不代表 HaplotypeCaller 或完整 GATK pipeline 已经获得同等加速。

## 架构、AVX-512 与 bit-identical 决策

### Rust host 与 C++/Kokkos 边界

Kokkos 是 C++ 原生 API，因此本项目的推荐热路径是：

```text
GATK-compatible launcher / Nextflow
        -> C++ fast-gatk host/runtime
        -> Kokkos execution + memory space
        -> PairHMM/SW/k-mer kernels
```

Rust 调用一个 C ABI shim 本身通常只有很小的固定成本；真正危险的是把 FFI 放在 read、base 或 DP cell 粒度。若一次调用只提交一个 region/shard 的连续 batch，数据复制、调度和 kernel 时间会把这点成本摊薄；若每条 read 都跨一次边界，边界和分配成本会主导。故第一版建议 **C++ host + Kokkos core**，保留稳定的 C ABI 作为外部插件/测试边界；Rust 若保留，只做批次级 orchestrator 或独立 CLI，不进入内层循环。这样也能减少 Kokkos/CMake、CUDA/HIP/SYCL toolchain 与 Rust build-script 的组合复杂度。性能是否超过某个阈值仍应以批次级 FFI microbenchmark 验证，不能凭语言名称推断。

### 当前是否使用 Kokkos SIMD API

是。当前主路径的 allocation/copy、`View`、execution policy 和 SIMD 都是 Kokkos API；`pairhmm_kokkos.cpp` 不含 `<immintrin.h>` 或 `_mm*`。`Kokkos::Experimental::simd<double>` 的 ABI/宽度由 `KOKKOS_ARCH` 决定，因此 generic、AVX2、AVX-512 是同源码的不同构建产物，而不是算法内的 ISA 分支。

GPU 后端下 Kokkos 5.2 的默认 `simd<T>` 是 scalar ABI，所以这份源码可以编译成“一线程一 pair”的正确基线，但不能据此声称 GPU 已优化。GPU 性能版仍需用 Kokkos hierarchical/team policy 做 wavefront/tile，并继续复用相同 `PairIndexBatch`、`View` 和 strict oracle。当前只在本机 OpenMP/x86 三个构建上实际运行；ARM SVE/NEON、CUDA、HIP、SYCL 需要各自 CI/实机验证。

### 当前 CPU 与 AVX-512路线

本机是 AMD Ryzen 9 7945HX；`lscpu` 暴露 `avx512f/dq/bw/vl`、`avx512_vnni`、`avx512_bf16` 等标志，GATK GKL 实际运行日志也显示 `Using CPU-supported AVX-512 instructions`。Kokkos 构建产物为 generic、ZEN3/AVX2、ZEN4/AVX-512；same-core benchmark 在启动前读取 host flags 并检查对应 binary，按 AVX-512 → AVX2 → scalar 降级。显式请求不可用 ISA 会 fail-closed，因此不会把非法指令误报成性能结果。

当前 Kokkos AVX-512 一次处理 8 个独立 FP64 pair。在三行原位 DP 和 `RangePolicy` 优化后，它是本机可用的高吞吐 CPU backend；最终选择仍应以同一 workload 的实测结果为准，不能把某一次 AVX2/AVX-512 排名外推到所有节点。

### bit-identical 的目标和边界

这里不是 FP32 导致的差异：当前 native 和 GATK Java 都以 `double`/FP64 计算。差异来自 IEEE-754 的舍入路径，包括运算/归约顺序、是否融合 FMA、编译器重排、次正规数的 FTZ/DAZ 模式、`pow/log10` 实现，以及常量表是否逐位相同。GATK GKL 还明确打开 FTZ，并采用块状 SIMD 算法，所以它相对纯 Java 本来就不是 byte-identical，而是 bio-identical/tolerance-level。

当前结果已经证明两种契约可以分开：

* `scalar-strict` 使用 GATK Java 导出的逐位 quality/transition 表、固定递推顺序、`-fno-fast-math -ffp-contract=off`，512 个 pair 与真实 GATK Java 的 `bit_different=0`。
* AVX2/AVX-512 全部仍是 FP64；在当前固定输入集上通过了 raw-bit 比较。此前出现的 5 个差异来自 SIMD 写成 `A+(B+C)` 而标量是 `(A+B)+C`，现已修正。

因此“根 GATK”应先明确目标：若目标是替换纯 Java GATK，发布 `deterministic` 契约（严格模式必要时回退 scalar）；若目标是替换 GKL/FastestAvailable，发布 `fast` 契约并给出 tolerance。严格模式应固定表和常量、禁止 fast-math/FMA 重排、固定 MXCSR（关闭 FTZ/DAZ）、固定归约顺序，并在 CI 对 GATK golden 做 raw IEEE bits 比较。跨所有 CPU/libm 永久保证逐位一致，需要进一步内置软件数学函数或预计算超越函数结果；不能把单台机器上的通过外推成所有架构的保证。
