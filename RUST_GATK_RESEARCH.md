# Rust-GATK（gatk-rs）源码调研

调研日期：2026-08-20  
源码目录：`/home/turing-agents/Documents/fast-gatk/gatk-rs-source`  
上游仓库：[IPNP-BIPN/gatk-rs](https://github.com/IPNP-BIPN/gatk-rs)

## 1. 结论先行

`gatk-rs` 很有价值，但当前不能直接替换生产环境中的 `gatk` 命令。它更准确的定位是：

1. 按 GATK 4.6.2.0 源码逐符号翻译的 Rust 语义/兼容性基线；
2. 用 Java oracle 做差分测试、收集字节级 golden 的验证工程；
3. 可供我们的 Kokkos 项目复用的 host 层数据模型、错误语义和 conformance 资产。

它不是 Kokkos 加速实现，也不是已经完成的 GATK native runtime：仓库没有 `main()`、统一 CLI dispatcher 或可执行 `gatk` 二进制；核心变异调用器和大部分集群型工具仍未开始。最适合的路线是“Rust 兼容层 + Kokkos 热点内核 + 兼容性 launcher”，而不是把 `gatk-rs` 当成最终 GPU/集群引擎。

## 2. 版本、来源和源码盘点

本地下载方式：

```text
git clone --depth 1 https://github.com/IPNP-BIPN/gatk-rs.git gatk-rs-source
```

当前源码：

| 项目 | 值 |
|---|---|
| 分支 | `main` |
| HEAD | `1c5b232b09931dcc1029a4206cb1855eb1cabe5c` |
| HEAD 提交 | 2026-08-20，合并 PR #459 |
| 参考 GATK | 4.6.2.0，`76edc75c26504da94bbaee66584e107e76ee15de` |
| 参考依赖 | Picard 3.4.0、htsjdk 4.2.0、Barclay 5.0.0、Intel GKL 0.8.11 |
| Rust | `rust-toolchain.toml` 固定 1.97.1，edition 2021 |
| 仓库状态 | WIP，明确声明不是 Broad 官方 GATK |
| 本地规模 | 6 个 crate、374 个 `.rs` 文件、185 个集成测试文件 |

本地分析环境没有安装 `cargo/rustc`，所以本次只完成源码、Cargo manifest、CI 配置和静态盘点；没有把“本地编译通过”写成结论。上游 CI 的测试步骤仍以仓库当前 workflow 为准。

README 的工具盘点存在两个口径：历史方法部分写 375 个工具/10,796 个参数，当前生成盘点写 311 个工具/13,130 个参数。以当前 `docs/STATUS.md` 为准：311 个工具，其中 89 个 oracle-backed、12 个 unchecked、210 个 not started；GATK-origin 工具为 202 个，当前约 58 个已经开始。

## 3. Crate 和模块结构

| crate | 当前职责 | 对我们方案的价值 |
|---|---|---|
| `gatk-engine` | interval、BAM/FASTA 查询、GATKRead 适配、CIGAR/剪切、pileup、assembly region、BQSR、Mutect 过滤、数值和 Java 兼容工具 | Rust host 语义层、错误/边界条件、分片模型 |
| `gatk-readfilter` | 56 个 read filter 和 JEXL filter | 最适合第一个 Kokkos 整数/字节 kernel |
| `gatk-annotation` | 变异 annotation，README 称 53/54 已测 | VCF/annotation 兼容性和 golden |
| `gatk-tools` | walker 与工具的 library 实现，当前约 65 个源模块 | 记录变换、BQSR、pileup、VCF 工具的参考实现 |
| `gatk-barclay` | Barclay/jopt-simple 参数值模型、arguments file、argument collection | 兼容 GATK 参数和拒绝语义 |
| `gatk-corpus` | 测试 corpus 解析支持 | 差分测试输入层 |

依赖策略很谨慎：htsjdk-rs 组件按 git revision 固定；`noodles-fasta` 和 `noodles-bam` 只承担 `.fai`/`.bai` 等格式 plumbing，GATK 对记录选择、参考序列大写/IUPAC 展开、interval 合并等语义仍在本仓库中实现。这种“格式依赖、语义自持”的边界值得保留。

## 4. 已有实现和未覆盖范围

### 已有、且有较强验证的部分

- 全部 56 个 read filter：README 声称 4,661 个决定与参考一致。
- `PrintReads`：6 个 BAM 和 `.bai` golden 在声明的 JDK deflater 条件下字节一致。
- `ReadWalker`、`IntervalWalker`、`LocusWalker`、`AssemblyRegionWalker` 的大量 traversal 语义。
- indexed FASTA、BAM/.bai 查询、interval `-L/-XL`、CIGAR、read clipping、pileup、downsampling。
- BQSR 的 report/table、`BaseRecalibrator`、`ApplyBQSR` 和相关数值逻辑。
- Mutect2 的过滤/统计子逻辑、CNV 的部分工具、VCF/annotation 工具。

### 明确未完成、不能视为替换能力的部分

根据当前 `docs/STATUS.md`/`ROADMAP.md`，以下关键工具仍在 not started 或 G3/GPU 里：

- `HaplotypeCaller`、`Mutect2` 主 caller、assembly graph、genotyping；
- PairHMM（ROADMAP 明确仍是待办）；
- `GenotypeGVCFs`、`CombineGVCFs`、`GenomicsDBImport`；
- `MarkDuplicates`/`MarkDuplicatesSpark` 及完整 Picard 依赖链；
- Spark 39 个工具、CRAM、完整 ML inference、SV/PathSeq 等；
- 统一命令行、usage text 和 plugin discovery。

因此，它还不能运行常见的 GATK germline/somatic 主流程，也不能作为现有 Nextflow DSL 中 `gatk HaplotypeCaller` 的透明替代。

## 5. 并行、硬件和 I/O 现状

### 硬件

源码没有 Kokkos、CUDA、HIP、SYCL、OpenMP 或 GPU kernel；没有 `main()`，也没有 `clap`/其他 CLI 框架。源码中的少量 `Spark`/`rayon` 命中主要是上游语义注释或间接依赖，不构成当前的分布式执行实现。

仓库自己的 GPU milestone 也明确写出：GATK 4.6.2.0 本身没有 GPU kernel，PairHMM 是 GKL/纯 Java 路径；GPU 只能作为第二实现，并且必须通过同一 byte-identity gate。这个判断与我们的设计一致：GPU 只能是可验证的加速后端，不能成为兼容性唯一实现。

### I/O 和内存

当前代码是 library 级、以 conformance fixture 为中心的实现：

- FASTA 使用 `noodles-fasta` indexed reader，语义层再做 upper-case 和 IUPAC→N；
- BAM 解码/BGZF 使用 htsjdk-rs 组件，`.bai` plumbing 使用 `noodles-bam`；
- 没有面向 WGS 的 bounded pipeline、异步预取、NUMA placement、GPU staging、spill manifest 或外部排序框架；
- 临时目录主要出现在测试 fixture；仓库没有 MarkDuplicates/GenomicsDB 级别的生产内存压力实现。

ROADMAP 的速度 milestone 也写明“目前没有任何性能数字”，并要求先建立同输入、同 conformance fixture 的 wall time/CPU time/峰值 RSS 基线。因此，不能仅凭“Rust”宣称优于 Java GATK。

## 6. 验证方法是它最值得借鉴的部分

`gatk-rs` 的强项不是完成度，而是验证纪律：

- 固定 GATK 参考 commit、JDK 17、真实 x86-64 oracle container；
- t-wise covering arrays（`t=2`，关键路径 `t=3`）；
- coverage-guided differential fuzzing；
- 每个 golden 记录 raw bytes 与 canonicalization 规则；
- 对 Java `Math.log/exp/pow`、deflater、Double formatting、集合遍历顺序等可观察差异做显式隔离；
- CI 运行 `cargo fmt`、`clippy -D warnings`、release workspace test、外地 locale/TZ/TMPDIR/单线程 determinism，以及两次干净构建的可复现性检查。

我们的项目应直接吸收这套“oracle + golden + determinism + provenance”框架，而不是只吸收 Rust 代码。

## 7. 对 GATK+SLURM+Nextflow 直接替换的判断

### 现在不能直接替换

缺口不在 SLURM 或 Nextflow，而在运行时合同：

1. 没有 `gatk` launcher 和命令注册表；
2. 没有统一参数解析、`@arguments_file`、`--java-options` 兼容入口；
3. 关键 caller、joint genotyping、duplicate marking、GenomicsDB 未完成；
4. 没有声明稳定的 exit code、stderr、header、BGZF/VCF 输出合同；
5. 没有生产级资源/临时目录/并发控制。

### 可实现的兼容外壳

建议我们做一个兼容 launcher，先不改变 Nextflow/SLURM 的调用方式：

```text
gatk <tool> [原有参数...]
  └─ 解析 GATK 参数和 --java-options
  └─ FAST_GATK_BACKEND=java|rust|kokkos|auto
  └─ 按工具能力选择 native 实现；未覆盖时回退 Java GATK
  └─ 保持 stdout/stderr、退出码、header、索引和失败信息合同
```

Nextflow 仍然负责样本级/interval 级 DAG，SLURM 仍然负责 job、CPU、内存、GPU、scratch 和容错；native runtime 内部只做线程/队列，不嵌套 `sbatch`。这样可以在一个流程中逐工具切换，不需要重写 DSL。

## 8. 与 Kokkos 方案的推荐组合

### 分层

```text
兼容 gatk launcher（参数/退出码/日志/回退）
        |
Rust host compatibility layer
(BAM/VCF/FASTA、interval、walker、资源预算、spill、manifest)
        |
Kokkos kernel layer
(read filter/pileup → BGZF → PairHMM → assembly/genotyping)
        |
CPU SIMD / OpenMP / CUDA / HIP / SYCL
```

### Rust 负责什么

- CLI/arguments-file/配置、工具注册和 feature flags；
- BAM/VCF/FASTA record 生命周期、错误和资源 RAII；
- interval/shard/assembly-region 的确定性调度；
- memory budget、bounded queue、scratch/spill manifest、checkpoint；
- SLURM cgroup/环境读取，以及 Nextflow job 的单进程资源合同；
- oracle/golden/determinism harness。

### Kokkos 负责什么

- 结构化、批量、热点计算：read filters、pileup、BGZF block、PairHMM、SW、局部 assembly 和后续 genotyping；
- CPU SIMD、OpenMP、CUDA、HIP、SYCL 的统一 execution/memory space；
- team/vector parallelism、scratch memory 和异步 host/device copy；
- 所有浮点 kernel 固定遍历顺序、禁止 fast-math/FMA 重排；不能满足 byte identity 时必须回退 CPU 并标成 bio-identical。

### Rust↔Kokkos 接口原则

- 用稳定 C ABI 或极薄 `cxx`/bindgen shim，不把每条 read 逐条跨 FFI；
- 以 shard/assembly-region 为批次，传递结构化 contiguous buffers；
- 返回 opaque execution/event handle，显式 `wait/fence`；
- host 侧保留 canonical CPU reference，GPU 结果先过同一 golden comparator；
- Kokkos build 选择 backend 时，镜像/模块名和 `FAST_GATK_BACKEND` 一起写入 provenance。

## 9. 建议的落地顺序

| 阶段 | 目标 | 结果 |
|---|---|---|
| R0 | 引入 `gatk-rs-source`、GATK 4.6.2.0 reference、oracle/golden harness | 先锁兼容性事实，不宣称性能 |
| R1 | 兼容 launcher + `PrintReads`/read filters/interval 查询 | 可在 Nextflow 中对少数工具透明替换 |
| R2 | Rust bounded I/O、scratch/spill、SLURM resource adapter | 形成生产级单节点 native runtime |
| R3 | CPU SIMD/Kokkos read filter、pileup、BGZF | 先做低风险整数/字节热点 |
| R4 | PairHMM CPU reference，再做 Kokkos CPU/GPU | 以 LoglessPairHMM byte golden 为门槛 |
| R5 | HaplotypeCaller/Mutect2 assembly、genotyping | 进入真正的 caller 替换 |
| R6 | GenotypeGVCFs/GenomicsDB/MarkDuplicates 和 scatter/gather | 覆盖常见 germline/somatic pipeline |
| R7 | 双架构 GPU determinism、性能基线、逐工具 rollout | 逐步扩大 `auto` native 覆盖率 |

## 10. 风险和决策

1. **完成度风险：高。** `gatk-rs` 的字节级声明是已覆盖 slice 的声明，不是全 GATK 声明。
2. **性能风险：高。** 仓库明确没有性能基线；Java GATK/GKL 的 steady-state 可能已很强，必须实测。
3. **数值风险：高。** Rust、Kokkos、CPU/GPU libm 和 FMA/reduction 顺序可能改变 VCF 数值或 BGZF 字节。
4. **兼容性风险：中高。** 不仅要复制算法，还要复制参数拒绝、header、索引、日志和失败顺序。
5. **许可证风险：可控但必须审计。** `gatk-rs` 使用 Apache-2.0，LICENSE 明确包含 GATK 派生部分；复制代码前仍需保留归属、来源和第三方依赖审计。
6. **仓库同步风险：中。** 目前是滚动 `main`，本地应固定 commit，并将参考 GATK/htsjdk/Picard revision 写进构建 provenance。

最终建议：把 `gatk-rs` 纳入我们的“兼容性参考和测试资产”，不要把它当作已经可生产替换的运行时；Kokkos 项目从 Rust host + C ABI Kokkos kernels 开始，先实现 dispatcher、read filter、pileup、BGZF 和 PairHMM，再逐步覆盖 caller 和 joint-genotyping。

## 11. 只吸收有价值的部分

### 直接吸收

| 吸收项 | 进入我们的项目位置 | 原因 |
|---|---|---|
| 固定 GATK reference commit、依赖 revision 和 toolchain | `provenance.json`、容器标签、构建 manifest | 防止 Java/GATK、Rust、Kokkos 版本漂移导致结果不可追溯 |
| Java oracle + golden corpus | `conformance/` 和 CI | 每个 native 工具先证明 CLI、输出、错误和字节合同，再谈性能 |
| t-wise covering array + 差分 fuzzing | 参数测试和回归测试 | GATK 参数组合不可能穷举，覆盖数组比手工 happy path 更可靠 |
| determinism gate | CI 的 CPU/GPU 双路径 | 固定排序、归约、随机数、压缩和 locale/TZ，避免集群结果漂移 |
| “格式依赖、语义自持”边界 | HTS I/O 层 | 可以使用成熟 BAM/FASTA plumbing，但 interval、filter、header、错误语义必须由我们控制 |
| Java 兼容的错误/边界语义 | Rust host API | GATK 兼容不只是正确结果，还包括拒绝条件、退出码和错误文本 |
| walker/interval/assembly-region 的分层模型 | runtime scheduler | 对应 Nextflow scatter 和 SLURM allocation，便于按 contig/shard 调度 |
| byte-identical 与 bio-identical 分级 | tool capability matrix | 无法保持位级一致时明确降级，禁止模糊地宣称“兼容” |
| license/provenance audit | CI guard | 吸收实现思路时保留 Apache/GATK 归属和第三方依赖审计 |

### 选择性吸收

- `gatk-engine` 的 interval、pileup、read adapter、BQSR 数据结构：只把行为和测试向量纳入我们的 Rust host；不逐行复制成另一套 runtime。
- `gatk-readfilter`：优先作为 Kokkos 第一个 kernel 的语义参考和 golden 来源；GPU 只处理批量 predicate，过滤顺序和保留顺序由 host 固定。
- `PrintReads` 的 header/`@PG`/索引字节测试：作为最小 drop-in 工具和 launcher 验证样板。
- `gatk-barclay` 的参数模型：吸收参数优先级、`@arguments_file`、collection flattening 和拒绝行为；统一 CLI 仍由我们实现。

### 明确不吸收

- 不引入 `gatk-rs` 当前未完成的 HaplotypeCaller、Mutect2、PairHMM、GenomicsDB、MarkDuplicates 作为生产依赖；这些由我们的 Kokkos/runtime 路线重新实现。
- 不照搬其 library-only 形态、分散的 per-tool crate 和 branch 命名；生产项目必须有一个稳定 `gatk` 入口和工具注册表。
- 不把 `noodles` 的具体实现当成最终 HTS I/O 方案；先以兼容性测试决定 htsjdk-rs、htslib、noodles 或自研 codec。
- 不采用“Rust 天然更快”的假设；所有性能结论必须来自同输入、同输出合同、同资源配额的基准。
- 不把其 x86-64/JDK17 oracle 直接当成 GPU 或 ARM 证明；Kokkos 需要额外的跨设备 determinism 和性能测试。

### 项目决策

`gatk-rs` 在本项目中只承担三种角色：**兼容性行为参考、oracle/golden 测试资产、Rust host 设计参考**。不作为运行时依赖、不作为 GPU 实现、不作为完成度或性能依据。
