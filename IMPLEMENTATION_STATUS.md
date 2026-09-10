# fast-gatk 当前实现状态

更新时间：2026-09-10

## 目标与不可变约束

目标是实现与 **GATK 4.6.2.0** 行为对齐的核心调用能力，采用不可改变的分层：

```text
C++ Host
  ├─ BAM/VCF/FASTA I/O、ActivityRegion 调度、CIGAR、EventMap、变长 haplotype
  ├─ GATK 规则的候选集、等位基因与 gVCF 输出结构决策
  └─ 向 kernel 提供扁平、确定的数值输入

Kokkos kernel
  ├─ activity profile、图计数、PairHMM、Smith-Waterman 评分
  ├─ genotype GT/GQ、AF/QUAL、ReferenceConfidence 数值归约
  └─ Serial / OpenMP（以及可用的其他 execution space）的一致数值实现
```

该边界没有被改成纯 CPU 实现，也没有把 CIGAR、HTSlib 或不规则图状态塞进 kernel。

“1:1 完成”目前**尚未得到全范围证明**；文档中的完成度仅用于项目管理，不可替代逐项 oracle。

## 完成度快照

| 范围 | 当前估计 | 说明 |
| --- | ---: | --- |
| 全仓库工具面（含长尾） | 约 78% | 长尾小工具按项目目标暂不优先。 |
| 目标核心：HC / Mutect2 / BQSR | 约 93% | 核心路径已有大量 GATK oracle；仍有完整回归和边界覆盖未完成。 |
| HaplotypeCaller 多等位 / gVCF 收尾项 | 进行中 | 本轮解决了两个已复现的源码语义差异，仍需扩展回归。 |

## 本轮已完成并有当前证据的工作

### HaplotypeCaller 多等位基因上限

- `--max-genotype-count` 已按 `HaplotypeCallerGenotypingEngine.removeAltAllelesIfTooManyGenotypes()` 的含义处理：Host 使用 K-best haplotype score 的最高/次高分选择 ALT，发生在 PairHMM 边缘化之前。
- `--max-alternate-alleles` 已与 GATK 的 GL 驱动排序及 Number=G PL 重映射对齐，不能用前者替代后者。
- 四倍体、多 ALT、`max-alternate-alleles=1`、`max-genotype-count=5` 的 VCF 和 gVCF 候选行均有 GATK 精确 oracle。

主要代码：

- `fastgatk-native/src/calling_pipeline.cpp`
- `fastgatk-kernels/src/kmer_graph.cpp`
- `fastgatk-kernels/include/fastgatk/kernels/kmer_graph.hpp`
- `fastgatk-native/scripts/verify_hc_complex_multiallelic_oracle.py`

### gVCF spanning deletion 与符号等位基因

- gVCF 在上游缺失覆盖当前位置时，按 `REF, ALT…, *, <NON_REF>` 的原始矩阵处理 PL。
- 当 `--max-alternate-alleles` 截断掉 `*` 时，不会在输出层因“被上游缺失覆盖”而把它错误加回。
- 实测 GATK 行：`AT,*,<NON_REF>` 在 ALT 上限为 1 时变为 `AT,<NON_REF>`；QUAL、MLEAC/MLEAF、GT、AD、PL、PGT/PID/PS、SB 均精确一致。
- 流式 gVCF 输出路径已将实际 read batch 传给注释重算，而非使用空 read 指针。

当前 OpenMP oracle：

```bash
python3 fastgatk-native/scripts/verify_hc_chr20_max_mnp_gvcf_gatk_oracle.py \
  --native fastgatk-native/build/fastgatk-hc-call
```

结果：通过（含原始跨度删除、max-ALT、普通 VCF 抑制伪 `*`）。

### `--alleles` / GenotypeGivenAlleles 的注入边界

本轮发现并修复的差异不是“零 read 支持 ALT 过滤”，而是更具体的 EventMap 时序差异：

1. GATK 的 `AssemblyResultSet.addGivenAlleles()` 选择参考 haplotype 加最高分的最多 5 个 base haplotype。
2. 仅向 EventMap 不重叠的 base haplotype 注入 feature event。
3. 重新生成 EventMap，再进行 PairHMM 与 genotyping。

此前 Native 在 EventMap 已生成后直接补候选，遇到 feature ALT 与既有 MNP 重叠时会生成 GATK 不存在的 ALT 列。

现在 Native 在 C++ Host 的同一 AssemblyResult 阶段注入人工 haplotype；然后使用既有 Kokkos PairHMM/AF/GT-GQ 数值路径。若该 event 已存在，或所有受限 base haplotype 都与其重叠，则不再事后制造候选。

新增的 chr20 oracle 覆盖：feature `CA -> AT,CG` 与已有 MNP/EventMap 重叠时，GATK 和 Native 都只保留 `AT,*,<NON_REF>`，不会输出伪造的 `CG`。

当前 OpenMP 验证：

```bash
python3 fastgatk-native/scripts/verify_hc_alleles_gatk_oracle.py
python3 fastgatk-native/scripts/verify_hc_chr20_max_mnp_gvcf_gatk_oracle.py \
  --native fastgatk-native/build/fastgatk-hc-call
```

两项均已通过。前者覆盖多 ALT、锚定 indel、过滤 feature、两个 filtered-feature 兼容开关、空 pileup 与 gVCF。

## 最近验证状态（当前二进制的实测证据）

证据目录 `.diag/regression/`（含各次全量日志与自动生成的证据块）。

| 范围 | OpenMP | Serial | 备注 |
| --- | --- | --- | --- |
| 全量回归（重建后复跑） | **280/280 通过**（975s） | **280/280 通过**（968s） | 二进制重建后的当前证据，见下方「关于二进制的更正」 |
| HC/Mutect2 子集 | **72/72 通过**（261s） | **72/72 通过**（268s） | `omp-hc-mutect2.log` / `serial-hc-mutect2.log` |
| `verify_hc_alleles_gatk_oracle.py` | 通过（86s） | 通过（112s） | 原「待本轮复跑」项已结清 |
| `verify_hc_complex_multiallelic_oracle.py` | 通过（45s） | 通过（33s） | 四倍体 / max-ALT / max-genotype-count |
| `verify_hc_chr20_max_mnp_gvcf_gatk_oracle.py` | 通过（47s） | 通过（53s） | MNP + gVCF `*` + max-ALT |

### 关于二进制的更正（重要）

首轮回归确实在 `fastgatk-hc-call`（17:03 / 17:08，比全部源文件新）上取得 72/72 与 280/280。
但随后发现 **`fastgatk-mutect2` 是陈旧的**：`fastgatk-native/build/fastgatk-mutect2` 时间戳为
15:00，而共享的 `fastgatk-native/src/calling_pipeline.cpp` 是 17:03 —— 也就是说，首轮
Mutect2 的测试结果**并没有覆盖本轮共享 calling pipeline 改动**。

已据此重建两个后端（`fastgatk-mutect2` 重新链接为 18:17，证明该改动确实影响 Mutect2，
因为共享库变化会触发重链），并**重跑全量**：

- OpenMP 280/280（975s）、Serial 280/280（968s），证据目录
  `.diag/regression/20260910-181746/`，日志 `omp.log` / `serial.log`。

因此当前可用的结论是：**共享 calling pipeline 改动未使 Mutect2 回归**（在该测试集覆盖范围内）。
这也是「必须核对二进制新鲜度、不能只看测试是否全绿」的一个实例——
`run_regression.sh` 现已内置源码比二进制新的告警。

> 注意：本节的绿色只对本节记录的二进制与 git 版本有效。任何生产代码改动都必须重新取得
> 证据，不得沿用本节数字。

已核销的历史条目：`fastgatk-mutect2-gvcf-reference-blocks-gatk-contract` 曾被
`LastTestsFailed.log`（11:58）记为失败，实为陈旧记录；当前单跑 7.65s 通过。

> 注意：本节的绿色只对本节记录的二进制与 git 版本有效。任何生产代码改动都必须重新取得
> 证据，不得沿用本节数字。

## 未完成事项

### 近期必须完成

1. ~~完成最新改动后的 Serial 编译与点名 oracle~~
   **已完成（2026-09-10）**：`fastgatk-native/build-serial/fastgatk-hc-call` 已比全部源文件新；
   三个点名 oracle 在 OpenMP 与 Serial 双后端全部通过（见「最近验证状态」）。

2. ~~运行最新二进制的完整 HC/Mutect2 回归~~
   **已完成（2026-09-10）**：双后端 HC/Mutect2 子集 72/72；并额外取得全量 280/280 的双后端通过证据。
   原文要求的命令现由 `run_regression.sh` 统一封装：

   ```bash
   fastgatk-native/scripts/run_regression.sh -R 'fastgatk-(hc|mutect2)'
   fastgatk-native/scripts/run_regression.sh          # 双后端全量
   ```

3. 对 gVCF 流式/分区路径扩大 oracle 覆盖，特别是多个重叠 indel、符号等位基因、注释
   read-overlap margin 与重分块的组合。**进行中**，见「本轮并行工作」。

### 仍需审计的核心一致性边界

- `--alleles` 的更多复杂注入情形：重复序列中的 indel、跨 AssemblyRegion 边界、多个相邻或重叠 feature event、base-haplotype 排名并列。
- 复杂 EventMap 与图路径的完整性：Native 当前以已物化的 K-best 路径和 Host CIGAR 作为受限注入母体；需要用更多 GATK fixture 证明其在高复杂区域等价。
- Mutect2：需要在最新 HC 共享 calling pipeline 改动后重跑完整 oracle，确认正常样本重放、联合 AssemblyResultSet 和 somatic gVCF 不受影响。
- BQSR 与其余核心路径的扩大真实数据/极端参数验证；当前已有实现和覆盖，但尚不足以声称全参数空间 1:1。

### 暂不优先的范围

- 长尾小工具和不影响 HC / Mutect2 / BQSR 加速库核心价值的工具面。
- 这些工具不是“已完成”；只是依据项目目标被明确排在核心一致性之后。

## 本轮并行工作（2026-09-10）

三路审计并行进行，全部为「只读审计 + 新增 oracle」模式，**不改生产代码**，
产出统一为差异清单（最小复现 + GATK/Native 差异行 + 首个差异字段 + 源码级假设）：

| 轨道 | 范围 | 产出 |
| --- | --- | --- |
| Track A | gVCF 流式/分区：重叠 indel、符号等位基因 × max-ALT、read-overlap margin、重分块组合 | **发现 P1 偏离** → `fastgatk-native/evidence/gvcf-stream-by-region-divergence-20260910.md` |
| Track B | `--alleles` / GenotypeGivenAlleles 复杂注入：重复序列 indel、跨 AssemblyRegion、相邻/重叠 feature、排名并列 | `.diag/track-b-alleles-findings.md` |
| Track C | Mutect2 独立复核：共享 calling pipeline 改动后是否回归、正常样本重放、联合 AssemblyResultSet、somatic gVCF、TLOD 精度 | `.diag/track-c-mutect2-recheck-findings.md` |

### 执行状态与阻塞

三路审计以子代理并行启动，但**全部失败于子代理提供方的配额上限**
（`minimax-cn` 返回 `429 rate_limit_error: 已达到 Token Plan 用量上限`，见各子代理
session 记录；三次尝试均在 turn 1 即失败）。因此第 1 轮改为由主会话接管单线执行，
并行度受此限制。子代理已完成的部分产物（三个脚本）被保留并复用：

- `verify_gvcf_stream_overlapping_indels_gatk_oracle.py`（Track A）→ 已跑通，
  并据此**独立复现**出 P1 偏离（见下）。
- `verify_mutect2_recheck_normal_replay.py`、`verify_mutect2_recheck_assembly_resultset_joint.py`
  （Track C）→ 尚未执行。

| 轨道 | 状态 |
| --- | --- |
| Track A | **已产出 P1 结论**（`--stream-by-region` 偏离 GATK 与非流式路径），证据已归档 |
| Track B（`--alleles` 复杂注入） | **未开始**（子代理无产物，配额阻塞） |
| Track C（Mutect2 独立复核） | **部分完成**：已发现 mutect2 二进制陈旧并重建复跑（280/280）；两个复核脚本未执行 |

### Track A 结论：`--stream-by-region` 偏离（P1）

在 `fixtures/chr20/mnp.bam` + `20:10019901-10020710` 上，`fastgatk-hc-call
--stream-by-region N`（N < 区间长度）与**同命令去掉该选项**的输出不一致，也与
pinned GATK 4.6.2.0 不一致：

- **参考块粒度改变**：非流式与 GATK 均为 48 条记录（把 10020230–10020428 拆成 29 个
  `<NON_REF>` 块）；tile=500 时流式只有 20 条，28 条被合并成单块。
- **注解读数错误**：20:10020680 行 `RAW_MQandDP` 由 `28800,8` 变为 `97200,27`，
  `SB` 由 `0,0,3,3` 变为 `0,0,0,0`；而 GATK 与非流式 native 在该行**逐字段相同**。
  `QUAL/PL/GT/GQ/PGT/PID/PS` 不受影响。
- tile 扫描显示**只有 tile 覆盖整个区间时才逐字节相同**；OpenMP 与 Serial 偏离一致，
  指向 Host 路径而非 kernel。
- 现有测试 `fastgatk-hc-region-streaming-contract` 断言 tile=500 的 key 集合与非流式相同，
  但该断言在本次夹具上不成立 —— 说明该契约当前**依赖夹具**，需补回归守卫。

完整证据、最小复现命令与根因分析见
`fastgatk-native/evidence/gvcf-stream-by-region-divergence-20260910.md`。

### 第 2 轮定位进展（修复**尚未开始**）

已把基准钉死：**非流式 native 与 GATK 在该区间上全部 48 条数据行逐字段相同**
（含 29 块碎片结构），因此非流式是对的、流式是偏离方。

三类偏离（第 3 类为本轮新发现）：

1. **参考块粒度改变**：非流式/GATK 把 10020230–10020428 拆成 29 个小 `<NON_REF>` 块，
   流式合并成 1 块（tile=100/500 时少 28 条记录）。
2. **注解读数错误**：10020680 处 `RAW_MQandDP` `28800,8 → 97200,27`、`SB`
   `0,0,3,3 → 0,0,0,0`。
3. **凭空产生 phasing（新）**：10020228/10020229 处，非流式/GATK 输出未定相 `GT=0/1`，
   流式输出已定相 `0|1`/`1|0` 并新增 `PGT/PID/PS` —— 属输出语义改变，不只是数值偏差。

**只与非流式等价的条件**：tile × `--assembly-region-padding` 扫描（差异行数，越低越好）

| tile | pad=100 | pad=200 | pad=300 |
| ---: | ---: | ---: | ---: |
| 300 | 1 | 1 | 1 |
| 405 | 43 | 43 | 3 |
| 500 | 43（仅 20 条记录） | 3 | 3 |
| 810（=区间长度，单 tile） | **0（逐字节相同）** | 1 | 1 |

结论：**默认 pad=100 下只有单 tile 才逐字节相同；加大 padding 不是修复**
（它修回 tile=500 的块结构，却弄坏原本逐字节相同的单 tile 情形）。

**一处自我更正**：先前记录的「拼接器 `can_merge` 已被证伪」作废 —— 该推理假设各 tile 产出的块
与非流式相同，前提不成立（各 tile 产出的块本身就是退化的 `PL=0,0,0`）。从坐标看拼接器必然参与
（记录跨 tile=500 的 core 边界 10020401），但它只是"把两个已退化的块接起来"的机制。

**仍未定位到具体代码**：候选为 `build_reference_blocks`（精确）vs
`build_profile_local_reference_blocks`（近似，`calling_pipeline.cpp:13111/15206/15672`），
以及 `run_region_streaming` 的 `input_halo_intervals`（`hc_call.cpp:6463`）的窗口口径。
`FASTGATK_DEBUG_GVCF_EMISSION=1` 显示 `symbolic-pl=unavailable` 在**两种模式下都出现**，
不是判别依据（先前把它当作线索，一并更正）。D1（phasing）与 D2（注释作用域）是否同源未判定。

**这同时暴露了一个方法学教训**：状态文档曾把「流式路径已传实际 read batch（非空指针）」
当作已修正的证据，但"非空"不等于"范围正确"。修复必须由 oracle 判定，不能由改动描述判定。

### 并行原则（本仓库硬约束）

**oracle/夹具/审计可并行；生产代码改动必须单线串行**，因为核心路径的全部修复都落在
同一批文件（`calling_pipeline.cpp`、`kmer_graph.cpp`、`hc_call.cpp`、`mutect2_tool.cpp`），
并行编辑必然互相覆盖。所有 oracle 都使用 `tempfile.TemporaryDirectory`，因此测试进程
可安全并发。

## 工程基线

本轮补齐的四项工程前提（此前缺失）：

1. **版本控制**：本仓库此前不是 git 仓库，导致只能写「本轮改动前通过」这类不可验证基线。
   现已建立基线提交 `180c2bf`（569 文件 / 16MB），跟踪源码、脚本、文档与小型夹具；
   忽略 `third_party/`、`gatk-source/`、`gatk-rs-source/`、根目录 `testdata/`、`fixtures/`、
   各构建树、`work/`、`.ab4/`、`.diag/` 等大体积内容与产物。规则见 `.gitignore`。
   注意：`fastgatk-native/tests/fixtures/` 下的小型夹具**必须入库**（被
   `verify_overlapping_quality_correction.py` 等引用），因此忽略规则只锚定仓库根目录。

2. **回归入口**：`fastgatk-native/scripts/run_regression.sh` 一条命令跑双后端，
   存档日志到 `.diag/regression/<时间戳>/` 并生成可直接粘贴的证据块。
   它同时做陈旧性检查（源码比二进制新则告警）、工具链前置检查，并在
   **过滤器未命中任何测试时判为「无证据」而非通过**（ctest 本身此时返回 0，会伪造绿色）。

3. **权威构建目录**：只有两个后端是权威的 ——
   `fastgatk-native/build`（OpenMP）与 `fastgatk-native/build-serial`（Serial）。
   仓库根目录存在一份 **孤儿 in-source 构建配置**（`CMakeCache.txt`、`Makefile`、
   `CMakeFiles/`、`CTestTestfile.cmake`，9 月 4 日配置），它只注册 30 个 hc/mutect2 测试、
   且没有任何二进制产物，并已把生成物散落到 `fastgatk-core/`、`fastgatk-kernels/`、
   `fastgatk-runtime/`。**不要用它**；引用回归结果时务必写明是哪一份构建目录。

4. **工具链解析**：`java` 不在 PATH 上；oracle 依赖 vendored 的
   `third_party/jdk17/bin/java` 与
   `third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar`。
   历史脚本的解析方式不一致（部分认 `JAVA` 环境变量、部分硬编码）。现已提供
   `fastgatk-native/scripts/oracle_toolchain.py`：优先级为「显式参数 > `JAVA`/`GATK_JAR`
   环境变量 > vendored 默认路径」，缺失时给出明确错误。**新写的 oracle 请使用它**，
   存量脚本的批量迁移属后续工作。

### 待处理的口径问题

- `progress_score.json` 的 `global_score` 为 0.757、`hc_mutect2_bqsr` 为 0.915，但
  `tool_audit_expected_entries` 声明 48 条而 `tool_audits` 实际只有 2 条（Mutect2 0.833、
  FilterMutectCalls 0.912）。本文件的 78% / 93% 与该分数目前**没有同源证据链**，
  在补齐 tool audit 之前，两个数字都只能当作项目管理估计。
- **静默 skip 风险（已实测；当前环境未发生）**：117 个 `verify_*.py` 在缺少 GATK/oracle
  输入时会打印 `{"status":"skip"}` 并返回 0；注册测试中只有 134 个设置了
  `FASTGATK_REQUIRE_GATK_ORACLE=1`，因此共有 **32 个「可静默跳过且未设强制」的候选**。
  已对这 32 个候选做定向 verbose 复跑（`.diag/regression/skip-audit-verbose.log`）：
  32/32 通过，全部产出了 `status` 字段，**skip 计数为 0**（日志中唯一的 "skip" 字样来自
  `fastgatk-mutect2-tlod-formula` 的 `multiallelic_loci_skipped: 0` 字段）。
  结论：本机当前的「全绿」是真实执行结果，不是静默跳过。但这 32 个测试在缺少 bundled
  GATK/testdata 的机器上会**静默变绿**，属可移植性/CI 隐患，应逐一补上
  `FASTGATK_REQUIRE_GATK_ORACLE=1`。

  说明：`ctest -O`（output log）**不会**记录通过用例的 stdout，只有 `-V` 才会；因此
  「统计静默跳过」必须用 `-V`。这是上一轮误判「0 skip」的原因，已修正。

## 下一步执行顺序

1. ~~结束并确认 Serial 构建，完成 Serial oracle。~~ 已完成，见「最近验证状态」。
2. ~~跑完整 HC/Mutect2 回归。~~ 已完成；并已取得全量双后端 280/280 证据。
3. 第 1 轮差异清单已完成第一条（Track A 的 `--stream-by-region` P1）；
   Track B / Track C 的剩余部分见「执行状态与阻塞」。
4. **单线串行**落地修复，优先修 Track A 的 P1（`--stream-by-region` 的注释 read batch
   与参考块边界）：一次一个改动 → 双后端增量构建 → 定向 oracle + 72 子集 → 全量回归 →
   commit（有 git 后「改动前/后」才可验证）。修复后把
   `verify_gvcf_stream_overlapping_indels_gatk_oracle.py` 的断言由 diagnostic 改为
   「流式 == 非流式 == GATK」并注册进 CTest。
5. 补齐 Track B（`--alleles` 复杂注入）与 Track C 的两个复核脚本；子代理配额恢复前由主会话
   执行，或等配额恢复后重新并行下发。
6. 补齐静默 skip 的 32 个候选测试（补 `FASTGATK_REQUIRE_GATK_ORACLE=1`）与 tool audit 口径，
   使百分比与证据同源。
7. 只有当每个核心选项都有匹配范围的当前证据时，才讨论「1:1 完成」。

