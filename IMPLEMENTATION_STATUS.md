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
| 全量回归（第 10 轮复验） | **280/280 通过**（1019s） | **280/280 通过**（1017s） | 对 commit `03b02b7`、二进制由 pristine 源重建后的复验；证据目录 `.diag/regression/20260911-010134/` |
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

首次尝试时三路审计的子代理**全部失败于提供方配额**（`minimax-cn` 返回
`429 rate_limit_error: 已达到 Token Plan 用量上限`；三次尝试均在 turn 1 即失败）。
第 3 轮改用 **`deepseek-official` / `deepseek-v4-flash`** 重新并行下发，
三路全部完成。**结论：本仓库的子代理委派必须显式指定 `deepseek-official`**；
`subagent` 工具不暴露 provider 参数，需经 `workflow` 的 provider/model 覆盖来指定
（注意 provider 名是 `deepseek-official`，写成 `deepseek` 会直接失败）。

| 轨道 | 状态 | 报告 |
| --- | --- | --- |
| Track A（gVCF 流式/分区） | **完成** | `fastgatk-native/evidence/gvcf-stream-by-region-divergence-20260910.md` |
| Track B（`--alleles` 复杂注入） | **完成，发现 5 处分歧** | `fastgatk-native/evidence/2026-09-10-parallel-audit/track-b-alleles-findings.md` |
| Track C（Mutect2 独立复核） | **完成** | `fastgatk-native/evidence/2026-09-10-parallel-audit/track-c-mutect2-recheck-findings.md` |
| Track D（流式根因定位，第 3 轮新增） | **完成** | `fastgatk-native/evidence/2026-09-10-parallel-audit/track-d-streaming-rootcause.md` |

### 第 3 轮关键更正：D2 不是流式缺陷，而是 `-L` 窗口依赖缺陷

Track D 发现并经主会话**复核确认**：不加任何 `--stream-by-region`，
仅把 `-L` 起点改到 `20:10020381`，native 就复现出错误的
`RAW_MQandDP=97200,27` / `SB=0,0,0,0`：

| `-L` 窗口（均为非流式） | GATK | native | 数据行差异 |
| --- | --- | --- | --- |
| `20:10019901-10020710` | `28800,8` / `0,0,3,3` | 同 GATK（48 条） | **0 行** |
| `20:10020381-10020710` | `28800,8` / `0,0,3,3`（13 条） | `97200,27` / `0,0,0,0`（13 条） | **2 行** |

即 **GATK 的该注释与 `-L` 窗口起点无关，native 的有关**。含义有两点：

1. 先前「非流式正确、流式偏离」的框架**只对首个窗口成立**，不能外推（已更正两份文档）。
2. 这是一个**比流式问题更广的 parity 缺陷**：普通非流式 HaplotypeCaller 只要换 `-L` 窗口
   （如 scatter 分区）就会复现。**修复优先级应高于流式路径本身。**

### Track B / Track C 要点

- **Track B：`--alleles` 注入边界不是 1:1**（5 处分歧）。已证最小反例两例：
  ① 重叠 feature event —— GATK 输出 2 行（`chr1 900` 与 `chr1 904`），native 只输出 `900`，
  且用 `--drop-alleles` 对照证明两行确为 `--alleles` 驱动；native 调试轨迹显示两个强制事件
  都注入到 `calling_class=1`，其中一个在 `EVENTMAP_REGION` 与 `EVENTMAP_CALL` 之间丢失。
  ② **退出码分歧** —— 某 `--alleles` 记录经 Event 最小化后 ALT 为空，GATK 以 exit 3 中止
  （`IllegalArgumentException: Null alleles are not supported`），native 却 exit 0 输出空 VCF。
  另有 2 处观察项经 `--drop-alleles` 对照判定**并非** `--alleles` 分歧（gVCF indel `END`
  与 `--max-genotype-count 2`，后者 GATK 自身也会崩）。
- **Track C：未观察到 Mutect2 回归**。`ctest -R 'fastgatk-mutect2'` 双后端 **30/30**
  （144.8s / 142.2s）；`verify_mutect2_recheck_normal_replay.py` 双后端通过。
  另有两处**对既有记录的更正**：① 先前记录的 Δ3.5@69368 **已不可复现**（GATK 与 native
  同为一个值，旧 AS_SB_TABLE 缺陷已消失）；② 引擎 parity 头号数字不是 1.4e-14，现测为
  max **4.69e-13**（双后端一致，仍远低于其 <1e-8 门限）；③ `FASTGATK_TLOD_FULL` 并非环境通道，
  实际门控是 `FASTGATK_DEBUG_TLOD=1`。
  交回的 `verify_mutect2_recheck_assembly_resultset_joint.py` 原稿有缺陷（matched normal
  未接入且断言本身对 GATK 也不成立），已由该轨道改写。

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

### 第 2 轮定位进展（**D2 已修复并落地**；D1/D3 经判别实验证明非缺陷）

**最新状态（第 5 轮起）：D2 —— 唯一被证明「同一 `-L` 窗口下 native ≠ GATK」的真 bug —— 已修复。**
改动仅在 `fastgatk-native/src/hc_call.cpp`（+53/-2）：新增 `owner_has_pairhmm_context()`，
并在注释时优先选用**真正携带该等位基因 PairHMM context 的 AssemblyRegion owner**
（twin owner），历史顺序保留为回退。根因是同一等位基因出现在多个 owner 中，
被注释用的那个 owner 从未为其请求 PairHMM，故 context ordinal 是哨兵 `UINT32_MAX`、
likelihood 行为 `-inf`，导致 MQ 证据门落到几何回退谓词（27 条读而非 GATK 的 8 条）
且 strand 证据归零。

门禁（已注册进 CTest，strict 形态，当前绿）：
`fastgatk-hc-window-invariance-gatk-oracle`（12 窗口逐窗口比对 GATK，
外加 `POS 10020680` 的 `RAW_MQandDP`/`SB` 门控，约 50–90 s）。

**分类收敛结论（重要，取代此前的「三类偏离 = 三个待修缺陷」）**：判别实验证明
**D1（凭空 phasing）与 D3（参考块粒度）都不是缺陷** —— pinned GATK 自己在相应窗口下
就产生同样的 phased 行与同样的退化块，native 逐字节相同。它们与
`--stream-by-region` 的关系属**范围/语义问题**：GATK 的 gVCF 输出本身依赖 `-L` 窗口
（块粒度、注释值、是否定相都会变），而流式路径按 tile 窗口分别求值再拼接，
**原理上无法**复现「对用户 `-L` 跑一次」的结果。处置选项：显式定义其语义并记录为
已文档化分歧（GATK 无该开关），或放宽内存上界、放弃分块。
连带影响：`fastgatk-hc-region-streaming-contract` 断言「流式 key 集合 == 非流式」
在当前实现上原理性不可满足，应重新定义或降级为语义契约。

以下为 D2 的定位过程（保留作为证据链）。

三类偏离（第 3 类为本轮新发现）：

1. **参考块粒度改变**：非流式/GATK 把 10020230–10020428 拆成 29 个小 `<NON_REF>` 块，
   流式合并成 1 块（tile=100/500 时少 28 条记录）。
2. **注解读数错误**：10020680 处 `RAW_MQandDP` `28800,8 → 97200,27`、`SB`
   `0,0,3,3 → 0,0,0,0`。已确认该缺陷**是「分块相关」的**：默认 pad=100 下
   tile=200/300/600 取错误值，tile=405/500/700/810 取正确值；`--assembly-region-padding`
   也会翻转结果（tile=500：pad=100 正确、pad=300 错误）。即**流式路径的注释不满足分块不变性**，
   而非流式路径按构造满足。
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

**已定位到具体代码（第 4–8 轮，取代下方旧表述）**：

- D2 的真实性质是 **`-L` 窗口依赖**，**不是流式缺陷**：不加任何 `--stream-by-region`，
  仅把 `-L` 起点改为 `20:10020381` 即复现错误的 `RAW_MQandDP=97200,27` / `SB=0,0,0,0`；
  GATK 在两个窗口下都恒为 `28800,8` / `0,0,3,3`（已用 pinned jar 的
  `--debug-genotyper-output` 独立复核：两窗口均为
  `Event ... alleles=[CA*, *, AT] ... with 8 reads and 20 disqualified`）。
- 机制链（插桩实测）：窗口 B 在 0-based `10020679` 处 `result.candidates` 里有**两份完全相同的
  `CA→AT` 候选** —— 一份已配对（context ordinal 0），一份哨兵（`UINT32_MAX`）。
  渲染采用哨兵那份 → `have_context_mapping_evidence=0` → MQ 证据门落到几何回退分支
  → 放行 27 条读（应 8 条），同时 strand 证据归零。
- 重复候选的**产生点已证实**是 `calling_pipeline.cpp:13738` 的无去重 append（partition merge）：
  窗口 B 的重复位置正是 `10020679`，窗口 A 的重复在其它位置、不含该位点。
- **方案 1（在 merge 处成组去重）已实施并被验收判据否决**：补丁使窗口 A 从 0 行差异变成
  **4 行差异**（回归），且窗口 B **未被修复**。已 `git checkout` 回退并复验基线。
  两条硬信息：① 其它位置的重复候选是**承重**的，删除会改变输出；
  ② 删掉 `10020679` 的重复并未让渲染改用已配对的那份 ⇒ 渲染的选择不由 merge 先后决定，
  需插桩确认两份候选在 merge 前后的相对顺序与各自 ordinal。
- **下一步（唯一推荐方向）**：**方案 2** —— 渲染/注释处优先采用已配对孪生。
  注意第 6 轮已记录的限制：**只换 context ordinal 治不了 `SB`**（strand 读的是该候选自己的
  likelihood 行，重复候选那几行为 `-inf`），故方案 2 必须让 **strand 路径一并采用孪生的
  likelihood 行**，属中等改动而非小修。
- 验收判据（可机器判定）：两个窗口下该行都须为 `RAW_MQandDP=28800,8` 且 `SB=0,0,3,3`，
  且窗口 A 的 48 行保持与 GATK 逐字节相同；另需复核 `10020421/10020431` 的第三种错误状态。

**投入产出提示**：D2 自第 4 轮起已消耗 5 轮，结论是「小修不够」。若继续不划算，
Track B 的**退出码分歧**（某 `--alleles` 记录经 Event 最小化后 ALT 为空时，GATK 以 exit 3
中止而 native exit 0 静默成功）形状最收敛、验收只需比对退出码，是更划算的替换目标。

**一处方法学教训**：状态文档曾把「流式路径已传实际 read batch（非空指针）」当作已修正的证据，
但"非空"不等于"范围正确"。修复必须由 oracle 判定，不能由改动描述判定。
同样地，第 5 轮我引用的 agent 推荐修法（靶子为 `*` 等位基因）未经独立验证即被采纳进计划，
第 6 轮证明靶子错了 —— **agent 的结论必须先用插桩/实验独立复核再作为行动依据**。

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

## 下一步执行顺序（第 12 轮更新；修复仍**未落地**）

1. ~~结束并确认 Serial 构建，完成 Serial oracle。~~ 已完成。
2. ~~跑完整 HC/Mutect2 回归。~~ 已完成；最新为第 10 轮对 commit `03b02b7` 的
   双后端全量 **280/280**（`.diag/regression/20260911-010134/`），套件规模在第 11 轮
   注册诊断 oracle 后为 281（OpenMP）。
3. ~~第 1 轮三路审计~~ 全部完成，报告在
   `fastgatk-native/evidence/2026-09-10-parallel-audit/`（Track B / C / D）
   与 `fastgatk-native/evidence/gvcf-stream-by-region-divergence-20260910.md`（Track A）。
4. **落地修复（唯一未完成的交付）**。当前状态与建议：

   - **D2（`-L` 窗口依赖的 `RAW_MQandDP`/`SB`）**：根因、产生点、验收判据都已定案，
     但**方案 1（merge 处成组去重）已被验收判据否决**（窗口 A 0 → 4 行回归、
     窗口 B 未修复），不应再试。**唯一推荐方向是方案 2**：渲染/注释处优先采用
     **已配对孪生**，且 **strand 路径必须一并采用孪生的 likelihood 行**
     （只换 context ordinal 治不了 `SB`）—— 属**中等改动**，不是小修。
     验收判据：两个窗口（`20:10019901-10020710`、`20:10020381-10020710`）下该行都须为
     `RAW_MQandDP=28800,8` 且 `SB=0,0,3,3`，且窗口 A 的 48 行保持与 GATK 逐字节相同；
     另需复核 `10020421/10020431` 的第三种错误状态。详细交接见上述 evidence 文件
     「第 6/7/8 轮」各节。
   - **投入产出提示**：D2 已消耗 5 轮且结论是「小修不够」。若继续不划算，
     **Track B 的退出码分歧**（`--alleles` 记录经 Event 最小化后 ALT 为空时，
     GATK exit 3 而 native exit 0 静默成功）形状最收敛、验收只需比对退出码，
     是更划算的替换目标；注意仓库内**没有现成的 `makeMinimalRepresentation` helper**。
   - 优先级 2：流式 D1（凭空 phasing，`hc_call.cpp:4085-4195`）与 D3（参考块粒度，
     tile 内逐 locus RCM 取值）。
   - 优先级 3：Track B 的 `--alleles` 注入缺失（重叠 feature event 丢行）。
   - 每步都要：改 → 双后端增量构建 → 定向 oracle + 72 子集 → 全量回归 → commit。
     任一判据不过立即 `git checkout` 回退（第 8 轮已示范该纪律）。
5. **诊断 oracle 已注册**：`fastgatk-hc-gvcf-stream-overlap-diagnostic`（第 11 轮）。
   它是 diagnostic 形态（永远 exit 0）；**修复落地后必须把它的比较改成硬断言
   （streamed == non-streamed == GATK）并加 `FASTGATK_REQUIRE_GATK_ORACLE=1`**，
   否则该缺陷会一直以「绿色但错误」的形式存在。
6. 补齐静默 skip 的 32 个候选测试（补 `FASTGATK_REQUIRE_GATK_ORACLE=1`）与 tool audit 口径，
   使百分比与证据同源（见「待处理的口径问题」）。
7. **委派注意**：子代理必须显式用 `provider=deepseek-official`
   （写 `deepseek` 会直接失败；`subagent` 工具不暴露 provider，需经 `workflow` 覆盖）。
   且 **agent 结论必须先用插桩/实验独立复核再作为行动依据**（第 5 轮已因违反此条而靶子出错）。
8. 只有当每个核心选项都有匹配范围的当前证据时，才讨论「1:1 完成」。


## 当前状态与下一步（第 8 轮更新，取代上面更早的「下一步执行顺序」）

### 已落地并上锁

| 项 | 状态 | 守护 |
| --- | --- | --- |
| D2（`-L` 窗口依赖的 `RAW_MQandDP`/`SB`）双倍体路径 | **已修复** | `fastgatk-hc-window-invariance-gatk-oracle`（strict，12 窗口） |
| D2 同类第二实例：非双倍体（`sample_ploidy != 2`）分支 | **已修复** | `fastgatk-hc-ploidy-window-invariance-gatk-oracle`（strict，7 窗口，ploidy 3） |
| 双后端全量 | **283/283**（OpenMP 1039s / Serial 1024s） | `.diag/regression/round7-full.log` |

改动位置：仅 `fastgatk-native/src/hc_call.cpp`。机制：记录被从「从未为该等位基因请求
PairHMM、context ordinal 为哨兵 `UINT32_MAX`、likelihood 行为 `-inf`」的 AssemblyRegion owner
注释；修复为优先选用携带同 Allele 且 context 非哨兵的 **twin owner**（`owner_has_pairhmm_context`）。

### 分类收敛：D1/D3 不是缺陷

判别实验证明 pinned GATK 自己在相应窗口下就产生同样的 phased 行与同样的退化块，
native 逐字节相同 —— 属 `--stream-by-region` 的**范围/语义问题**（GATK 输出本身依赖 `-L`，
流式按 tile 分别求值再拼接，原理上无法等价）。待办：显式定义其语义并记录为已文档化分歧；
`fastgatk-hc-region-streaming-contract` 的「流式 key 集合 == 非流式」断言原理性不可满足，需改写。

### 缺陷形态普查（已建立可判定判据）

**判据**：哨兵候选 ⇔ 从未进入 PairHMM 请求表 ⇒ 其 likelihood 行为 `-inf`。
因此**消费 likelihood 行的探针在该 owner 上必然失败（安全）**；
**只做下标越界检查的探针会「成功」并静默降级（易感）** —— 已知实例为
`calculate_output_variant_annotations` 的 `has_rows`（`calling_pipeline.cpp:13811-13815`）。

| 范围 | 站点数 | 易感 | 备注 |
| --- | ---: | --- | --- |
| HC（`hc_call.cpp`） | 15 | **仅剩 1 处**：`:3521` ordinary-VCF `annotations_from_owner` | 可达性未测得（该 guard 在当前夹具上不激活，需先造 multi-ALT ordinary-VCF 夹具）；另有 R-A2/R-A3 残差（已修站点偏好组内部的 first-success 递归） |
| Mutect2（`mutect2_tool.cpp`） | 14 | **6 处理论易感**（最重：`somatic_annotation_counts` @4116，影响 F1R2/F2R1/SB/MBQ/MFRL/MMQ/MPOS 并传播进 FilterMutectCalls） | **未能在任何 Mutect2 夹具上复现**（44/44 行与 GATK 一致、590 owners、跨 5 个 `-L` 起点与 `--force-active` 不变）——**属潜伏，非已证缺陷** |
| `filter_mutect_tool.cpp` | 0 | — | 无 owner 循环 |

报告：`fastgatk-native/evidence/2026-09-11-wave0/`（`round-g1`/`round-g2`/`round-g3`）。

### 下一步（按价值/成本）

1. **Track B 两项**（唯一剩下的已证真 bug）：`--alleles` 重叠 feature 丢行；
   退出码分歧（**必须加 region 感知门控** —— 文件级 eager 移植会在「空等位记录落在 `-L` 之外」
   时过度中止，制造出 GATK 没有的新分歧）。
2. **补 `hc_call.cpp:3521`**：先造能激活其 guard 的夹具（multi-ALT 的 ordinary VCF），否则无法验证。
3. **消掉 R-A2/R-A3 残差**（已修站点偏好组内部的 first-success 递归）。
4. **契约层**：定义 `--stream-by-region` 语义并改写那条原理性不可满足的断言（不碰生产代码，可并行）。
5. **Mutect2 易感站点**：仅在出现 Mutect2 分歧时优先查这 6 处；当前无证据支持改动。

### 顺带记录（G3 发现，与本形态无关）

`fastgatk-mutect2 --stream-by-region` 在 1 Mb contig 上对 100000/50000/20000/5000 全部 tile 尺寸
均以 `RESOURCE_EXHAUSTED: single-base Mutect2 tile exceeds safe memory budget`（exit 2）中止；
10 kb 窗口（size 1000）下流式与非流式一致。即 **Mutect2 的区域流式在 contig 规模上不可用**，
这是一条待评估的容量/边界问题（未定性为缺陷）。
