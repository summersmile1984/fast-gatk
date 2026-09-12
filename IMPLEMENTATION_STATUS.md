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
下表第 1 行是**当前树**；其余各行是历史记录（数字不同只是因为当时的树不同，不是回归）。

| 范围 | OpenMP | Serial | 备注 |
| --- | --- | --- | --- |
| 全量回归（第 57 轮，当前树） | **309/309 通过**（1473.5s） | **309/309 通过**（1480.8s） | commit `9a2e3a5`（QD 分子用未取整 double），两后端均已重建；零陈旧告警；证据 `.diag/regression/20260913-052529/` |
| 全量回归（第 56 轮） | 309/309 通过（1503.4s） | 309/309 通过（1515.2s） | `5302ea8`（PL 归一化修复） |
| 全量回归（第 54 轮） | 308/308 通过（1482.3s） | 308/308 通过（1458.1s） | `0e19bbe` + 块起始坐标 REF 修复 |
| 全量回归（第 53 轮） | 308/308 通过（1495.1s） | 308/308 通过（1466.0s） | `00a02d0` + dense `-L` 裁剪 |
| 全量回归（第 52 轮） | 308/308 通过（1504.4s） | 308/308 通过（1496.3s） | `3389e98` + dense 跨位点物化 |
| 全量回归（第 51 轮） | 307/307 通过（1458.9s） | 307/307 通过（1422.5s） | `a0a30b7` + star-only 拒绝修复 |
| 全量回归（第 50 轮） | 306/306 通过（1497.9s） | 306/306 通过（1538.3s） | `640678b` + depth-gate 修复 |
| 全量回归（第 49 轮） | 305/305 通过（1566.2s） | 305/305 通过（1531.5s） | `56e2ea6` + reverse-trim 修复 |
| 全量回归（第 48 轮） | 304/304 通过（1622.6s） | 304/304 通过（1444.9s） | commit `0862251`，工作树未提交变更 0 项 |
| 全量回归（第 44 轮） | 303/303 通过（1431.4s） | 303/303 通过（1453.3s） | commit `c971cf3`；此后生产代码又改了 2 次（`dce568f`、`f0a8277`），故必须重跑 |
| 全量回归（第 42 轮） | 302/302 通过 | 302/302 通过 | commit `b3cd293` |
| 全量回归（第 10 轮复验） | 280/280 通过（1019s） | 280/280 通过（1017s） | 对 commit `03b02b7`、二进制由 pristine 源重建后的复验；证据目录 `.diag/regression/20260911-010134/` |
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

## 待办：forced-alleles 多等位/`*` 发射分歧（第 10 轮已完整规格化，**刻意未修**）

**分歧**：`--alleles` 强制等位基因时，GATK 会发射一条含符号 `*` 的多等位行并只打 LowQual 过滤，
native 完全不发射该行。最小复现（已写成严格门禁，当前**按设计失败**）：

```text
GATK   : 697 TC>T 98.60 .  |  698 C *,A 0 LowQual  |  700 CA>C 144.77 .
NATIVE : 697 TC>T 98.60 .                          |  700 CA>C 144.77 .     (缺 698 行)
```

三条对照（drop-alleles / 仅强制删除 / 仅强制删除+drop-alleles）两侧均逐字节相同 → 归因确定。

**GATK 语义已从 pinned 源码证实**（非采信转述）：
- `GenotypingEngine.java:167-170`：`... && forcedAlleles.isEmpty()` —— 这是该文件里**唯一**的
  强制等位基因豁免；`forcedAlleles` 取自 `AssemblyBasedCallerUtils.java:1005-1009`。
- 不丢弃而是降级：`GenotypingEngine.java:183-186` → `QUAL=0` + `FILTER=LowQual`。
- 符号 `*` 的保留：`GenotypingEngine.java:304-327` 的 `calculateOutputAlleleSubset`
  （`isSpuriousSpanningDeletion` 314；`toOutput` 316 含 `forcedAlleles.contains(allele)`）。

**定位**：`calling_pipeline.cpp:17521-17533` 的无条件置信度门
（trace：`pos=697 ref=C alt=A support=0 qual=0 threshold=30 decision=confidence_suppressed`）。
native 该候选的 QUAL 已经是 **0**，即 GATK 的值——缺的只是发射。

**为何刻意未修（5 项必须协同的改动）**：
1. 发射豁免（依据上述 `forcedAlleles.isEmpty()`）；
2. **FILTER 列**：LowQual 无任何发射通路，7 处 writer 记录点全部硬编码 `'.'`
   （`hc_call.cpp:3615/3624/3862/4587/4598/4917/4945`），而 `GenotypeCall` 结构里没有 filter 字段；
3. **符号 `*` ALT**：ordinary diploid writer 明确关闭（`hc_call.cpp:3238/3248`），只有 gVCF 路径发射；
4. **Number=G/R 重排**到 GATK 的 `[REF,*,ALT]` 顺序（native 为 `[REF,concrete ALT,*]`）；
5. **对已发布等位基因划分做注释重算**——实测该候选**完全没有注释**
   （`FASTGATK_DEBUG_ANNOTATION_POSITION=697` 下 native 不打印 summary，
   因为 `calculate_variant_annotations` 位于被 `continue` 跳过的发射分支内），
   而该行需要 10 个精确 INFO 字段（含依赖 `*` 划分的 `ReadPosRankSum=3.523`）。

**结论**：这是**一个成规模的功能缺口**，不是最小修复——因此本轮只交付门禁与规格，不动代码。
门禁脚本 `fastgatk-native/scripts/verify_hc_forced_alleles_emission_gate_oracle.py` 已入库，
**尚未注册进 CTest**（strict 形态会故意使套件变红；如需在套件中可见，应以 diagnostic 形态注册）。

**同时记录一条新观察（已测、未修、未被门禁覆盖）**：不带 `--alleles`、用 12 条带真实
`chr1:698 C>A` 的读加 6 条删除读时，两侧都发射双等位行，但
`GATK QUAL/QD = 282.04/5.88` vs `NATIVE 291.07/6.06`（其余 ALT/INFO/FORMAT 相同，`chr1:697` 逐字节相同）。
它与本缺口**不同源**（不是符号 ALT 问题）；「与 native ordinary diploid 路径强制 concrete-only
导致 AF 排除 `*` 同根」只是**猜测**，需独立夹具与门禁验证。

## 已文档化分歧清单（契约层，第 13 轮）

以下两项**不是缺陷**，而是 GATK 语义与 native 扩展之间需要显式约定的边界。判定依据：
判别实验证明 pinned GATK 自己在相应窗口下就产生同样的输出，native 逐字节相同。

1. **`--stream-by-region` 不能与非流式逐字节等价**（原理性，非实现缺陷）。
   GATK 的 gVCF 输出本身依赖 `-L` 窗口（参考块粒度、注释值、乃至是否定相都会变），
   而流式按 tile 窗口分别求值再拼接，因此「对用户 `-L` 跑一次」的结果原理上无法复现。
   tile 内部计算是正确的：其输出与「对该 tile 的 halo 窗口跑一次非流式 GATK」逐字节相同
   （实测：`fixtures/chr20/mnp.bam` + `-L 20:10019901-10020500`，GATK 与 native 17 行、0 差异）。
   **处置**：`fastgatk-hc-region-streaming-contract` 中「流式 key 集合 == 非流式」
   已加注释限定为**夹具范围内的回归守卫**，并附反例与证据指针（48 → 20 行）。
   该断言未被削弱（对其夹具仍然成立且仍有守护价值），但不再被当作该开关的 parity 声称。
   后续若要把 `--stream-by-region` 变成可声称 parity 的特性，需先定义其语义
   （例如「等价于对每个 tile 窗口分别跑 GATK 后拼接」）或放宽内存上界、放弃分块。
2. **D1 凭空 phasing / D3 参考块粒度**：同属上述窗口依赖的表现，非独立缺陷（同上证据）。

> 判据（可复用）：**任何 parity 声称都必须连带 `-L`（以及全部其它输入参数）一起固定**；
> 「native 与 GATK 1:1」只能表述为「在相同输入下逐字节一致」。

## 测试注册债（第 17 轮发现，待分诊）

对 `fastgatk-native/scripts/verify_*.py` 与 `CMakeLists.txt` 做了一次机械比对：
**266 个脚本中 248 个已注册，18 个未注册**。这 18 个**多数是历史遗留**（基准/辅助/长期存在），
并非本会话产生；不应批量注册（有些不是测试，有些可能已在别处注册）。分诊原则：

- **必须注册**：本会话新增的**严格门禁**——已注册 9 道（window-invariance、ploidy-window-invariance、
  alleles-overlap、span-del-qual、gvcf-symbolic-prior、arbitrary-ploidy-span-del-prior、
  polyploid-gvcf-span-del-prior、spanning-prior-genotype-gq、**af-zero-format**）。
- **刻意不注册为 strict**：`verify_hc_forced_alleles_emission_gate_oracle.py`
  —— 它按设计**必须失败**（记录尚未修复的符号 ALT/LowQual 发射缺口）；
  若要在套件中可见，应以 diagnostic 形态注册，否则会把套件永久变红。
- **已注册（第 18 轮分诊完成）**：`verify_mutect2_recheck_normal_replay.py`、
  `verify_mutect2_recheck_assembly_resultset_joint.py`、
  `verify_hc_multialt_owner_annotation_fixture_oracle.py`（三者在双后端实测 exit 0 / status=pass）。
- **诊断形态、禁止注册为 strict**：`verify_hc_alleles_deep_boundary.py` 与
  `verify_hc_alleles_deep_limits.py` —— 实测 **exit 1 / status="diverge"**，
  即这两份脚本仍在记录**尚未修复的 `--alleles` 分歧**（Track B 当初报了 5 处，
  本会话只修掉了其中的重叠强制事件那一处）。这两份脚本是这些残留分歧的**活证据**，
  应作为后续 `--alleles` 轨道的起点，而不是被注册成会变红的门禁。
- **很可能是辅助/基准而非测试**：`verify_native.py`、`verify_kernel_benchmark.py`、
  `verify_kokkos_backend_matrix.py`、`verify_bam_intervals.py`、
  `verify_pairhmm_results_oracle.py`、`verify_somatic_*.py`、`verify_fragment_aggregation_gatk_oracle.py`、
  `verify_reblock_gatk_multisample.py`、`verify_hc_multiallelic_gatk_oracle.py`
  —— 需逐个确认是「漏注册」还是「本就不属于套件」。

> 教训：本会话有两次新建了严格门禁却忘了注册（AF-format 直到本轮才补上；
> 另有一批由子代理创建的脚本同样没注册）。**门禁不注册等于没有门禁** ——
> 每次新增 oracle 后应立刻做一次上面的机械比对。

## 测试契约审计（第 20 轮）：6 处**被测试钉死的 parity 分歧** + 一个隐性失真机制

对 273 个 add_test 条目做程序化枚举（**第 20 轮当时的数字**；此后本会话又新增 24 个
条目，当前 `fastgatk-native/CMakeLists.txt` 共 283 个）：94 个名字不含 `oracle`；其中约 30 个按构造确属 oracle，
**约 64 个是真·native 内部契约**。全部 94 个都进了分类表（测试名 / CMake 行 / 脚本 /
是否做记录级断言 / 类别 / 行号）。报告：`fastgatk-native/evidence/2026-09-11-wave0/round-testcontract-audit.md`。

### 新测出 6 处 class-(c)（测试把**偏离 GATK 的行为**当契约；均以 pinned GATK 跑测试自带夹具实测确认）

| 测试断言 | 钉住的行为 | GATK 实测 |
| --- | --- | --- |
| `verify_reblock_gvcf.py:64/:82-83` | `--drop-low-quals` 下多出 POS 20 的参考块（第 3 条记录） | 只输出 2 条，**整条丢弃**该位点（`ReblockGVCF.java:542-545`） |
| `verify_genotype_gvcf.py:391` | 存在 `*,G` 跨接删除记录（`len(star_records)==1`） | 默认 **0 条**；`-all-sites` 下 ALT='.' GT='./.'（`GenotypeGVCFs.java:326-330`） |
| `verify_reblock_gvcf.py:283-287` | 三倍体保持 `0/1/1` + 压缩成 10 项 PL | 重分块为 `0/0/0`、`END=30`、4 项 PL（`ReblockGVCF.java:515-531`） |
| `verify_select_variants.py:820` | `ref_only_records == []`（丢弃全 hom-ref 记录） | **保留**，ALT='.' |
| `verify_variant_recalibrator.py:84/422/437/456` | `culprit=full-covariance-gmm` | 该字符串在 pinned jar 里**出现 0 次**；GATK 写 `culprit=MQ` |
| `verify_variant_filtration.py:504-505` | `AS_FilterStatus=LowASQD,PASS` | `AS_FilterStatus=SITE\|SITE`，且不施加等位基因过滤 |

另有实测分歧：`verify_filter_mutect_calls.py`（MBQ/MMQ 元数不符；`FILTER=FAIL`/`AS_FilterStatus=low_tlod`，
而 **`--min-tlod` 根本不是 GATK 选项**、`low_tlod` 不是 GATK 常量）、
`verify_depth_of_coverage_multisample.py`（多出 per-sample 列与 1.00 分母 vs GATK 单一列与 1.50）、
`verify_mutect2.py:650`（缺 OTHER_NORMAL 样本列）、以及 4 处输入校验包络分歧
（负数 `-ip`/`-ixp` 拒绝、接受未索引 gzip 区间、3 字段 `.interval_list`）。

### 更大的问题：一个**让绿色失真的机制**

**11+ 个脚本**把 GATK 比对包在 `if java.exists() and jar.exists():` 里，**却从不断言 jar 存在**。
一旦 oracle 缺失，这些 (a) 字面量就**静默降级**为「未经验证的 native 自身期望」，而
**测试仍然保持绿色**。这与 verify_indel.py:173 那次事故同形，但**不可见**——因为它不会变红。

> 结论：**当前的 294/294 绿色并不能支撑 parity 声称**。其中包含至少 6 处被钉死的分歧，
> 且有一批「看起来是 GATK oracle、实际不是」的测试在静默退化。
> 这条应优先于「再修一个 bug」处理，因为它影响的是**所有后续结论的可信度**。

### 处置顺序（建议）

1. **修隐性失真**：让这些脚本在声称做 GATK 比对时**必须**确认 jar 存在（缺失则明确 skip 或失败），
   使隐性降级变为显性。此项不动生产代码，但会暴露一批需要决断的用例。
2. **逐条处理 6 处 class-(c)**：每处都是真 parity 缺陷，需与对应测试的过期断言一起修
   （参照第 19 轮 `verify_indel.py` 的正确处置方式）。
3. **处置 class-(b) 高危项**：`verify_hc_genotype_priors.py:128-184` 的精确 PL/GQ 字面量
   （native PairHMM 输出被当作期望，无 GATK 比对）、以及多处**按位置读 FORMAT 列**的脆弱写法
   （`verify_reblock_gvcf.py:77,80`、`verify_select_variants.py:116,790-795`、
   `verify_left_align.py:111-113,326-327` 等）。

> 注：`verify_indel.py` 经本轮实测确认**已正确**——该夹具上 native 与 GATK 三个 gVCF 记录逐字节相同，
> 故第 19 轮改成 GATK 契约的断言是对的。

## Oracle 完整性（第 21 轮）：消除「缺失 oracle 被误当证据」的机制

审计的估计被实测**放大了 16 倍**：**176 个已注册脚本 / 294 个已注册测试中的 187 个（63.6%）**
把 GATK 比对包在一个**从不断言 jar 存在**的存在性守卫里；其中 **51 个连
`FASTGATK_REQUIRE_GATK_ORACLE` 逃生口都没有**。（审计原文估计「11+」，实测 176。）

**关键限定：今天没有一个是真未验证的。** 每个脚本计算出的 oracle 路径都解析到
`third_party/jdk17/bin/java` 与 pinned jar，`paths_missing` 全空 —— 也就是说
**风险是潜在的（latent），不是已发生的（active）**。这一点很重要：它意味着此前的
294/294 并没有被实际污染，但**随时可能**。

处置：
- 新增共享助手 `fastgatk-native/scripts/oracle_guard.py`（复用已有的 `oracle_toolchain.py` 解析器），
  策略为：oracle 在 → 行为不变；oracle 缺失且 `FASTGATK_REQUIRE_GATK_ORACLE=1` → **响亮失败**
  （非零退出并指名缺失路径）；oracle 缺失且未设该变量 → 打印
  `[NOT VERIFIED AGAINST GATK]` 横幅并 exit 0（不得被误读为已验证）。
- 应用于全部 176 个已注册脚本（`+369/-42`；被删的 42 行全是纯存在性守卫表达式，
  **无任何期望值或断言被改写**）。7 个未注册脚本未动（CTest 不运行它们）。
- **`run_regression.sh` 现在默认设置 `FASTGATK_REQUIRE_GATK_ORACLE=1`** —— 这是关闭
  「横幅在正常运行中不可见」这一残余缺口的关键：`--output-on-failure` 会丢弃通过用例的输出，
  故横幅只在 `-V` 或直接运行时可见；改为由强制门禁兜底。

验证（均实测）：
- 套件在**设与不设** `FASTGATK_REQUIRE_GATK_ORACLE=1` 两种情况下均
  **294/294 双后端通过**（1195.7s/1198.2s 与 1194.7s/1183.8s）——
  运行时证明没有任何已注册测试缺 oracle，也证明该变量今天即可用作 CI 门禁。
- 「咬合」演示：把一个已注册脚本的 oracle 隐藏（临时目录镜像仓库布局、jar 目录置空），
  (i) 设变量 → exit 1 并报 `[ORACLE REQUIRED]` 与缺失路径；
  (ii) 不设 → exit 0 但 stderr 打印 `[NOT VERIFIED AGAINST GATK]`，stdout 标 `"java_oracle": false`；
  (iii)/(iv) oracle 存在时两种情形输出逐字节相同（即改动是惰性的）。
  并**经真实 CTest 沙盒复验**：`JAVA=bogus` + 变量 → 测试 FAILED；仅 `JAVA=bogus` → Passed（横幅需 `-V` 才可见）。

仍未证明（报告 §4）：本轮的保证是「**缺失 oracle 会响亮**」，**不是**「其余断言已由 GATK 验证」——
审计 §4.4 的发现（native-only 字面量、按位置读 FORMAT 列、SelectVariants 作为读取器）仍待处理；
只扫描了 java/jar 存在性守卫，未覆盖其它静默跳过机制（裸 except、吞掉非零返回码）；
`FASTGATK_REQUIRE_GATK_ORACLE` 运行只能证明**到达了守卫**的脚本找到了 oracle，
被更早 return 挡住的守卫不会被触及。

## 会话累计（第 26 轮更新）：15 个已证真 bug 已修并各自上锁

双后端全量 **299/299**（在 `run_regression.sh` 默认强制 `FASTGATK_REQUIRE_GATK_ORACLE=1`、
零陈旧告警的前提下取得）。已注册严格 GATK 门禁 **18 道**。

| # | 缺陷 | 守护门禁 |
| --- | --- | --- |
| 1 | D2 双倍体：`-L` 窗口依赖的 `RAW_MQandDP`/`SB` | `hc-window-invariance-gatk-oracle` |
| 2 | D2 非双倍体实例 | `hc-ploidy-window-invariance-gatk-oracle` |
| 3 | Track B：`--alleles` 重叠强制事件被抑制 | `hc-alleles-overlap-gate-oracle` |
| 4 | spanning-deletion AF 先验（`*` 应取 SNP 先验） | `hc-span-del-qual-gatk-oracle` |
| 5 | gVCF 参考置信度 `*` 先验 | `hc-gvcf-symbolic-prior-gatk-oracle` |
| 6-8 | `*` 长度先验同类三处 | `arbitrary-ploidy-span-del-prior` / `polyploid-gvcf-span-del-prior` / `spanning-prior-genotype-gq` |
| 9 | AF/MLEAF 零值格式化（`0.000` vs `0.00`） | `af-zero-format-gatk-oracle` |
| 10 | gVCF 具体变异记录误带 INFO/END | `gvcf-indel-end-gatk-oracle` |
| 11 | BQSR `culprit` 应为「最差注释名」而非模型来源串 | `variant-recalibrator-culprit-gatk-oracle` |
| 12 | SelectVariants ref-only 记录保留（依 `--exclude-non-variants`） | `select-variants-refonly-gatk-oracle` |
| 13-14 | VariantFiltration 等位基因过滤语义（含 flag 路径的 FILTER 列） | `variant-filtration-asfilterstatus` / `variant-filtration-flag-only` |
| 15 | ReblockGVCF `--drop-low-quals`：再基因分型须**先于** `--rgq-threshold` 转换 | `reblock-gvcf-droplowqual-gatk-oracle` |

### 测量完整性（本会话新增的保障）

- 静默 skip 审计：117 个可 skip 脚本中的 32 个候选实跑，0 skip。
- **oracle 完整性（第 21 轮）**：176 个已注册脚本 / 187 个测试曾把 GATK 比对包在**从不断言 jar 存在**的
  守卫里（其中 51 个无逃生口）；实测当天均未失真（路径全部解析到 vendored jar），属**潜在**风险。
  已以共享助手 `oracle_guard.py` 改为 fail-closed，并让 `run_regression.sh` 默认强制
  `FASTGATK_REQUIRE_GATK_ORACLE=1`。
- `run_regression.sh` 陈旧性检查原先与 `fastgatk-hc-call` 单文件比较，会把 HC 不链接的源文件
  误报为陈旧（假告警会训练人忽略告警）；已改为与构建目录中**最新**产物比较。

### 已确认仍未修（按优先级）

1. **`verify_reblock_gvcf.py` 三倍体重分块**（case B）：**分歧确认为真**（GATK 重分块为
   `0/0/0 END=30` 4 项 PL，native 为 `0/1/1` 10 项），修复**已实现并用 pinned GATK 验证过**，
   但**刻意回退**——因为落地它会使同一测试文件中另一块（trim/gap/NON_REF-AD，`:441-480`）**无法满足**，
   而那块是三个 native 特性（reverse trimming、deletion gap、NON_REF AD）**唯一**的覆盖。
   建议路径：与其一起修正该夹具的 PL 向量（其 `PL[0]=0` 使 GATK 也重分块）与 INFO-DP 规则，
   并重新推导约 25 个值。门禁 `verify_reblock_gvcf_triploid_gatk_oracle.py` 已存在但**刻意未注册**
   （strict 会红）。
2. **`verify_genotype_gvcf.py:391`**：native 发射 `*,G` 跨接删除记录，GATK 默认 0 条
   （`-all-sites` 下 `ALT='.'`）——仍未修。
3. **forced-alleles 符号 `*`/LowQual 发射缺口**：需 5 项协同改动（FILTER 列通路、符号 `*` ALT、
   Number=G/R 重排、注释重算等），属成规模功能缺口；门禁
   `verify_hc_forced_alleles_emission_gate_oracle.py` 已存在但**刻意未注册**。
4. **`gvcf-max-alt-alleles-1` 参考块 PL/GQ**：定位到强制等位基因下 plausible-indel 信息量分类，
   4 个对照已建立；因涉及 RCM 与 haplotype realignment 管道而非最小改动，未修。
5. **零杂合度参数校验分歧**：GATK 接受 `--heterozygosity 0` 而 native exit 2；
   放宽校验非最小改动（GATK 自身在该值下退化），且会激活已知未修的 `estimate_mle_allele_counts` 回退分歧。
6. VariantFiltration 后续：flag + 无表达式 + 无 mask 被 GATK 接受但 native CLI 拒绝；
   mask/cluster 区间保真（indel/符号未证）；JEXL 裸标识符与 null 方法语义；
   native 的 `START` 别名、`CHROM`/`FILTER` 属性、`QUAL='.'` 按 `-10` 处理等既有 flag 无关缺口。
7. 各处已测但**刻意不门控**的子分歧：SelectVariants 的 FORMAT/INFO 键序与 AF 精度、
   ReblockGVCF 的 QUAL/键序/ref-block FILTER、BQSR 60 记录夹具的模型数值差（55/60 行逐字节相同）。
8. **Mutect2**：14 个 owner 选择站点中 6 处理论易感，但**未能在任何夹具上复现**（44/44 行与 GATK 一致、
   590 owners、跨 5 个 `-L` 起点稳定）——属潜伏风险，非已证缺陷。
9. ~~**GenotypeGVCFs 缺 GATK `regenotypeVC` 的 `INFO/DP > 0` 前置条件**~~
   **已修（第 50 轮，第 32 个已修 bug）**：默认模式不再输出 GATK 不会重新基因分型的位点；
   门禁 `fastgatk-genotype-gvcf-depth-gate-gatk-oracle`。dense 模式的直通形状仍属未修的物化缺口。
   同根的**反向**残差（输入无 `INFO/DP` 时 GATK 会发布算出的 `INFO/DP` 与 QD，native 两者都不写）
   **仍未修、无门禁**，见第 50 轮节末。
10. **GenotypeGVCFs 对具体记录透传 `INFO/END`**（第 49 轮新测得，已亲自复测）：GATK 因
   `stop == start + REF.length() - 1` 而**不写** END，native 原样透传；trim 生效时该 END 还会过期
   （span 2-4 而 `END=5`）。属**既有**分歧（裁/不裁/守卫三组对照均分歧）。无门禁。

### 方法论沉淀（可复用）

1. **判据化**：哨兵候选 ⇔ 未进 PairHMM 请求表 ⇒ likelihood 行 `-inf`；故「消费行 ⇒ 免疫，只查下标 ⇒ 中招」。
   据此普查 HC 15 站点 / Mutect2 14 站点。
2. **根类化**：GATK「按等位基因长度选常量」vs native「固定值」——一个根类贡献 5 个 bug；
   从「一处」推到「一类」再推到「修干净」。
3. **测试契约审计**：发现 6 处测试把**偏离 GATK 的行为**当契约（已修 4 处）。
   处置判据：修复 parity 缺陷若使既有测试变红，先判断该断言固定的是 **GATK 行为**还是 **native 自身行为**；
   后者是测试债，应随修复一起纠正——**不要因为红了就放弃修复**（第 19 轮 verify_indel.py 实例）。
4. **门禁不注册等于没有门禁**：每次新增 oracle 后立刻做一次「脚本 vs CMakeLists」机械比对。
5. **verifier 自身也会错**：本会话多次由委派方或我自己的复核推翻先前结论
   （拼接器假设、agent 的 `*` 靶子、我标反的数值方向、D2 的 twin 在另一个 owner、R2 的 fix_site）。
   故：**agent 结论必须先用插桩/实验独立复核再作为行动依据**；处方同样要被验证。
6. **夹具优先**：7 处「不可达」分歧靠夹具构造打开 3 处；关键是非显然的输入条件
   （`*` 需在 AF 矩阵但不在被判基因型；GQ 先验需 `USE_POSTERIOR_PROBABILITIES`）。
7. **bug-for-bug 对齐 + 标注意图推测**：GATK 的若干行为疑似其自身 bug（如 VariantFiltration 的
   无 INFO 拆分上下文、ReblockGVCF 的 5 参数 builder 副作用）。parity 项目应**逐 bug 对齐**，
   但**必须在文档中标明「这是测得行为、不是我们认同的设计」**。

## 第 27–35 轮增量（累计 23 个已证真 bug 已修并上锁）

第 26 轮之前见上一节；以下为之后新增的 8 个修复（编号 16–23）：

| # | 缺陷 | 要点 |
| --- | --- | --- |
| 16 | GenotypeGVCFs「输出等位基因子集塌缩到 REF」**崩溃** | 非字段差异而是中止（`invalid genotype PL remap dimensions`）；Kokkos 内核拒绝 `target_allele_count < 2`，与 `*` 无关（普通 `A G,<NON_REF>` 同样崩） |
| 17 | 孤立 `*` 剪除后存活具体 ALT 的**基因型投影** | 上一轮「按 PL 重推」假设**被证伪**：GATK 走 PREFER_PLS 第二支，投影**源基因型**（`bestMatchToOriginalGT`），保留拷贝顺序与相位、不赋 GQ |
| 18 | 跨接删除**归属判定边界** | `span.begin <= pos` 应为 `<`（`GenotypingEngine.java:369`）；native 的归属集合取自**输入记录**而非**已发射等位基因** |
| 19 | 「仅跨接删除」记录：已归属的 `*` 被**豁免 AF 阈值** | GATK 对 `*` 与具体 ALT 一视同仁施加 `passesThreshold()` |
| 20 | dense 模式 REF-only 行 `QUAL=Infinity` | GATK 值**确为 `+Infinity`**（Java `%f` 打成字面 `Infinity`）；native 算出了 `-Infinity` 补数却被 `std::isfinite` 守卫丢弃 |
| 21 | dense FILTER 列 `RGQ` vs `.` | **两层**：native 从不重建 FILTER 状态（继承源叶子）；且继承的字典索引对着**输出** header 解析，HTSlib 单字典 + 解析期自动注册 → 陈旧索引指向首个追加 id |
| 22 | 输出 header 缺 `##FILTER=LowQual` | GATK **无条件**声明（`GenotypeGVCFsEngine.java:416`） |
| 23 | 输入 `##FILTER=PASS` 行被剥离 | **真实原因是 HTSlib 注入合成 PASS 行并把输入那行当重复销毁**；纯朴修复（直接删剥离）**经实测否决**——会对未声明 PASS 的输入也输出 PASS |

### 新暴露的**整类**边界：header 级分歧（此前无任何门禁覆盖）

门禁只比数据行，因此整类 header 差异长期不可见。第 22/23 轮的门禁**先只隔离出 header**
（修前 22 处 header-only violation、数据行全绿）——这个现象本身就是盲区的证据。已完整记录：

- `##FILTER` 组的**绝对位置**（GATK 2 vs native 3），根源是**整体 header 排序**分歧
  （GATK = htsjdk 排序，native = 输入顺序 + 自有 INFO 排序）；
- `##contig` 位置（native 1 vs GATK 28）；
- `##GATKCommandLine`（GATK 写、native 不写）；
- GATK 额外写的标准行：`##INFO=BaseQRankSum/MQRankSum/ReadPosRankSum`、第二个 `##INFO=AD`/`##INFO=DP`、第二个 `##FORMAT=AD`；
- htsjdk 给每个 `Description` 加引号并改写 MLEAC/MLEAF 措辞。
- `##fileformat`(VCFv4.2) 与 `##source` 两侧一致。

> 意义：**「数据行一致 ≠ 文件一致」**。若要声称 1:1，header 层是必须单列的一条战线。

## 第 36–57 轮增量（累计 **39** 个已证真 bug 已修并上锁）

第 27–35 轮的 16–23 号见上一节。以下 24–31 号在此前各轮已修复并上锁，
但**只有提交信息与证据文件、没有进这份交接文档**——本节补上（交接债）。

| # | 缺陷 | 要点 / 证据 |
| --- | --- | --- |
| 24 | GenotypeGVCFs 输出 header **内容**分歧整类 | 15 缺 / 10 多 / 10 文本不符 → 0/0/0（`d4faf64`） |
| 25 | NDA（`--annotate-with-num-discovered-alleles`）注释缺失 | 一处守卫阈值（`628b5e4`） |
| 26 | 畸形 FORMAT 记录**静默截断整个输出** | 本会话最严重：`bcf_read` 的 -2（解析失败）被当成 -1（EOF）（`fe953b5`） |
| 27 | 同一根因的系统性普查与修复 | 27 文件 / 51 处（34 `bcf_read` + 2 `sam_read1` + 15 `hts_getline`）（`b3cd293`） |
| 28 | GATK 兼容 header 的整体顺序 | htsjdk 排序规则；改动前先评估影响面（`c971cf3`） |
| 29 | GenotypeGVCFs 自身的 `FILTER=LowQual` | 阈值判在**原始 double**上，不是打印出的 QUAL 令牌（`dce568f`） |
| 30 | 删除归属改为追踪**已发射**等位基因 | 判定为状态记账；顺带修出两处此前未知分歧（含反向）（`f0a8277`） |
| 31 | GenotypeGVCFs **反向 trim**（第 49 轮） | 见下节 |
| 32 | GenotypeGVCFs 缺 `INFO/DP > 0` 前置条件（第 50 轮） | 见下文「第 50 轮」一节 |
| 33 | 只含跨接删除的位点在默认模式必须被拒绝（第 51 轮） | 见下文「第 51 轮」一节 |
| 34 | dense 模式跨位点记录物化（第 52 轮） | 见下文「第 52 轮」一节；本会话最大的一处结构性缺口 |
| 35 | dense 记录供给未裁剪到 `-L` 区间（第 53 轮） | 见下文「第 53 轮」一节 |
| 36 | 块起始坐标的 REF 取自 FASTA 而非记录（第 54 轮） | 见下文「第 54 轮」一节 |
| 37 | 等位基因频率后验在**线性空间**下溢，导致真实调用被丢（第 55 轮） | 见下文「第 55 轮」一节；本会话最严重的一处丢调用 |
| 38 | 子集投影后的 PL 未按 htsjdk 归一到最小值 0（第 56 轮） | 见下文「第 56 轮」一节 |
| 39 | QD 分子用了取整后的 QUAL 而非未取整 double（第 57 轮） | 见下文「第 57 轮」一节；至此 chr20 语料三种遍历逐字节一致 |

## 第 49 轮：GenotypeGVCFs 反向 trim（第 31 个已修 bug）与三方独立复核

### 缺陷与修复

GATK 对每个重新基因分型的位点都会把等位基因裁到最小表示：
`GenotypeGVCFsEngine.java:160-169` 在 `:165 finalizeAnnotations` 之后、`:189` 站点注释之前调用
`:167 GATKVariantContextUtils.reverseTrimAlleles()` = `trimAlleles(vc, trimForward=false, trimReverse=true)`
（`:1443-1445`）。native 直接发布合并后的 REF/ALT，**未裁**，因此两个等位基因共享尾碱基的位点
会比 GATK 长若干碱基：

| fixture（模式） | GATK | native（修前） |
| --- | --- | --- |
| `2 AAAA AACA,<NON_REF>` | `chr1 2 . AAA AAC` | `chr1 2 . AAAA AACA` |
| `2 ACGTACGT ACGT,<NON_REF>`（ACGT 参考） | `chr1 2 . ACGTA A`（共同尾串 4，裁 3） | 原样 |
| `2 AAAA AACC,<NON_REF>` | 原样（无共享尾碱基） | 原样 |
| `2 AAA A,<NON_REF>` | 原样（ALT 长度 1 → `:1458` 守卫） | 原样 |
| dense `4 AA *,<NON_REF>` | `chr1 4 . A *`（`*` 原样） | `chr1 4 . AA *` |

修复：Host 侧新助手 `apply_gatk_reverse_trim()`（`genotype_gvcf_tool.cpp`），
在 `apply_gatk_output_allele_subset()` 的三个返回点调用。规则逐条对齐 GATK：
`:1458` 守卫判在**已发射**（子集之后）的等位基因表上；`:1462-1467` 候选 =
非符号且非 `*` 的等位基因（REF 也在内），尾串受**最短候选**封顶；
`:1469-1475` 封顶即意味着某个等位基因被吃空，此时**少裁一个**（等位基因永不为空）；
`:1497-1501` 符号与 `*` 原样复制；`:1515-1518` 重建记录时 **start 不变**、`stop = start + REF.length()-1`。

**放置约束（本轮的硬约束）**：必须在 `EmittedDeletions::record()` **之后**。GATK 在
`GenotypingEngine.java:178-179` 记录的删除尺寸取自 trim **之前**的等位基因，若顺序颠倒，
`AA`→`A` 会把 `*` 记录的 `deletionSize` 从 1 变成 0。

新门禁 `fastgatk-genotype-gvcf-reverse-trim-gatk-oracle`
（`scripts/verify_genotype_gvcf_reverse_trim_gatk_oracle.py`，7 个用例，双后端严格通过）：
覆盖 1 碱基 ALT 守卫、吃空少裁角、**与输出等位基因子集的先后顺序**（多把 1 碱基 ALT 先被子集剪掉，
守卫因此不生效）、长等位基因、以及 dense 下「唯一候选是 REF」裁到单碱基而 `*` 保持原样的用例。

### 三方独立复核（`workflow` + `deepseek-official`；报告在 `fastgatk-native/evidence/2026-09-13-round49/`）

1. **对抗性验证**（`reverse-trim-adversarial-verification.md`）：独立造了 104 个位点的宽度扫描 +
   约 35 个定向夹具 + 直接对 pinned jar 调用 `reverseTrimAlleles()` 的 Java 探针。
   **裁剪算法本身未被证伪**：0 个等位基因差异、0 处与独立 Python 转写不一致；
   守卫、吃空少裁、POS 不动、FORMAT 不变、符号/`*` 原样、放置约束（含刻意构造的破坏夹具）、
   以及**第二条遍历路径**（`--stream-by-locus`）全部成立。
2. **真实语料可达性**（`reverse-trim-reach-on-real-corpora.md`）：142 个真实 gVCF / 165,527 条记录中
   **该裁剪触发 0 次**（用 pinned jar 的 `reverseTrimAlleles()` 独立复核同为 0）。
   原因是 `:1458` 守卫：165,524 条记录带长度 1 的非 `*` 等位基因（其中 162,712 条是一碱基 REF），
   仅 3 条到达 `normalizeAlleles()` 且都不共享尾碱基；去掉守卫的反事实仍为 0。
   故既有语料**不可能**覆盖该分歧，本轮修复不会改变任何真实夹具的行——这是「零回归」的机制解释。
   该报告还用 5,000 条随机记录对真实 GATK 方法做了差分测试（2,990 条命中，0 处不一致）。
3. **多输入 header**（`multi-input-header-merge.md`）：见下节更正。

### 复核暴露的两处**既有**分歧（不是本轮修复引入的，已由主会话亲自复测确认）

| 现象 | GATK | native | 判定 |
| --- | --- | --- | --- |
| 输入 `INFO/DP=0`、`FORMAT/DP=20` 的变异位点 | 默认模式**不输出**；dense 模式输出 REF-only 直通行（`AAAA AACA . . DP=0`，`GT:AD ./.:0,20`） | 输出**已定型的变异行**（`AAA AAC 92.64 … 0/1`） | **既有**：native 缺 GATK `regenotypeVC` 的 `INFO/DP > 0` 前置条件（`GenotypeGVCFsEngine.java:160`）。行的**有无**与 trim 无关，故修前同样分歧（修前该行是未裁的 `AAAA AACA`） |
| 输入 `INFO/END=5` 的具体记录 | 一律**不写** END（`stop` 等于隐含末端） | 原样透传 `END=5` | **既有**：在「不裁」（`AAAA AACC`）与「守卫不裁」（`AAA A`）两个对照上同样分歧，故与本轮 trim 无关；仅当 trim 生效时该 END 才**过期**（span 变 2-4 而 END=5） |

> 方法论：这两条先是委派方**对抗性**发现的，但它们的机制叙述（DP 门控来自 merger
> 的 DP 取值路径）是**推理**；主会话用三组对照（裁 / 不裁 / 守卫不裁，以及 `INFO/DP` 缺席）
> 亲自复测，确认「与 trim 无关」这一结论，并把「DP 门控机制」标为待独立验证的假设。
> **委派方的结论在进入文档前必须被这样分开处理：测得的部分可用，推理的部分要另证。**

### 关于「多输入 header 合并」的**前提更正**（重要）

原文写「GATK 合并全部 `-V` 的 header」。**这个前提是错的**，已实测：
`VariantLocusWalker.java:33-35` 把 `-V` 声明为单个 `String drivingVariantFile`，
pinned GATK 4.6.2.0 对第二个 `-V` 直接报
`A USER ERROR has occurred: Illegal argument value: Argument 'V/variant' cannot be specified more than once.`（exit 1，18 个用例全部如此）。
native 的「可重复 `-V`」因此是**超集**，对齐目标只能是 `CombineGVCFs`/`MultiVariantDataSource`
的合并语义（`VCFUtils.smartMergeHeaders` + `VcfUtils.getSortedSampleSet`）。据此实测出的真实分歧：

1. **样本列顺序**：GATK 按**字典序排序**样本列，native 按**输入顺序**（仅当样本名恰好有序时两者相同——属「假通过」）；
2. **后续输入的声明被丢弃**：native 的 header 是首个输入的 `bcf_hdr_dup`，后续输入只合并 contig 与样本名，
   故后续输入独有的 `##INFO`/`##FORMAT`/`##FILTER`/`##ALT` 全部丢失（对调 `-V` 顺序可对称复现）；
3. 输入独有样本在无数据位点：GATK 写 `./.`，native 写 `./.:.:.:.,.,.:.`。

**已被证伪（native 本来就对）**：`##source` 两边都不做并集（首个值胜出）；contig 声明**会**合并；
contig 顺序两边都跟随首个输入；contig 长度冲突两边都报错拒绝；同名样本两边都塌缩成一列；
共有声明的相对顺序 34/34 一致；输入 header 相同时最终 header 完全一致。

### 方法论补充（第 27–35 轮）

1. **纯朴修复必须先实测**：第 23 轮的「直接删剥离」看起来完全合理，实测却对未声明 PASS 的输入
   也输出 PASS。**能被测量否决的假设，必须在合入前测量**。
2. **门禁的覆盖面本身是风险**：「只比数据行」漏掉整类 header 差异（第 22 轮）。
   新增门禁时应问「它看不见什么」。
3. **一次修复可能解除另一处的掩盖**：第 19 轮修好 `*` 阈值后，第 18 轮标记的
   「输入记录 vs 已发射等位基因」结构性差异**立刻显现**——修复会暴露被掩盖的缺陷。
4. **`javap` 是可用证据源**：多轮用 `javap -c/-constants` 从 pinned jar 证实 GATK/htsjdk 行为
   （5 参数 builder 不拷贝 filters、`%f` 渲染无穷、字典单例）。这是介于「读源码」与「跑实验」
   之间的第三种证据。

## 第 50 轮：GenotypeGVCFs 的 `INFO/DP > 0` 前置条件（第 32 个已修 bug）

第 49 轮的对抗性复核发现的「native 把 GATK 根本不重新基因分型的位点也输出了」，
在本轮先建门禁、再按 GATK 规则修掉。

### 规则（源码逐行核对，非推理）

`GenotypeGVCFsEngine.java:157-174` 的**整个**重新基因分型块（等位基因子集、`:167` 反向 trim、
以及全部站点注释）只在「记录是变异 **且** 合并后 `INFO/DP > 0`」时才进入：

```java
if ( originalVC.isVariant() && originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY,0) > 0 ) {
    ... calculateGenotypes / finalizeAnnotations / reverseTrimAlleles ...
} else {
    result = originalVC;                       // :174 直通
}
```
且 `:181` 再次施加同一深度判据后才注释/发射，故在**非 dense** 模式下直通记录走到 `:198 return null`
——**该位点一条记录都不输出**。深度取**合并后**记录的值：
`ReferenceConfidenceVariantContextMerger.calculateVCDepth()`（`:352-360`）在存在 `INFO/DP` 时**只用它**
（不回退到基因型），否则对每个样本取 `getBestDepthValue()`（有 `MIN_DP` 用 `MIN_DP`，否则用 `DP`）求和，
并在 `:382` 仅在 `depth > 0` 时发布。

### 实测（主会话亲自跑 pinned GATK 4.6.2.0）

| fixture（默认模式） | GATK | native（修前） |
| --- | --- | --- |
| 变异位点 + `INFO/DP=0`、`FORMAT/DP=20`、PL 强杂合 | **无记录** | `chr1 2 . AAA AAC 92.64 … 0/1` |
| 同上但无 `INFO/DP`、`FORMAT/DP=0` | **无记录** | `chr1 2 . AAA AAC 92.64 … 0/1` |
| `INFO/DP=0` 而 `MIN_DP=20`（`calculateVCDepth` 优先 INFO/DP） | **无记录** | `… DP=0 … 0/1` |
| 边界：`INFO/DP=1`（对照） | 输出并（反向）裁剪 | 一致 |
| `INFO/DP=0` 且 PL 为纯合参考（对照） | 无记录 | 一致（子集已剪掉 ALT） |

dense 模式同一输入下 GATK 的直通形状是 `chr1 2 . AAAA AACA . . DP=0`（`FORMAT` 仅 `GT:AD`、`./.:0,20`）
外加覆盖位点 3-5 的 `ALT='.'` 行——该臂与未修的「覆盖位点物化」缺口纠缠，**本轮刻意不动**
（门禁 docstring 里写明了这个边界）。

### 修复与门禁

- 助手 `gatk_merged_record_depth()`（复刻 `calculateVCDepth`）+ `gatk_skips_regenotyping()`，
  在**两条遍历**的 compute 阶段**最前面**判定（必须在 `apply_gatk_output_allele_subset` 之前：
  GATK 在该路径上根本不会执行 `recordDeletions()`，故也不能登记删除状态）；
- 新严格门禁 `fastgatk-genotype-gvcf-depth-gate-gatk-oracle`（5 用例：3 个分歧用例 + 2 个对照，
  含 `INFO/DP=1` 边界与纯合参考对照）；
- 既有语料回归风险评估：注册夹具中**没有** `DP=0` 的记录，真实 HC gVCF 的 `INFO/DP` 恒 > 0
  （真实语料可达性见上一轮报告：165,527 条记录中 0 条命中该形状）。

### 仍然存在的**相邻**分歧（已测，未修，无门禁）

输入**没有** `INFO/DP` 而 `FORMAT/DP=20` 时：GATK 会按 `calculateVCDepth` 算出 20 并**发布** `INFO/DP=20`
（且 QD=4.63），native 既不发布 `INFO/DP` 也不写 QD。与深度前提条件同根（同一深度计算），
但方向相反（不是多输出而是少注释）；真实 HC gVCF 恒带 `INFO/DP`，故可达性同样接近 0。


## 第 51 轮：只含跨接删除的位点必须被拒绝（第 33 个已修 bug）

### 规则（源码逐行核对）

`GenotypingEngine.calculateGenotypes()` 在把已归属的 `*` 留在输出等位基因子集里之后，
**在登记任何删除状态之前**直接拒绝该位点：

```java
// return a null call if we aren't forcing site emission and the only alt allele is a spanning deletion
if (! emitAllActiveSites() && outputAlternativeAlleles.alleles.size() == 1
        && Allele.SPAN_DEL.equals(outputAlternativeAlleles.alleles.get(0))) {
    return null;                                   // :172-175
}
final List<Allele> outputAlleles = outputAlternativeAlleles.outputAlleles(vc.getReference());
recordDeletions(vc, outputAlleles);                // :178-179
```

`emitAllActiveSites()` 就是 `EMIT_ALL_ACTIVE_SITES` 遍历，在 GenotypeGVCFs 里等价于
`--include-non-variant-sites`。故规则**按模式分裂**：默认模式拒绝、dense 模式发射；
而 ORPHAN `*`（无覆盖删除）根本走不到这里——它已在 `:314` 作为 spurious spanning deletion 被剪掉。

### 实测（主会话亲自跑 pinned GATK 4.6.2.0）

| fixture | GATK | native（修前） |
| --- | --- | --- |
| 默认：`2 AAA A,<NON_REF>` + `3 A *,<NON_REF>`（`*` 被覆盖，纯合） | **仅 1 行**（位点 2） | 2 行（多出 `chr1 3 . A * 0 LowQual …`） |
| 默认：同上但 `*` 为杂合 | **仅 1 行** | 2 行 |
| dense：同上（`--include-non-variant-sites`） | 发射 `3 A *`（例外） | 发射（仅 QD 符号差异，见下） |
| dense：`3 A *,G,<NON_REF>`（有具体 ALT 存活，对照） | 一致 | 一致 |
| 默认：orphan `*`（无覆盖删除，对照） | 无记录 | 无记录（本来就对） |

### 修复与门禁

在 `apply_gatk_output_allele_subset()` 中、**任何 `upstream_deletions.record()` 之前**判定
`!include_non_variant_sites && output_alleles.size()==2 && output_alleles[1]=="*"` → 丢弃记录。
顺序是硬约束：GATK 的拒绝发生在 `:178-179` 之前，所以被拒绝的位点**不得**登记删除状态。

门禁 `fastgatk-genotype-gvcf-star-only-locus-gatk-oracle`（5 用例，双后端通过）。

### 同轮被**当场否证**的修复（记录在案）

同一轮我先尝试修 dense 下 `*`-only 行的 `QD=-0.00`，判据是「GATK 的 QD 分子是未取整的
`-10*log10PError`（`QualByDepth.java:78/:86`），而 native 拿四舍五入后的 QUAL 当分子」，
于是改了两处：QD 分子**与**浮点文本格式器（对 QD 保留负零）。

**注册门禁当场把它否证**：同一个 fixture 里 GATK 在位置 3 写 `QD=-0.00`、在位置 4 写 `QD=0.00`
（两行 QUAL 都渲染成 `0`）。也就是说符号来自 AF 计算器 **约 1e-16 的舍入方向**，不是呈现规则。
实测 native 该行的 `call_confidence = +4.82e-16`（p 略小于 1），GATK 落在另一侧（p 略大于 1，
Java 的 `log10PError` 因此为正 → 分子 `-4.8e-16`）。要逐字节复现这个符号必须复刻
GATK 的浮点运算顺序（内核级），Host 呈现层做不到。

处置：两处改动**全部回退**（`git diff` 现在只剩 `*` 拒绝那一段），该分歧记为
「已测量、机制已定位、未修」，证据留在未注册的
`verify_genotype_gvcf_dense_materialize_gatk_oracle.py` 的 REPORTED ONLY 用例里；
新门禁中这一行只比前 7 列（CHROM..FILTER）并在 `why` 里写明原因。

> 教训：**符号零这类 1e-16 级差异必须用「多于一处的位点」验证**，否则会把
> 「换一个例子又反了」的伪修复当成修复。这次是**注册门禁自己抓到的**——
> 反向 trim 门禁的位置 4 行从 `QD=0.00` 变成了 `QD=-0.00`。


## 第 52 轮：dense 模式跨位点记录物化（第 34 个已修 bug，本会话最大的一处结构性缺口）

### 规则

GATK 的位点遍历访问**每一条输入记录跨过的每一个参考坐标**，而不只是记录起点的坐标
（`VariantLocusWalker.java:150-176`），`--include-non-variant-sites` 随后把这些位点发布出来。
在这样被跨过的位点上，merger 把跨接事件的删除等位基因换成符号 `*`
（`ReferenceConfidenceVariantContextMerger.java:150-151`、`:222-245`），样本数据则是它的
reference-confidence 投影（`mergeRefConfidenceGenotypes()`、`:575-612`）。

native 的位点集合过去只由「记录起点 ∪ 纯参考块逐坐标展开」构成，因此**被跨接记录覆盖的位点整行缺失**。
这正是上一轮以前一直记录的「结构性、刻意未修」项，也是唯一有真实语料证据的缺口：
GATK 自带的 chr20 HaplotypeCaller gVCF（1291 条）里有 **139 个**这样的坐标。

### 实测（gate 的 4 个 gated 用例，pinned GATK 4.6.2.0）

| fixture（dense） | GATK | native（修前） |
| --- | --- | --- |
| `2 AAA A,<NON_REF>`（span 2-4） | 3 行（2/3/4） | **1 行**（缺 3、4） |
| 同上 + `5 A *,G,<NON_REF>` | 4 行 | 2 行（缺 3、4） |
| `2 AAA A` + `4 AA *,<NON_REF>` | 4 行 | 3 行（缺 3） |

### 实现

新助手 `materialize_spanning_loci()`，**只在 `--include-non-variant-sites` 下**、
在 `split_reference_blocks_at_variants()` 之后、排序/分组之前运行（聚合路径，即默认路径）：

- 只在「被某条记录跨过、且自身没有记录起点」的位置合成记录（GATK 保留起点记录：
  `GenotypeGVCFsEngine.java:339-354`，故 `starting-record-wins-dense` 对照不变红）；
- 仅当**源记录的基因型确实调用了那个删除等位基因**时才合成 `*` 行。hom-ref 的跨接记录贡献 NO_CALL，
  GATK 那时发布的是 REF-only 行——另一种形状，**仍未修**（`pruned-deletion`、`no-deletion-alt`
  两个对照用例正是在钉这一点：它们必须保持「没有任何 `*` 行」）；
- 样本数据复用既有的 allele-union 重映射（`remap_record_to_allele_union`）投影到
  `[新参考碱基, *]`：对单一删除等位基因的源记录这就是等位基因下标上的恒等映射
  （GT/AD/PL 逐项与 GATK 实测一致）；
- 新参考碱基取自 `-R` 索引（与既有分块路径同一助手），并把 `END` 清掉。

### 顺带收窄的 QD 符号零规则（与第 51 轮那次被否证的尝试不同）

合成的跨位点行有 **4/4** 实测为 `QD=-0.00`（QUAL 恒为 0），与 regenotyped 行的 1e-16 不确定性
（同一 fixture 内 `-0.00` 与 `0.00` 并存）**性质不同**：合成位点没有自己的输入记录，
其站点置信度是**恰好 0**。故新增 `Record::materialized_spanning_locus` 标志，
仅在「该行由本 pass 合成 且 QUAL==0」时把 QD 分子取为 `-0.0`，
并让浮点文本格式器只对 `QD` 保留负零（第 51 轮那次是**整体**放宽，被注册门禁否证并回退；
这次收窄到合成行，四个基因型门禁同时通过）。

### 门禁

`fastgatk-genotype-gvcf-dense-spanning-loci-gatk-oracle` —— 即此前**刻意未注册**的
`verify_genotype_gvcf_dense_materialize_gatk_oracle.py`：本轮起严格通过（0 violations）并注册。
其中原「REPORTED ONLY」的 `owned-star-only-locus-default-mode-refused` 用例
（默认模式 `*`-only 拒绝）已由第 51 轮修复，本轮**升格为 gated**；
只剩 `owned-star-only-locus-dense-negative-zero-qual`（regenotyped 行的 QD 符号零，1e-16 级）
仍标为 REPORTED ONLY、不参与退出码。

### 真实语料验证（chr20，GATK 自带 gVCF）

见 `fastgatk-native/evidence/2026-09-13-round52-dense/real-corpus-chr20-verification.md`。
窗口 `20:10000000-10003000`（dense，两侧同 `-L`）：

| 指标 | 结果 |
| --- | --- |
| GATK 有而 native 没有的位点 | **0**（修复前这一整类坐标整行缺失） |
| `*` 行逐字节不一致的位点 | **0**（本轮新增的行与 GATK 完全一致） |
| 其余不一致 | 69 个 REF-only 行 + 9 个 other，均为**既有**分歧（下） |
| native 多出的位点 | 20 个，全在窗口末端之后，`ALT='.'`（既有分歧） |

同一语料暴露的两处**既有**分歧（与本轮无关，已测量、未修）：

1. **dense 的参考块展开不遵守 `-L`**：窗口 `20:10050000-10051000`（1000 bp、无删除等位基因，
   本轮 pass 不触发）GATK 恰好输出 1001 行（每坐标一行），native 输出 4809 行
   （界内 4362 行 = 同一位点多行，界外另 447 行）。窗口 1 多出的 20 行是同一机制。
2. **块展开行的 REF 取错来源**：窗口 1 中 69 个位点两侧都是 REF-only 行但 REF 字母不同——
   **GATK 用该位点自身记录的 REF（`T`/`C`/`A`），native 用 FASTA 碱基**（该参考在这些位置确为 `N`）。
   例：`20:10000000` 的 gVCF 记录为 `T <NON_REF>`，GATK 输出 `T`，native 输出 `N`。

### 本轮未覆盖（下一步）

1. **流式路径**：pass 目前只接在聚合路径（默认路径）。`--stream-by-locus` 的等价接入点未做，
   故该选项下仍缺这些行（既有行为，未回归）。
2. **非删除跨接的 REF-only 形状**：`2 AAAA AACC,<NON_REF>` 这类无删除等位基因的跨接记录，
   GATK 在被覆盖位点发布 REF-only 行（`GT:AD ./. :0` / `QUAL 192.21` 等两种形状）；
   native 仍不产生这些行（两个对照用例只钉「不得出现 `*` 行」）。
3. **多 ALT 源记录**：本 pass 用「保留 REF + 被调用的那个删除等位基因」表达投影；
   源记录有多个 ALT（其中若干需变成 NO_CALL）的更一般情形未做。

## 第 53 轮：dense 记录供给必须裁剪到 `-L` 区间（第 35 个已修 bug）

第 52 轮的真实语料验证顺带暴露了这条既有分歧，本轮先用最小夹具钉住、再修。

### 规则与实测

GATK 只遍历**被请求的区间**，`--include-non-variant-sites` 的逐坐标展开因此在区间末端停止。
最小夹具：一条 `END=11` 的纯参考块记录 `2 A <NON_REF>` + `-L chr1:5-7`：

| | 行数 | 位点 |
| --- | --- | --- |
| GATK | 3 | 5, 6, 7（裁剪到区间） |
| native（修前） | **10** | 2..11（把整块展开） |

真实语料上同机制且量大得多：`-L 20:10050000-10051000`（1000 bp）时
native 输出 **4809** 行、位点范围 10046639-10051447（整个外层块），GATK 输出 1001 行、范围正好落在区间内。

### 修复

在聚合路径里，`split_reference_blocks_at_variants()` 与本轮的跨位点物化之后，
当 `--include-non-variant-sites` 且确实给了 `-L` 时，用既有的
`starts_in_regions()` 把记录供给裁剪到区间（`Record` 逐个判 rid/pos）。
同时把区间传给 `materialize_spanning_loci()`：**起点在区间外的源记录不再参与合成**——
因为那个位点根本不会被发射，GATK 也就不会登记它的删除，被覆盖坐标上的 `*` 会变成 orphan，
与遍历状态不一致。

### 真实语料复验（chr20，GATK 自带 gVCF）

| 窗口 | GATK 行/位点 | native 行/位点 | 仅 GATK 有 | 仅 native 有 |
| --- | --- | --- | --- | --- |
| `20:10000000-10003000` | 3001 / 3001 | **3001 / 3001** | 0 | 0 |
| `20:10050000-10051000` | 1001 / 1001 | **1001 / 1001** | 0 | 0 |

（修前分别是 3021 与 4809 行。）剩余差异只剩下面这一类。

### 同一验证暴露的、仍未修的下一个目标（已精确定位）

**块展开行的 REF 取错来源**：两个窗口里剩余差异分别是 82 与 3 个位点，全部是
REF-only 行、且**只有 REF 字母不同**——**GATK 用该位点自身记录的 REF，native 用 FASTA 碱基**。
例：`20:10000000` 的 gVCF 记录是 `T <NON_REF>`，而 `GRCh37.chr20.fa` 在该处是 `N`；
GATK 输出 `REF=T`，native 输出 `REF=N`。窗口 2 的 3 个例外同理（块起始坐标）。
即：**块的起始坐标应保留记录自身的 REF**（那个位点就是记录本身，不是合成行），
只有内部坐标才用 FASTA 碱基。这条已可单独立轮修。

## 第 54 轮：块起始坐标必须保留记录自身的 REF（第 36 个已修 bug）

第 53 轮真实语料验证留下的 82 / 3 个「只有 REF 字母不同」的位点，本轮修掉。

### 规则与实测

locus 的合并 REF 取自**在该坐标起始的记录**，只有在没有记录起始的坐标才回落到参考碱基
（merger 收到 `ref.getBase()`，`ReferenceConfidenceVariantContextMerger.merge()`）。
最小夹具（全 `A` 参考 + 记录 `2 C <NON_REF> END=5`，dense）：

| | POS 2 | POS 3/4/5 |
| --- | --- | --- |
| GATK | **C**（记录的 REF） | A（参考碱基） |
| native 修前 | A ✗ | A |
| native 修后 | **C** | A |

修复两处：解码阶段的逐坐标展开只在 `position != start` 时用 FASTA 碱基；
`split_reference_blocks_at_variants()` 中起点等于块起点的分段保留 `original.alleles.front()`。

门禁用例 `dense-block-start-keeps-the-record-ref`（修前确认失败）已加入
`fastgatk-genotype-gvcf-dense-spanning-loci-gatk-oracle`。

### 真实语料复验

| 窗口 | 差异位点（修前 → 修后） |
| --- | --- |
| `20:10000000-10003000` | 82 → **13**（9 个只有 QD 不同；4 个为既有的「变异行退化」） |
| `20:10050000-10051000` | 3 → **0**（该窗口现已逐字节完全一致） |

两窗口的「仅 GATK 有 / 仅 native 有」位点仍均为 0。

### 剩余两类（已定位，未修，列入下一步）

1. **QD 抖动顺序（9 个位点）**：等位基因相同、只有 `QD` 不同。`QualByDepth.fixTooHighQD()`
   在 `QD >= 35` 时把它映射为 `30 + nextGaussian()*3`，两边的**随机数消耗顺序/次数不同**
   （实测：`20:10001298` GATK `27.24` vs native `30.97`，而 `20:10000758` GATK 恰为 `30.97`）。
2. **变异行退化为 REF-only 行（4 个位点，既有）**：`20:10000758` 输入 `T A,<NON_REF>`（同型变异），
   GATK 输出变异行，native 输出 `T . ... MLEAC=.;MLEAF=.`（ALT 被剪掉、走了 REF-only 物化路径）。
   已确认**与第 54 轮改动无关**（修前/修后该行逐字节相同）。

## 第 55 轮：等位基因频率后验的**下溢**（第 37 个已修 bug，本会话最严重的一处丢调用）

### 症状与定位路径

真实语料验证（第 54 轮）留下 13 个差异位点，分两类：9 个只有 QD 不同、4 个 native 发 REF-only 行
而 GATK 发变异行。追查第 4 类时发现**更严重的事实**：native 在**非 dense** 模式下的整份 chr20 语料
只输出 234 行，而 GATK 输出 252 行——**18 个真实变异调用被整条丢掉**（占 GATK 输出位点的 7%）。

用单记录夹具复现（`20:10000758`，同型变异 `1/1`、AD `0,95,0`、PL `3867,286,0,…`），
再用 `FASTGATK_DEBUG_SUBSET=1` 插桩（新增的环境变量开关，默认不输出）看到关键证据：

```
[SUBSET] pos=10000757 nalleles=2 pruned=1 output=[T,] absent=[0,0,]      <- 缺陷
[SUBSET] pos=10000116 nalleles=2 pruned=0 output=[C,T,] absent=[-126.138,-92.664,]   <- 正常
```

`absent=[0,0]` 表示核函数**完全没有产出信息**：`calculate_allele_frequency_kokkos()`
在**线性空间**归一化后验

```cpp
probability = pow(10, term - maximum) / denominator
```

一个高置信同型变异样本的 `P(纯合参考)` 远小于最小的正 double（PL 差动辄数千），
于是该式**下溢为 0**，`log10(p0)` 变成 `-inf`，而宿主把它当作「该样本没有似然值」的哨兵跳过；
所有样本都被跳过后 `log10_p_allele_absent` 保持初值 `0.0`，标准置信度子集据此判定 ALT 不可信而剪掉它，
默认遍历下整条记录消失。

18 个丢失位点的形态完全一致（GT `1/1`/`1|1`、AD 参考深度 0、PL 跨度 ~2600-4000），
而 79 个存活的同型变异位点 PL 跨度较小 —— 与「下溢」判据吻合。

### 修复（在内核内完成，Host/kernel 边界未变）

`fastgatk-kernels/src/genotype.cpp`：把最终后验的两项累积改为 **log 空间（log-sum-exp）**，
不再物化任何线性概率；`device_sample_p0` 与 `device_sample_absent` 直接写 log10 值，
故宿主原有的 `std::min(0.0, …)` 求和与哨兵判断全部照旧。

### 效果（同一份 chr20 语料，主会话实测）

| 场景 | 修前 | 修后 |
| --- | --- | --- |
| 非 dense 全语料行数 | 234（GATK 252） | **252 = GATK**，仅 GATK 有/仅 native 有均为 **0** |
| dense `20:10000000-10003000` 差异位点 | 13 | **1** |
| dense `20:10050000-10051000` 差异位点 | 0 | **0** |
| 单记录夹具 `20:10000758` | 0 行 | **与 GATK 逐字节一致**（含 `QD=25.36`） |

**一个根因、两类症状**：QD 抖动的错位是「少了一次调用」的连带结果——修复后
`QualByDepth.fixTooHighQD()` 的 `nextGaussian()` 抽取次数与顺序重新对齐，9 个 QD 差异位点随之消失。

### 门禁

新严格门禁 `fastgatk-genotype-gvcf-confident-hom-alt-gatk-oracle`（3 用例：高置信同型变异必须发射且逐字节一致、
杂合对照、同型参考下 ALT 仍必须被剪掉且不发射），双后端 ctest 通过。

### 语料上仍剩的 3 个差异位点（已定位，未修，下一轮）

| 位点 | 差异 |
| --- | --- |
| `20:10002458`（`G>GTTT`）、`20:10024300`（`C>CTT`） | FORMAT/PL **整体平移**：GATK `2090,167,0` vs native `2172,249,82`（差恒为 +82）；另一处 +45 |
| `20:10068160`（多等位 `*` 记录） | 只有 `QD`：GATK `9.33` vs native `9.34` |

## 第 56 轮：子集投影后的 PL 未归一化（第 38 个已修 bug）

第 55 轮之后语料上只剩 3 个差异位点，其中 2 个是 FORMAT/PL「整体平移」，本轮修掉。

### 机制（用真实 PL 行逐值核对）

源记录是多等位的，输出等位基因子集把 ALT 剪到 1 个之后，**投影出的 PL 行保留了源行的偏移**；
htsjdk 重建子集基因型时会经过 `GenotypeLikelihoods`，它把每个样本的 PL 行**归一到最小值为 0**。

`20:10002458`（源 `G GT,GTTT,GTTTTT,<NON_REF>`，PL `2172,655,493,249,0,82,352,33,9,214,876,512,167,254,730`，
最优基因型 (1,2) 在子集中不存在）：

| | PL |
| --- | --- |
| GATK | `2090,167,0` |
| native（修前） | `2172,249,82`（= 源行在投影坐标 {(0,0),(0,2),(2,2)} 上的取值） |
| 差值 | 恒为 **+82**（= 投影行的最小值） |

`20:10024300` 同形，恒差 **+45**。

### 修复

新助手 `normalize_sample_pl()`：按样本把 PL 行的有效项减去其最小值（不触碰 missing / vector_end），
在 `apply_gatk_output_allele_subset()` 的**两条返回路径**上调用（`pruned == 0` 路径与投影后的
发布点）。它与子矩阵投影可交换，故放在投影之后即可与 htsjdk 一致；AF/注释在投影之前完成，
均匀平移对后验与 GQ 差值不变量，因此不影响 QUAL/GQ。

### 效果（同一 chr20 语料，主会话实测）

| 场景 | 修前差异位点 | 修后 |
| --- | --- | --- |
| 非 dense 全语料 | 3 | **1** |
| dense `20:10000000-10003000` | 1 | **0** |

五个基因型门禁（reverse-trim / depth-gate / star-only / dense-spanning / confident-hom-alt）全部严格通过；
其中 `confident-hom-alt` 门禁新增用例 `multi-allelic-projected-pl-is-normalized`
（夹具直接取自语料的真实 PL 行，修复前 native 会输出 `2172,249,82`）。

### 全语料仅剩的 1 个差异位点（已定位，未修）

`20:10068160`（多等位 `G,*` 记录）：只有 `QD` 不同，GATK `9.33` vs native `9.34`。
QUAL 56.01 与深度 6 在两侧相同，故这是 2 位小数渲染上的**浮点 tie-break** 差异
（与既有 rank-sum 的 epsilon 微调同类），需按 Java `Formatter` 对精确 double 的行为单独核对。

## 第 57 轮：QD 的分子必须用**未取整的 double**（第 39 个已修 bug）

第 56 轮之后语料上只剩 1 个差异位点（`20:10068160`，QD `9.33` vs `9.34`），本轮修掉 ——
**至此 GATK 自带 chr20 语料的三种遍历全部逐字节一致**。

### 机制（插桩实测，非推理）

`[QDDBG] pos=10068159 qual_float=56.009998321533203 cc=56.007546683221719 depth=6
qd_from_float=9.3349997202555333 qd_from_cc=9.3345911138702871`

- GATK 的 `QualByDepth` 用的是**未取整的 double**：`qual = -10.0 * vc.getLog10PError()`
  （`QualByDepth.java:78`），再 `QD = qual / depth`（`:86`）；
- native 用的是**已取整为 2 位小数并降为 float32 的 QUAL**：9.335 在 float32/2 位小数往返后
  落在十进制 9.335 的**上侧**，于是渲染成 `9.34`；GATK 由 9.3345911 渲染成 `9.33`。
- 即：这不是格式化器的问题，而是**分子取错来源**。native 内部早就有那个未取整的 double
  （`record.call_confidence` = `-10*log10PError`）。

### 修复

QD 分子改为优先使用 `record.call_confidence`（未取整 double），**但排除与 0 无法区分的情形**
（`|cc| <= 1e-9`）：那里的符号零来自 AF 计算器 ~1e-16 的舍入方向，GATK 自身都不自洽
（第 51/53 轮实测），保持用取整后的 QUAL 以免回归；合成跨位点行的 `-0.0` 规则（第 53 轮）保留。

### 效果（同一 chr20 语料，主会话实测）

| 场景 | 修前差异 | 修后 |
| --- | --- | --- |
| 非 dense 全语料（252 行） | 1 | **0**（逐字节一致） |
| dense `20:10000000-10003000`（3001 行） | 0 | **0** |
| dense `20:10050000-10051000`（1001 行） | 0 | **0** |

三个场景的「仅 GATK 有 / 仅 native 有」位点数也均为 0；五个基因型门禁全部严格通过。

> 说明：这一处的 tie-break 目前**没有合成夹具门禁**（该位点的 REF 是 `GTATATATATATGTA`，
> 需要专门的重复序列参考），它由上面这份语料级证据锁定；下一轮可补一个专用夹具。

## 收尾基线（第 42 轮起持续更新，主会话亲自运行）

**最新基线（第 57 轮，主会话亲自运行）：commit `9a2e3a5`（QD 分子用未取整 double；两后端均已重建）上
OpenMP 309/309（1473.5s）、Serial 309/309（1480.8s），零陈旧告警，
运行器默认强制 `FASTGATK_REQUIRE_GATK_ORACLE=1`。**

证据块（可直接复核）：`.diag/regression/20260913-052529/summary.txt`

> 流程教训（已在本轮踩到）：改动源码后**必须两个后端都重建**再跑全量。
> 本轮第一次跑全量时只重建了 omp，serial 用陈旧二进制跑 dense 门禁而失败
> （失败信息 = 新用例 `record count differs: GATK=3 native=10`，正是修前行为）。
> 该次运行作废（`.diag/regression/20260913-024535/`），重建 serial 后重跑才是本基线。

（更早的基线，均由主会话亲自测得，非委派方代跑：
第 56 轮 `5302ea8` 上 309/309（1503.4s / 1515.2s，PL 归一化）；
第 55 轮：AF 核日志域累积修复上 308/308（1490.3s / 1497.2s）；
第 54 轮 `0e19bbe` 上 308/308（1482.3s / 1458.1s）；
第 52 轮 `3389e98` 上 308/308（1504.4s / 1496.3s）；
第 51 轮 `a0a30b7` 上 307/307（1458.9s / 1422.5s）；
第 49 轮 `640678b` 上 305/305（1566.2s / 1531.5s）；
`0862251` 上 304/304（1622.6s / 1444.9s，工作树干净）；
`c971cf3` 上 303/303；`b3cd293` 上 302/302。
第 49/50 轮之所以必须重跑，都是因为生产代码改了 `genotype_gvcf_tool.cpp`。）

这条基线的意义：此前数轮的全量结果由委派方运行、我只做了 md5/时序核对；
自 `c971cf3` 起补上了「最终提交树上由主会话亲自测得」的那一步，因此
**「308/308 在强制 oracle 存在下成立」这一宣称有同源证据。**

配套的可信度条件（均已在本会话建立）：
1. `run_regression.sh` 默认要求 GATK oracle 在场（缺失即响亮失败），
   并只与**最新产物**比较陈旧性（消除假告警）；
2. 176 个脚本经 `oracle_guard.py` 改为 fail-closed（原可静默降级为「与自身比较」）；
3. 本会话新注册 28 道严格 GATK 门禁（另有 1 道 `-diagnostic` 门禁以 exit 0 记录差异、
   刻意不使全量变红），覆盖本会话 39 个修复中的关键行为；
4. 5 道刻意未注册（`verify_hc_forced_alleles_emission_gate_oracle.py`、
   `verify_reblock_gvcf_triploid_gatk_oracle.py`、
   `verify_genotype_gvcf_dense_materialize_gatk_oracle.py` 等），因其**按设计必须失败**——
   它们是尚未修复分歧的活证据，不应被注册成会永久变红的测试。

**仍未达成 1:1**（按剩余体量排序，均已在正文各节记录并可复现）：
1. **该 tie-break 尚无合成夹具门禁**（第 57 轮）：`20:10068160` 的 REF 是 `GTATATATATATGTA`，
   需要专门的重复序列参考才能做成夹具；目前由语料级证据锁定；
2. 第 52 轮记下的 **dense 跨位点物化的剩余两半**：流式 `--stream-by-locus` 的等价接入点、
   非删除跨接记录的 REF-only 形状（`GT:AD ./. :0` / `QUAL 192.21` 两种）；
3. 多输入 header（对齐 `CombineGVCFs` 语义：样本列字典序、后续输入声明并集）；
2. dense 跨位点物化的**剩余两半**（第 52 轮已修默认路径）：流式 `--stream-by-locus` 的等价接入点、
   以及非删除跨接记录的 REF-only 形状（`GT:AD ./. :0` / `QUAL 192.21` 两种）；
2. 多输入 header：**对齐目标已更正**——GATK 拒绝多个 `-V`，故应对齐 `CombineGVCFs` 合并语义；
   实测分歧为样本列顺序（字典序 vs 输入序）与后续输入声明（INFO/FORMAT/FILTER/ALT）丢失；
3. ReblockGVCF case B（修法已验证但会使现有 trim/gap/NON_REF-AD 断言块不可满足，
   需先重做其约 25 个派生值）；
4. 退出码类分歧（空等位基因 `--alleles` GATK exit 3 vs 本实现 exit 0；
   SAM 文本路径受 htslib `sam_read1_sam` 折叠解析失败为 -1 所限，不可修）；
5. 第 50 轮记下的同根反向残差：输入无 `INFO/DP` 时 GATK 发布算出的 `INFO/DP` 与 QD，native 两者都不写；
6. 第 51 轮记下的符号零残差：dense 下 `*`-only 行的 `QD` 符号取决于 AF 计算器 ~1e-16 的舍入方向，
   GATK 自身在同一 fixture 内不一致（位置 3 为 `-0.00`、位置 4 为 `0.00`）；已测量、机制已定位，
   需在内核层复刻 GATK 的浮点运算顺序才可能对齐，Host 呈现层已证明修不好（伪修复被注册门禁否证）。

此外还有一批**已测量但尚未设门禁**的残差（多 contig 顺序、非 ASCII Description、
NaN 补集、「同一位点两条记录」、跨位点替换的 REF-only 物化、默认模式 `*`-only 行、
缺样本 `'./.'` vs `'./.:.:.:.,.,.:.'`、QD `-0.00` 渲染、
具体记录的 `INFO/END` 透传），以及
Mutect2 6 处「推理上脆弱但从未复现」的所有者位点（44/44 行与 GATK 一致）。
