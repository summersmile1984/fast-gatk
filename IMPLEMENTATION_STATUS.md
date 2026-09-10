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
| 全量回归 | **280/280 通过**（1002s） | **280/280 通过**（1000s） | 双后端并行执行 |
| HC/Mutect2 子集 | **72/72 通过**（261s） | **72/72 通过**（268s） | `omp-hc-mutect2.log` / `serial-hc-mutect2.log` |
| `verify_hc_alleles_gatk_oracle.py` | 通过（86s） | 通过（112s） | 原「待本轮复跑」项已结清 |
| `verify_hc_complex_multiallelic_oracle.py` | 通过（45s） | 通过（33s） | 四倍体 / max-ALT / max-genotype-count |
| `verify_hc_chr20_max_mnp_gvcf_gatk_oracle.py` | 通过（47s） | 通过（53s） | MNP + gVCF `*` + max-ALT |

上述结论对应的二进制：`fastgatk-native/build/fastgatk-hc-call`（17:03）与
`fastgatk-native/build-serial/fastgatk-hc-call`（17:08）。两者都比当时**全部**源文件新
（最新源改动为 `fastgatk-native/src/calling_pipeline.cpp` 17:03），因此本轮 Host 注入与
gVCF 改动确实被这轮回归覆盖。

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
| Track A | gVCF 流式/分区：重叠 indel、符号等位基因 × max-ALT、read-overlap margin、重分块组合 | `.diag/track-a-gvcf-streaming-findings.md` |
| Track B | `--alleles` / GenotypeGivenAlleles 复杂注入：重复序列 indel、跨 AssemblyRegion、相邻/重叠 feature、排名并列 | `.diag/track-b-alleles-findings.md` |
| Track C | Mutect2 独立复核：共享 calling pipeline 改动后是否回归、正常样本重放、联合 AssemblyResultSet、somatic gVCF、TLOD 精度 | `.diag/track-c-mutect2-recheck-findings.md` |

并行原则（本仓库的硬约束）：**oracle/夹具/审计可并行；生产代码改动必须单线串行**，
因为核心路径的全部修复都落在同一批文件（`calling_pipeline.cpp`、`kmer_graph.cpp`、
`hc_call.cpp`、`mutect2_tool.cpp`），并行编辑必然互相覆盖。所有 oracle 都使用
`tempfile.TemporaryDirectory`，因此测试进程可安全并发。

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
- **静默 skip 风险**：117 个 `verify_*.py` 在缺少 GATK/oracle 输入时会打印
  `{"status":"skip"}` 并返回 0；注册的 267 个测试中只有 134 个设置了
  `FASTGATK_REQUIRE_GATK_ORACLE=1`。因此「全绿」可能包含静默跳过。已启动一次全量
  verbose 复跑统计实际 skip 数，结果与处置追加到本节。

## 下一步执行顺序

1. ~~结束并确认 Serial 构建，完成 Serial oracle。~~ 已完成，见「最近验证状态」。
2. ~~跑完整 HC/Mutect2 回归。~~ 已完成；并已取得全量双后端 280/280 证据。
3. 汇总三路审计的差异清单，按影响核心输出正确性的程度排序。
4. **单线串行**落地修复：一次一个改动 → 双后端增量构建 → 定向 oracle + 72 子集 →
   全量回归 → commit（有 git 后「改动前/后」才可验证）。
5. 补齐静默 skip 审计与 tool audit 口径，使百分比与证据同源。
6. 只有当每个核心选项都有匹配范围的当前证据时，才讨论「1:1 完成」。

