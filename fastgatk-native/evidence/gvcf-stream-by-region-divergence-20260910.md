# `--stream-by-region` 与 GATK / 非流式输出不一致（2026-09-10）

状态：**已复现，P1（流式路径正确性问题）**。发现于第 1 轮 Track A 审计，
第 2 轮完成定位（**修复尚未开始，本报告不含代码改动**）。

## 结论摘要

`fastgatk-hc-call --stream-by-region N`（当 tile 小于区间长度时）的输出与
**同一条命令去掉 `--stream-by-region`** 的输出不一致，也与 pinned GATK 4.6.2.0 不一致。
第 2 轮已验证出三类不一致：

1. **参考块粒度改变（结构性）**：非流式路径与 GATK 都把 20:10020230–10020428 拆成
   29 个小 `<NON_REF>` 块；流式路径把它合并成 **1 个** 块（`POS=10020230 END=10020428`）。
   在 tile=100/500 时，流式输出比非流式输出**少 28 条记录**。
2. **注解读数错误（数值性）**：在 20:10020680 这一行，流式路径给出的
   `RAW_MQandDP`（`28800,8 → 97200,27`）与 `SB`（`0,0,3,3 → 0,0,0,0`）与非流式、与 GATK 都不同。
3. **凭空产生 phasing（语义性，第 2 轮新发现）**：在 20:10020228/10020229，非流式与 GATK 输出
   未定相的 `GT=0/1`，流式输出为已定相 `0|1`/`1|0` 并新增 `PGT/PID/PS`。

判定基准：非流式路径与该区间上的 GATK **全部 48 条数据行逐字段相同**，所以非流式是对的、
流式是偏离方。第 2 轮另已确认 **只有在 tile 覆盖整个输入区间时才逐字节相同**，
且加大 `--assembly-region-padding` **不是修复**（见「已验证」一节）。

## 重要更正（第 3 轮）：D2 不是流式缺陷，而是 `-L` 窗口依赖缺陷

本报告先前的定位框架是「非流式正确、流式偏离」。**该框架对 D2 是错的**，由独立审计发现、
并经我复核确认：

```bash
# 完全不加 --stream-by-region
fastgatk-hc-call -R fixtures/chr20/ref20mnp.fasta -I fixtures/chr20/mnp.bam \
  -L 20:10020381-10020710 --emit-ref-confidence GVCF --max-mnp-distance 1 --threads 1 \
  --add-output-vcf-command-line false -O n.g.vcf
```

在该窗口下，**没有任何流式**，native 仍然给出错误的
`RAW_MQandDP=97200,27` / `SB=0,0,0,0`。关键对照（同一夹具，两条 `-L`）：

| `-L` 窗口 | GATK | native | GATK vs native 数据行差异 |
| --- | --- | --- | --- |
| `20:10019901-10020710` | `28800,8` / `0,0,3,3`（48 条） | `28800,8` / `0,0,3,3`（48 条） | **0 行（逐字段相同）** |
| `20:10020381-10020710` | `28800,8` / `0,0,3,3`（13 条） | `97200,27` / `0,0,0,0`（13 条） | **2 行** |

即：**GATK 的该注释与 `-L` 窗口起点无关，native 的却有关**。因此

- D2 的真实性质是 **`-L` 窗口依赖**（interval-dependent annotation evidence），
  而非流式路径特有的错误；流式只是通过选择 tile 窗口**触发**了它。
- 这是一个**比流式问题更广的 GATK parity 缺陷**：普通非流式 HaplotypeCaller 只要用不同的
  `-L` 窗口（例如 scatter 分区）就会复现。本报告原标题中的「与 GATK 不一致」对 D2 而言，
  根因不在 streaming。
- 我先前「非流式是正确的一方」这一判断**只对 `20:10019901-10020710` 这一个窗口成立**，
  不能外推。

D1（凭空 phasing）与 D3（参考块粒度）的流式相关性见
`fastgatk-native/evidence/2026-09-10-parallel-audit/track-d-streaming-rootcause.md`：
该审计把 D1 钉到 `gvcf()` 的 phasing 块（`hc_call.cpp:4085-4195`），
把 D3 钉到 tile 内部的逐 locus RCM 取值（tile=500 在 10020230 算出
`depth=8 informative=0 gq=0 pl=0,0,0`，而非流式为 `depth=12 informative=12 gq=36 pl=0,36,494`）。

## 第 3 轮：D2 的 trace 级证据（仍未到行级根因）

用已有的位置门控 `FASTGATK_DEBUG_ANNOTATION_POSITION=<0-based pos>` 对同一夹具跑两个非流式窗口，
得到注释边界（`calculate_variant_annotations`）的逐证据 trace。

**跨接删除位点的锚点（0-based 10020678 / 1-based 10020679）**：

| 窗口 | 该候选的注释 pass | ref_f/ref_r/alt_f/alt_r | rank sums |
| --- | --- | --- | --- |
| `20:10019901-10020710` | candidate=4，identity 映射（source 100..105 → physical 100..105） | `3/1/0/2` | `-6.96e-17 / 0.623 / -0.253` |
| `20:10020381-10020710` | candidate=4，**先出现一次全零 pass** | `0/0/0/0` | `nan / nan / nan` |
| `20:10020381-10020710` | 之后 candidate=0，identity 映射（source 48..53 → physical 48..53） | `3/1/0/2` | 同左 |

即：**在偏移窗口下，跨接删除记录拿到了一次「零证据」的注释**（链计数全 0、rank sum 为 nan），
这与输出的 `SB=0,0,0,0` 直接吻合；而同一 implement 在原始窗口下拿到的是真实证据。

**同一窗口下另一个候选（0-based 10020679 / 1-based 10020680）**：原始窗口的 render-time pass
走 **外部 ordinal 映射**（`likelihood_source_records_are_external = true`，source 100..105 →
physical 235,236,237,238,239,**241**，注意 240 被跳过），并把原本判为 REF 的读改判为 ALT
（`alt_f/alt_r = 3/3`，与输出的 `SB=0,0,3,3` 吻合）；偏移窗口下该候选的两次 pass 均为全零。

**已排除**：`likelihood_source_records_are_external` 并非「从未传 true」的死参数 ——
`calculate_output_variant_annotations`（`calling_pipeline.cpp:13816`）明确按第 9 个位置参数传 `true`，
所以 render-time 重算确实走外部 ordinal 解析。

**仍未证明**：全零 pass 是「in-run 期间算出并被沿用」还是「render-time 重算时证据查找失败」。

## 第 4 轮：机制已定位到具体布尔量（决定性证据）

按上一节写明的实验，临时在 `calling_pipeline.cpp` 的 `debug_annotations` 块内加一行
`[FASTGATK_ANNOTATION_GATES]` 打印（**已还原，仓库未留改动**），重建 OpenMP HC 后对两个窗口取 trace。
对该分歧行（0-based 10020679 / 1-based 10020680，candidate=5）：

| 字段 | `20:10019901-10020710` | `20:10020381-10020710` |
| --- | --- | --- |
| `have_context_mapping_evidence` | **1** | **0** |
| `mapping_count` / `mapping_square_sum` | 8 / 28800 | **27 / 97200** |
| `strand_total` | 4（in-run）/ 6（render） | **0 / 0** |
| `reads_records` / `likelihood_source_count` | 113 / 113（in-run） | 116 / 116（in-run） |
| 输出 `RAW_MQandDP` / `SB` | `28800,8` / `0,0,3,3` ✓ 同 GATK | `97200,27` / `0,0,0,0` ✗ |

**机制**：`calculate_variant_annotations`（`calling_pipeline.cpp:4293-4305`）里 MQ 证据的门是

```cpp
const bool mapping_evidence = have_context_mapping_evidence
    ? context_mapping_evidence[record] != 0
    : (read_overlaps_annotation_interval(...) && annotation_read_survives);
```

- 原始窗口：`have_context_mapping_evidence = 1` → 走 context 分支 → 计数 8（= GATK）。
- 偏移窗口：`have_context_mapping_evidence = 0` → 落到几何回退分支 →
  `read_overlaps_annotation_interval` 放行了 **27** 条读（应 8 条），
  同时 `strand_total` 归零（`retained_for_hc_allele_annotations` 那条链没有产出任何 BestAllele 证据）。

**关键旁证**：同一位置、同一窗口下的 `candidate=1` 却有
`have_context_mapping_evidence = 1`、`mapping_count = 8`（正确）。所以缺失是**逐候选的**，
不是全局的 —— context 映射证据在偏移窗口下没有为跨接删除候选生成/保留。

### 已证伪

- **几何谓词带窗口依赖**：`read_overlaps_annotation_interval` 与
  `read_overlaps_hc_genotyping_interval`（`calling_pipeline.cpp:6694` / `6725`）都只依赖
  读的 `position`/`reference_end` 与候选的 `position`/REF 长度/margin，**不含任何窗口参数**。
  窗口依赖来自 `have_context_mapping_evidence` 这个逐候选布尔量，而非几何判定本身。

### 已进一步定位：哨兵值 `context_ordinal == UINT32_MAX`

`have_context_mapping_evidence` 的置位点在 `calling_pipeline.cpp:4152`，位于
`for (const auto index : likelihood_indices)` 循环**内部**；而该循环在
`calling_pipeline.cpp:4129` 会对哨兵值直接 `continue`：

```cpp
const auto context_ordinal =
    likelihood_result->likelihood_candidate_read_context_ordinals[index];
if (context_ordinal == std::numeric_limits<std::uint32_t>::max()) continue;   // 4129
...
have_context_mapping_evidence = true;                                          // 4152
```

由于**同一窗口、同一位置**的 `candidate=1` 有 `have_context_mapping_evidence = 1`，
说明外层的 `if`（4112-4114，要求 `likelihood_candidate_read_context_ordinals` 与
`likelihood_candidate_read_realignments` 均非空）是通过的。因此偏移窗口下
`candidate=5` 的情形只能是：**该候选的
`likelihood_candidate_read_context_ordinals[candidate_index]` 是哨兵 `UINT32_MAX`**，
即 PairHMM 证据收集阶段没有为这个跨接删除候选记录读 context。

于是链路是：**PairHMM 阶段未给该候选分配 read context（哨兵值）
→ `have_context_mapping_evidence = 0` → MQ 证据门落到几何回退分支
（`read_overlaps_annotation_interval`）→ 放行 27 条读（应 8 条）
→ `RAW_MQandDP = 97200,27`；同时 strand 证据归零 → `SB = 0,0,0,0`。**

### 第 5 轮：语义已定案 + 验收判据（GATK 侧已独立复核）

把「native 该用哪个读集合」升级为已定案的语义问题，并给出可机器判定的验收判据。

**GATK 的规则**（引自 pinned Java 源码，详见
`fastgatk-native/evidence/2026-09-10-parallel-audit/gatk-retain-evidence-semantics.md`）：
`HaplotypeCallerGenotypingEngine.java:194/197` 的 `retainEvidence(target::overlaps)`，
其中 `target = new SimpleInterval(mergedVC).expandWithinContig(informativeReadOverlapMargin=2)`
→ 本题为 `[10020678, 10020683]`；作用对象是
`realignReadsToTheirBestHaplotype/changeEvidence` **之后**的读集合
（`HaplotypeCallerEngine.java:959-966`），且同一集合被复用于注释
（`HaplotypeCallerGenotypingEngine.java:598-601`）。

**我独立复核了关键经验事实**（用 pinned jar 的 `--debug-genotyper-output`，两个窗口各跑一次）：

```text
Event at: [VC HC0 @ 20:10020680-10020681 Q. of type=MIXED alleles=[CA*, *, AT] attr={} GT=[] filters= with 8 reads and 20 disqualified
```

**两个窗口输出完全相同（同样 8 条 retained / 20 条 disqualified）**，且
`alleles=[CA*, *, AT]` 证实该记录含 `*`。因此：

- GATK 的 retained 集合**与 `-L` 窗口起点无关**，正确值恒为 **8**。
- native 在窗口 A 的 8 是**命中 context 分支**得来的正确值；窗口 B 的 27 是
  **几何回退分支**引入的错误值。结论：**(a) context 映射population 才是 GATK 的模型，
  (b) 几何回退是缺陷。**

**验收判据（可机器判定）**：该行在两个窗口下都必须是
`RAW_MQandDP=28800,8` 且 `SB=0,0,3,3`；同时窗口 A 的 48 行必须与 GATK 保持逐字节相同。
另注意 Track D 报告指出窗口 `10020421/10020431` 还存在第三种错误状态（`SB=1,2,6,16`），
修复后应一并复核。

**推荐修法（未落地）**：根缺陷是哨兵值本身 —— 跨接删除候选从未进入 `candidate_reads`，
故 `candidate_read_context_ordinals` 保持 `missing_context`（`calling_pipeline.cpp:10548`）。
推荐在 `calling_pipeline.cpp:9152-9153` 构造 `pairhmm_event_candidates` 时，
为已装配候选补上**同 locus 的符号兄弟等位基因**（同 tid/position/reference），
使其获得 context ordinal。该轨道论证 impact 小（`request_indices` 按 (read_id, haplotype)
去重 @9315-9322；`marginalization_row_ids` 按 (group_ordinal, locus_id, source_record) 键控
@11519-11532），但**未经编译与测量**。

**其余四个候选修法已被逐一否决**（同一报告）：
(i) 把回退换用 `read_overlaps_hc_genotyping_interval` —— 只是把源 CIGAR 区间平移 1bp，
仍会放行约 27 条，两个窗口都修不好；
(ii) 哨兵时令 `mapping_count=0` —— 能做到窗口无关但**值错**（GATK 是 8）且治不了 SB；
(iv) 在 `hc_call.cpp:4530` 只保留非空重算 —— **是空操作**，因为窗口 B 的重算结果并不空
（它就是 `97200,27`）；
(iii) 上游配对修复是正确方向，推荐修法即其收敛形式。

### 第 6 轮：根因更正 —— 未配对的是**重复候选**，不是 `*` 等位基因

上一轮引用的推荐修法（「为已装配候选补上同 locus 的**符号兄弟等位基因**」）建立在
「未配对的候选是 `*`（跨接删除）」这一假设上。本轮用临时插桩打印**候选 → 等位基因 →
context ordinal** 的映射（**已还原，仓库无改动**），该假设被直接推翻：

```text
-L 20:10019901-10020710
  [FASTGATK_ANNOTATION_CANDIDATES] index=5 ref=CA alt=AT context_ordinal=0

-L 20:10020381-10020710
  [FASTGATK_ANNOTATION_CANDIDATES] index=1 ref=CA alt=AT context_ordinal=0
  [FASTGATK_ANNOTATION_CANDIDATES] index=5 ref=CA alt=AT context_ordinal=4294967295 SENTINEL
```

即：

- **窗口 A** 在该位置只有一个候选：`index=5`，`CA→AT`，ordinal **有效**（0）。
- **窗口 B** 在该位置有**两个完全相同的候选**：`index=1`（`CA→AT`，ordinal 0，已配对）与
  `index=5`（`CA→AT`，ordinal **SENTINEL**，未配对）。

因此真正的根因是：**窗口 B 产生了同一 `(tid, position, REF, ALT)` 的重复候选**，
其中未配对的那一份（index=5）被渲染/注释采用；而窗口 A 没有重复，故一切正确。
`*` 等位基因与 `have_context_mapping_evidence` 无关。

**关键约束**：早先的 gates 打印显示该调用的 `grouped_candidates=1`，
即**该候选所在的输出组里只有它自己** —— 已配对的孪生候选（index=1）**不在组内**。
所以任何「在 `likelihood_indices` 内寻找兄弟」的写法都不会生效，
必须跨 `result.candidates` 查找同一 `(tid, position, REF, ALT)` 的已配对孪生。

**修复空间的重新评估**：

1. 让渲染采用已配对的孪生（而不是未配对的重复）——最接近 GATK 语义
   （GATK 对同一 VC 只用一套读集合），但需要改输出组的选择/去重逻辑。
2. 仅复用孪生的 context ordinal —— 能修 `RAW_MQandDP`；但**治不了 `SB`**：
   strand 路径读的是**该候选自己**的 likelihood 行，而重复候选那几行是 `-inf`
   （早先实测：窗口 B 下 index=5 的 `strand_total=0`，index=1 的是 4）。
   因此这只是部分修复，输出行仍然错。
3. 上游让重复候选根本不产生 —— 需先定位重复是在哪一步（合并/汇总）产生的。

**结论**：上一轮的推荐修法**作废**（靶子错了）。下一步应先定位重复候选的产生点，
再决定是「阻止重复」还是「渲染时优先已配对孪生」。验收判据不变：
两个窗口下该行都须为 `RAW_MQandDP=28800,8` 且 `SB=0,0,3,3`，且窗口 A 的 48 行保持逐字节相同。

**重复候选的可疑产生点（未证实，需先验证）**：`calling_pipeline.cpp:13738` 是一处
**无去重的候选追加**（partition merge）：

```cpp
merged.candidates.insert(merged.candidates.end(),
    part.candidates.begin(), part.candidates.end());
```

### 第 7 轮：产生点**已证实**就是这处无去重追加

上一节把它列为「可疑」，本轮用临时插桩直接判定（**已还原，仓库无改动**）。在该 append
之前遍历 `part.candidates × merged.candidates`，对 `(tid, position, reference,
alternate, reference_allele, alternate_allele)` 全同者打点：

| 窗口 | 该窗口内的重复位置 |
| --- | --- |
| `20:10019901-10020710` | `10020227`、`10020228`、`10020428`、`10020430`、`10020433`、`10020437` |
| `20:10020381-10020710` | **`10020678`（AC→TA）**、**`10020679`（CA→AT）** |

关键点：

- 窗口 B 的重复位置 **`10020679`（CA→AT）正是那个分歧行**
  （1-based `POS=10020680`，ALT `AT,*,<NON_REF>`）。
- 窗口 A 的重复发生在**其它**位置（`10020227`/`10020228`/`10020428`/…），
  **不含** `10020679` —— 这正是窗口 A 该行正确、窗口 B 该行错误的原因。

因此：**重复候选确实由 `calling_pipeline.cpp:13738` 这处 append 产生**，
它只在该位置产生的那一份未进入 PairHMM 配对（哨兵），并被渲染采用。

### 修复约束（关键，避免下一轮踩坑）

该函数是**成组按位置追加**的，紧随其后还有多个**与候选一一对应的并行数组**：

```cpp
merged.reference_confidence_regions.insert(...part.reference_confidence_regions...);
merged.likelihoods.insert(...part.likelihoods...);
merged.candidate_prior_het.insert(...part.candidate_prior_het...);
merged.candidate_prior_hom_alt.insert(...part.candidate_prior_hom_alt...);
```

而 likelihood 行是**按候选下标索引**的。所以**不能只对 `merged.candidates` 去重**——
必须让这些并行数组同步跳过同一批条目，否则下标全部错位，会引入比原缺陷更严重的问题。
这是本轮没有直接动手改这一行的主要原因。

**下一步（两选一，均需先确认再动）**：
1. 在 merge 处做**成组去重**（候选与其并行数组一起跳过重复项），并确认保留的那一份是
   已配对的那个（本案例中 `merged` 里较早的那份 index=1 已配对，`part` 里较晚的那份
   index=5 是哨兵，故「跳过 incoming 重复」方向正确）；
2. 不改 merge，而在渲染/注释处**优先采用已配对孪生**（需跨 `result.candidates` 查找，
   因为已配对的孪生不在输出组内）。
两者都必须过本轮定的验收判据：两个窗口下该行 `RAW_MQandDP=28800,8` 且 `SB=0,0,3,3`，
且窗口 A 的 48 行保持与 GATK 逐字节相同。

## 最小复现

参考/输入夹具：`fixtures/chr20/ref20mnp.fasta`、`fixtures/chr20/mnp.bam`（仓库内已有）。

```bash
B=fastgatk-native/build/fastgatk-hc-call
COMMON=(-R fixtures/chr20/ref20mnp.fasta -I fixtures/chr20/mnp.bam \
        -L 20:10019901-10020710 --emit-ref-confidence GVCF --max-mnp-distance 1 \
        --threads 1 --add-output-vcf-command-line false)

$B "${COMMON[@]}"                          -O whole.g.vcf    # 48 条记录
$B "${COMMON[@]}" --stream-by-region 300   -O t300.g.vcf     # 48 条，但字段不同
$B "${COMMON[@]}" --stream-by-region 500   -O t500.g.vcf     # 20 条
```

## 证据 1：tile 尺寸扫描（OpenMP 与 Serial 结果一致）

| `--stream-by-region` | 记录数 | 相对非流式缺失 | 相对非流式多出 | 与非流式逐字节相同 |
| --- | ---: | ---: | ---: | --- |
| 100 | 20 | 28 | 0 | 否 |
| 200 | 41 | 8 | 1 | 否 |
| 300 | 48 | 0 | 0 | **否** |
| 405 | 41 | 8 | 1 | 否 |
| 500 | 20 | 28 | 0 | 否 |
| 810（=区间长度，单 tile） | 48 | 0 | 0 | **是** |
| 1000 | 48 | 0 | 0 | **是** |
| 1000000 | 48 | 0 | 0 | **是** |

要点：**只有 tile 覆盖整个区间时，流式输出才与非流式输出逐字节相同**。tile 变小时偏离出现，
且偏离程度随 tile 尺寸非单调变化。`--threads 1` 固定，排除线程数干扰。

## 证据 2：20:10020680 处注解读数

非流式（OpenMP 与 Serial 相同）：

```text
20  10020680  .  CA  AT,*,<NON_REF>  154.25  .  DP=8;ExcessHet=0.0000;MLEAC=1,1,0;MLEAF=0.500,0.500,0.00;RAW_MQandDP=28800,8  GT:AD:DP:GQ:PGT:PID:PL:PS:SB  1|2:0,4,2,0:6:84:1|0:10020679_AC_TA:249,84,99,165,0,159,252,99,168,264:10020679:0,0,3,3
```

流式 `--stream-by-region 300`：

```text
20  10020680  .  CA  AT,*,<NON_REF>  154.25  .  DP=8;ExcessHet=0.0000;MLEAC=1,1,0;MLEAF=0.500,0.500,0.00;RAW_MQandDP=97200,27  GT:AD:DP:GQ:PGT:PID:PL:PS:SB  1|2:0,4,2,0:6:84:1|0:10020679_AC_TA:249,84,99,165,0,159,252,99,168,264:10020679:0,0,0,0
```

差异字段：`RAW_MQandDP`（`28800,8` → `97200,27`）与 `SB`（`0,0,3,3` → `0,0,0,0`）。
`QUAL`、`PL`、`GT`、`GQ`、`PGT/PID/PS` 完全相同，即**基因型数值未受影响，受影响的是注释**。

## 证据 3：GATK 4.6.2.0 基准

```bash
third_party/jdk17/bin/java -Xmx1g -jar \
  third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar HaplotypeCaller \
  -R fixtures/chr20/ref20mnp.fasta -I fixtures/chr20/mnp.bam -L 20:10019901-10020710 \
  -ERC GVCF --max-mnp-distance 1 -O gatk.g.vcf --native-pair-hmm-threads 1 \
  --create-output-variant-index false --add-output-vcf-command-line false
```

- GATK 记录数 **48**，与非流式 native 相同。
- GATK 在 20:10020680 的行为
  `... RAW_MQandDP=28800,8 ... SB=0,0,3,3`，**与非流式 native 逐字段相同**。

结论方向明确：**非流式路径是对的（与 GATK 一致），`--stream-by-region` 是偏离的一方。**

## 与现有绿色测试的关系

`fastgatk-native/scripts/verify_hc_region_streaming.py`（CTest
`fastgatk-hc-region-streaming-contract`）断言：

- `--stream-by-region 1000000` 与普通输出**逐字节相同**；
- `--stream-by-region 500` 与普通输出**locus key 集合相同**。

该测试目前在双后端全绿，但它使用的是另一份夹具/区间；在
`fixtures/chr20/mnp.bam` + `20:10019901-10020710` 上，tile=500 的 key 集合
比非流式少 28 条。也就是说：**该契约在当前实现上并非普遍成立，而是依赖夹具**。
本次发现的用例应补进该 oracle，作为回归守卫。

## 根因分析（第 2 轮进展）

### 更正：先前「证伪拼接器」的推理有误

先前结论建立在**错误前提**上：我假设各 tile 产出的块与非流式的块相同，于是用非流式的 29 个块
去预测 `merge()`（`hc_call.cpp:5864`，PL 取逐元素最小值）的结果，得到
`PL=0,26,494 / GQ=26 / MIN_DP=12`，与实测的 `PL=0,0,0 / GQ=0 / MIN_DP=21` 不符，
据此判断「拼接器不是成因」。

该前提不成立：**各 tile 自己产出的块就已经是退化的（`PL=0,0,0`）**，两个退化块经 `merge()`
自然得到 `PL=0,0,0 / GQ=0`，与实测一致。而且从坐标上可判定拼接器**必然参与**：记录按 core
坐标输出，tile=500 的 core 边界在 10020401，而该记录为 `10020230–10020428`，跨越边界，
单个 tile 无法独立产出它。

正确表述：**拼接器是把两个已退化的块接起来的机制；真正的缺陷是 tile 在 core 边界附近算出的
参考置信度本身是错的**（DP 12–13 而非 31–32、`PL=0,0,0`）。

### 已验证：偏离集中在 tile 边界附近，且随 halo 变化（padding 不是修复）

同一夹具上做 tile × `--assembly-region-padding` 扫描，统计与非流式输出的差异行数：

| `--stream-by-region` | pad=100 | pad=200 | pad=300 |
| ---: | ---: | ---: | ---: |
| 300 | 1 行 | 1 行 | 1 行 |
| 405 | 43 行 | 43 行 | 3 行 |
| 500 | 43 行（仅 20 条记录） | 3 行 | 3 行 |
| 810（= 区间长度，单 tile） | **0 行（逐字节相同）** | 1 行 | 1 行 |

两个结论：

1. **默认 pad=100 下，只有单 tile（tile ≥ 区间长度）才与非流式逐字节相同**；任何真实分块都会偏离。
2. **加大 padding 不是修复**：它能把 tile=500 的块结构修回来（43 → 3 行），却把原本逐字节相同的
   单 tile 情形弄坏（0 → 1 行）。偏离只是被搬到了别处。

### 已验证：多 tile 下有两类不可消除的残留偏离（tile=500, pad=300）

**(D1) 流式路径凭空产生 phasing —— 本次新发现。** 非流式（= GATK）在 10020228/10020229 输出未定相：

```text
10020228  ...  0/1:4,8,0:12:99:302,0,129,314,153,467:3,1,5,3
```

流式输出为已定相，并新增 PGT/PID/PS：

```text
10020228  ...  0|1:4,8,0:12:99:0|1:10020228_G_A:302,0,129,314,153,467:10020228:3,1,5,3
```

`GT` 由 `0/1` 变为 `0|1`/`1|0` 并附带 PGT/PID/PS —— 这是**输出语义**的改变，不只是数值偏差。

**(D2) 注释 read 作用域错误。** 同前：10020680 处
`RAW_MQandDP` 由 `28800,8` 变为 `97200,27`，`SB` 由 `0,0,3,3` 变为 `0,0,0,0`。

进一步测量表明 **D2 是「分块相关」的，而非某个固定 tile 边界的局部效应**。
同一 locus 在不同分块下取值如下（默认 pad=100）：

| `--stream-by-region` | 10020680 的 `RAW_MQandDP` | `SB` | 与非流式一致 |
| ---: | --- | --- | --- |
| 非流式 | `28800,8` | `0,0,3,3` | 基准 |
| 200 | `97200,27` | `0,0,0,0` | 否 |
| 300 | `97200,27` | `0,0,0,0` | 否 |
| 405 | `28800,8` | `0,0,3,3` | 是 |
| 500 | `28800,8` | `0,0,3,3` | 是 |
| 600 | `97200,27` | `0,0,0,0` | 否 |
| 700 | `28800,8` | `0,0,3,3` | 是 |
| 810 | `28800,8` | `0,0,3,3` | 是 |

并且 `--assembly-region-padding` 也会翻转这个结果（tile=500 在 pad=100 正确、pad=300 错误）。
即：**流式路径的注释不满足「分块不变性」，而非流式路径按构造满足。**

已定位到的可疑代码位置（尚未证实，供修复用）：`run_region_streaming` 的每 tile 渲染调用

```cpp
auto result = fastgatk::calling::run(decoded.reads, references, tile_options);   // hc_call.cpp:6666
...
? gvcf(metadata_reader, result, ..., &decoded.reads,                              // hc_call.cpp:6674-6681
       tile_options.informative_read_overlap_margin, options.max_alternate_alleles)
```

`decoded.reads` 是 **tile 的 halo 读取批**，`tile_options` 携带该 tile 的
`interval_start/end` 与 `informative_read_overlap_margin`。注释证据集看起来由此二者共同决定，
因此随分块变化。注意不能简单归因为「batch 过大」：tile=810 的 batch 更大，却给出正确的 `28800,8`。

### 仍未完成：定位到具体代码

候选仍是参考块的两个构造器（精确的 `build_reference_blocks` vs 近似的
`build_profile_local_reference_blocks`），以及 `run_region_streaming` 的
`input_halo_intervals`（`hc_call.cpp:6463`）如何决定 tile 的读取窗口。

`FASTGATK_DEBUG_GVCF_EMISSION=1` 的输出显示 `symbolic-pl=unavailable` 在**两种模式下都出现**，
因此它**不是**判别依据（先前曾把这条列为优先线索，现更正）。

### 与状态文档旧记录的冲突

`IMPLEMENTATION_STATUS.md` 曾记录「流式 gVCF 输出路径已将实际 read batch 传给注释重算，
而非使用空 read 指针」，并据此视为已修正。本次结果表明：该改动**可能传入了非空但范围错误的
read batch（分区级而非 locus 级）**，因此"不再用空指针"并不等于"数值正确"。

## 已证明 / 未证明

**已证明**

- 非流式 native 与该区间上的 GATK 4.6.2.0 **全部 48 条数据行逐字段相同**（`diff` 无输出），
  包括 10020231–10020428 的 29 块碎片结构。这是本报告的基准：非流式路径是正确的一方。
- 在 `fixtures/chr20/mnp.bam` + `20:10019901-10020710` 上，`--stream-by-region N`
  （N < 区间长度）与非流式输出不一致，且**默认参数下不存在任何逐字节相同的多 tile 配置**。
- 偏离随 tile 与 `--assembly-region-padding` 变化：加大 padding 能修回块结构，却会弄坏
  单 tile 情形（见扫描表）——**padding 不是修复**。
- **两个可分别复现的残留缺陷**（tile=500, pad=300）：
  D1 在 10020228/10020229 凭空产生 phasing（`GT 0/1 → 0|1` 并新增 PGT/PID/PS）；
  D2 在 10020680 注释 read 作用域错误（`RAW_MQandDP`、`SB`）。
- OpenMP 与 Serial 的流式偏离完全相同（Host 路径问题，非 kernel）。
- 覆盖范围本身未丢失：tile=500 时 28 条被合并进 `10020230–10020428` 单块，
  10 条变异记录逐一保留。
- **更正**：先前「拼接器 `can_merge` 已被证伪」的结论作废（推理前提有误，见「更正」一节）。
  拼接器必然参与跨 core 边界的块合并；真正的缺陷在 tile 侧的参考置信度取值。

**未证明**

- 未定位到具体代码行；`build_profile_local_reference_blocks` 与
  `build_reference_blocks` 的取舍、以及 `input_halo_intervals` 的窗口口径都只是候选，
  未证实缺陷落在哪一处。
- **D1（凭空 phasing）与 D2（注释作用域）是否同源，尚未判定。**
- `RegionGvcfStitcher::can_merge` 只比较 `pl.size()` 而不比较 PL 值：本例中两个输入块
  都是 `PL=0,0,0`，故未暴露；在其他夹具上仍可能造成不正确的块合并，需单独验证。
- 未在该区间之外确定偏离边界（哪些夹具会触发、是否所有多重叠 indel 窗口都触发）。
- 未验证 `--stream-by-contig`、以及与 `--max-alternate-alleles`、read-overlap margin
  等选项的组合。
- 未验证下游影响（如 GenotypeGVCFs 对块粒度、`PL=0,0,0` 块、以及意外 phasing 的敏感度）。
- **修复尚未开始**：本报告不含任何代码改动。

## 复现脚本

`fastgatk-native/scripts/verify_gvcf_stream_overlapping_indels_gatk_oracle.py`
（第 1 轮 Track A 产出；当前状态为 **diagnostic**，即报告偏离而不判 pass/fail）。
尚未注册进 CTest；注册需在第 2 轮修复后把断言改为「流式 == 非流式 == GATK」。
