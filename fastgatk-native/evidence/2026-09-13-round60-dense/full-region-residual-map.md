# 第 60 轮：dense 模式在**整段真实语料**上的残差地图

此前只验过 3 kb 窗口（第 52–57 轮），窗口恰好全绿。本轮把比对扩到整个语料区间，
量出真实残差规模与分类。

## 方法与输入

| 项 | 值 |
| --- | --- |
| 输入 | GATK 自带 chr20 gVCF（1291 条，`.../haplotypecaller/expected.testGVCFMode.gatk4.g.vcf`） |
| 参考 | `testdata/chr20/reference/GRCh37.chr20.fa` |
| 区间 | `-L 20:10000000-10099999`（覆盖全部记录） |
| 模式 | `--include-non-variant-sites`（dense），两侧同参数，native 加 `--gatk-compatible-annotations` |
| 耗时 | GATK 5.3 s；native 25.9 s |

## 总量

| | GATK | native |
| --- | --- | --- |
| 行数 / 位点数 | 100000 / 100000 | 99993 / 99993 |
| 仅 GATK 有的位点 | — | **7** |
| 仅 native 有的位点 | — | 0 |
| 共有但行不同的位点 | — | **153** |

即整段 100 kb 上残差为 **160 个位点（0.16%）**，且没有任何「只存在于 native」的位点。

## 残差分类（按 GATK 形状 / native 形状 / 不同列）

| # | GATK | native | 不同的列 | 条数 | 典型位点 |
| --- | --- | --- | --- | --- | --- |
| A | `*` | `*` | 7（INFO） | **109** | 10004770 |
| B | `*` | `.`（REF-only） | 4,5,6,7,8 | **30** | 10097443 |
| C | `.` | `.` | 7（INFO） | 9 | 10008964 |
| D | `.` | `.` | 5,7（QUAL+INFO） | 3 | 10024301 |
| E | `.` | `.` | 7,8（INFO+FORMAT） | 2 | 10077008 |
| F | GATK 有、native 无 | — | — | **7**（上表「仅 GATK」） | 10062936 |

### 各类的实际差别

**A（109 条，占 71%）——native 在 `*` 行上多写了 rank-sum/MQ 注释。**
```
GATK   N  *  0  LowQual  AC=1;AF=0.500;AN=2;DP=80;ExcessHet=0.0000;FS=3.758;MLEAC=1;MLEAF=0.500;QD=-0.00;SOR=1.306
NATIVE N  *  0  LowQual  AC=1;AF=0.500;AN=2;BaseQRankSum=0.599;DP=80;...;MQ=52.80;MQRankSum=-6.625e+00;QD=-0.00;ReadPosRankSum=1.88;SOR=1.306
```
GATK 在合成/物化的 `*` 行上只发布它重算的注释（AC/AF/AN/DP/ExcessHet/FS/MLEAC/MLEAF/QD/SOR），
native 额外带上了 `BaseQRankSum`/`MQ`/`MQRankSum`/`ReadPosRankSum`。
**这一类单个修复即可消掉约 2/3 的残差**，是下一轮的首选。

**B（30 条）——GATK 保留 `*` 行，native 变成 REF-only 行。**
```
GATK   N  *      0    LowQual  AC=1;AF=0.500;AN=2;DP=153;...;FS=1.851;...;QD=-0.00;SOR=0.948  GT:AD:DP:GQ:PGT:PID:PL:PS
NATIVE N  .      .    .        DP=63                                                          GT:DP:RGQ
```
即该位点 native 走的是参考块/merge-null 物化路径，而 GATK 认为跨接删除仍被拥有（`*` 保留）。

**C（9 条）——REF-only 行上 native 多写注释**（同类 A 的 REF-only 版本）。

**D（3 条）——monomorphic-recovery 的 QUAL 不同**：GATK `163.67` vs native `Infinity`（同一形状的另一侧也出现过反向）。

**E（2 条）——REF-only 行上 native 多写了 `PGT:PID:PS`**。

**F（7 条）——GATK 输出 `QUAL=Infinity` 的 REF-only 行，native 完全没有该行**：
```
20  10062936  .  N  .  Infinity  .  DP=89;MLEAC=.;MLEAF=.  GT:DP:RGQ  0/0:23:16
```
与 D 同源（同一 QUAL/发射判定路径），7 个位点集中在两处（10062936–10062938、10087821 附近等）。

## 结论

- 整段语料的 dense 对齐度 ≈ **99.84%**（160/100000 位点残差），且**零**额外位点；
- 残差集中在 5 个可判定的类里，其中 A 类占 71% 且看起来是**单点修复**；
- 3 kb 窗口之所以全绿属采样偏差，**不应**用它代表整段（本轮更正）。

## 第 61 轮：A 类的**补充定位**（把修复点钉死）

对 109 行 `*`/`*` 逐行比较 INFO 键集合：

| 差 | 键 | 条数 |
| --- | --- | --- |
| native 有、GATK 无 | `MQ` | **109**（全部） |
| native 有、GATK 无 | `BaseQRankSum` / `MQRankSum` / `ReadPosRankSum` | 各 86 |
| GATK 有、native 无 | `QD` | 24 |

并确认这些行**确实由本轮新增的跨位点物化 pass 产生**：以 10004770 为例，附近唯一的输入记录是
`10004769 TAAAACTATGC > T,<NON_REF>`（跨 10004769-10004778），该位点没有自己的记录，
native 的 `*` 行来自 `materialize_spanning_loci()`；GATK 同样发布该行但不带源记录的
`MQ`/`BaseQRankSum`/`MQRankSum`/`ReadPosRankSum`（它只发布在该位点重算的注释）。

**一次被证否的修复（记录在案）**：在 `materialize_spanning_loci()` 里 `bcf_dup` 之后直接清除这四个
INFO 键——**实测无效**（整段残差仍为 153，109 行照旧）。原因是这些注释在后续的
merge/注释阶段被（从源记录的 Host 侧状态）重新写回。因此修复必须落在**输出边界**：
在编码/注释阶段对 `record.materialized_spanning_locus == true` 的行抑制这四个键，
而不是在 pass 里清除。该改动已回退，不留未经验证的代码。

## 第 63 轮：余下 40 行 `*`/`*` 的机制（多删除 ALT 源的塌缩）

以 10008953 为例：

```
输入 10008952:  CACACACACACACA > C,CCA,CCACACACACA,CCACACACACACA,<NON_REF>
                INFO: MLEAC=0,1,0,1,0   → 被调用的 ALT 是 #2(CCA) 与 #4(CCACACACACACA)，
                两者长度都 < REF(14) ⇒ **两个被调用的删除等位基因**
GATK   (10008953): N  *  0  LowQual  AC=2;AF=1.00;AN=2;…;QD=-0.00  GT:AD:DP:GQ:PL  1/1:1,11:34:62:1184,62,0
NATIVE (10008953): N  *  0  LowQual  AC=0;AF=0.00;AN=0;…（缺 QD）  GT:AD:DP:GQ:PL  ./.:1,14:34:0:1717,539,803
```

**机制**：`materialize_spanning_loci()` 的 `called_deletion_allele_index()` 只返回**第一个**被调用的删除
等位基因，投影目标因此是 `{REF, 那一个删除}`；而源 GT 是 `2/4`（两个删除各一份），
另一份落在目标列表之外，于是 `remap_record_to_allele_union()` 里既有的
`force_no_call_on_dropped_gt` 逻辑把该样本整体判成 **no-call** —— 这正是 native 输出
`./.`/`AC=0;AN=0`（并因 depth=0 缺 `QD`）的原因。

**GATK 的做法**：merger 把该跨接记录的**所有**删除等位基因都塌缩成**同一个**符号 `*`
（所以 GT 变成 `1/1` 两份都是 `*`），PL 则是对被塌缩的那些旧基因型组合做 log 域归并
（因此 GATK 的 `1184,62,0` 与 native 取单个子矩阵的 `1717,539,803` 不同）。

**修复配方（下一步）**：
1. 投影目标改为多对一：所有「删除类」ALT（长度 < REF）都映射到同一个 `*`，非删除 ALT 映射为 NO_CALL；
2. GT 逐拷贝按上述映射重写（不再是 1:1 的名称匹配）；
3. PL/AD 用**塌缩后的等价基因型**做 log 域归并（GATK 的 `newToOldGenotypeMap` 语义），
   而不是取某一个代表等位基因的子矩阵；
4. 完成后用整段语料复验（预期余下 `*`/`*` 行 40 → 0，整段残差 84 → 44）。

## 第 64 轮：塌缩模型被实测**逐个否证**——真正的原因是「要与覆盖该坐标的参考块合并」

按第 63 轮的配方先做模型验证（用 10008952 的真实数据反推 GATK 的 `1184,62,0`）：

```
源记录 10008952: alleles = [CACACACACACACA, C, CCA, CCACACACACA, CCACACACACACA, <NON_REF>]
                 GT 2/4 ; AD = 1,3,14,5,11,0
                 PL(21) = 1717,890,1216,539,769,803,739,496,154,687,595,342,0,443,533,1149,1060,712,761,607,1292
具体删除等位基因（长度 < REF 且非符号）: 索引 1,2,3,4
```

| 候选模型 | 预测 PL | 与 GATK `1184,62,0` |
| --- | --- | --- |
| 取单个代表等位基因的子矩阵 `{REF,CCA}` | `1717,539,803` | ✗（这正是 native 现状） |
| 取 `{REF,CCACACACACACA}` 子矩阵 | `1717,595,533` | ✗ |
| 把所有删除类等位基因做 log 域塌缩归并（log-sum-exp） | `1717,539,0` | ✗ |
| 同上但用线性域求和 | `1717,539,0` | ✗ |
| 上述两种再排除 `<NON_REF>` 基因型 | `1717,539,0` | ✗ |

**结论（更正第 63 轮的配方）**：GATK 在该位点的 PL **不是**源记录 PL 的任何塌缩，
而是**合并**的结果——该坐标同时被**同一张本的参考块记录**覆盖，GATK 的
`ReferenceConfidenceVariantContextMerger.mergeRefConfidenceGenotypes()` 会把
「跨接变异记录」与「覆盖该坐标的参考块记录」一起合成，PL/AD/GT 因而来自**两者**；
native 的 `materialize_spanning_loci()` 只克隆了变异记录，所以数值不同，
且我那次单等位基因投影又触发了 dropped-GT 逻辑把样本判成 no-call。

**因此该类的正确修法**（下一步，工作量明显大于配方 4 步）：合成位点时必须模拟
merger 的合并语义——把**所有覆盖该坐标的记录**（跨接变异记录 + 参考块记录）按
`mergeRefConfidenceGenotypes()` 的规则合成一份样本数据（含 `<NON_REF>` 维度的重映射），
而不是从单条源记录做投影。第 60 轮残差地图里的 B 类（30 条 `*` vs REF-only）大概率同源。

> 方法教训：先做**模型验证**（用真实数据反推目标值）再写实现，可以避免像第 61 轮那样
> 先写再被证否；本轮用 5 个候选模型换来一个更深的正确认识。

## 第 65 轮：第 64 轮的「与参考块合并」假设被**证伪**，投影模型仍未对上

**证伪**：直接枚举覆盖 10008953 的输入记录——该坐标**没有**参考块覆盖，唯一覆盖它的是
`10008952 CACACACACACACA > C,CCA,CCACACACACA,CCACACACACACA,<NON_REF>`（跨 10008952-10008965）。
周边记录为 `10008945/10008950/10008951` 的块（END 分别止于 10008947/10008950/10008951）与
`10008948 TA > T,<NON_REF>`（跨 10008948-10008949），**都不覆盖 10008953**。
故第 64 轮把该位点解释为「与参考块合并」是错的。

**模型仍不对上**：在该单记录投影下再试 `<NON_REF>` 的三种归并方式

| 模型（删除类 {C,CCA,CCACACACACA,CCACACACACACA} → `*`） | 预测 PL | 与 GATK `1184,62,0` |
| --- | --- | --- |
| `<NON_REF>` → 参考 | `1149,539,0` | ✗（首值接近、中值差很多） |
| `<NON_REF>` → `*` | `1717,539,0` | ✗ |
| `<NON_REF>` → 丢弃（NO_CALL） | `1717,539,0` | ✗ |

**结论**：GATK 在 `*` 位点的样本投影不是「按索引映射后做 log-sum」这么简单；
`ReferenceConfidenceVariantContextMerger.mergeRefConfidenceGenotypes()` 还会经
`getIndexesOfRelevantAllelesForGVCF` / `newToOldGenotypeMap` / `generateAD` 与
`GenotypeLikelihoods` 的归一化，其中对 NO_CALL 拷贝是否按参考处理、以及 AD 如何
在塌缩后的等位基因上合并，都还没有实测钉住。

**下一步（本轮给出的方法与判据）**：不要再靠反推猜模型，改为**直接对 pinned jar 里的 merger 做探针**
（沿用本会话第 49 轮用过的 Java 探针手法）：构造该条记录，调用
`ReferenceConfidenceVariantContextMerger.merge(...)` 与 `mergeRefConfidenceGenotypes(...)`，
打印返回 VC 的 alleles/GT/AD/PL，与 GATK 输出 `1184,62,0` 对齐后再写 native 实现。

## 第 66 轮：C/E 类（REF-only 行）的两次尝试与结论

目标是复用第 62 轮的抑制器消掉 C 类（9 条 REF-only 行的多余 INFO）与 E 类（2 条 REF-only 行
多余的 `PGT:PID:PS`）。两次尝试与实测：

| 尝试 | 结果 |
| --- | --- |
| 把抑制条件扩到 `finalized_monomorphic_ref`（只挂在两个「最终」出口） | 整段残差仍 84（C 类 9 条照旧）——**该路径不经过这两个出口** |
| 再把抑制器挂到两个 `finalized_monomorphic_ref` 早退出口，并额外清除 `PGT/PID/PS` | E 类的**FORMAT 部分确实修好**（那 2 行从「INFO+FORMAT 不同」变为「仅 INFO 不同」），但 C 类 INFO 仍不为所动（9 → 11，总数仍 84） |

**结论**：C/E 类所在的 REF-only 行**不是** `finalized_monomorphic_ref`，而是 dense 模式下
**参考块逐坐标展开**产生的行（`reference_block = true`，见解码阶段 `include_non_variant_sites && reference_only`
分支）。因此：
1. 对它们的抑制必须以「块展开行」为判据（例如新增 `Record::materialized_reference_only` 标志），
   而不是 `finalized_monomorphic_ref`；
2. 清除 `PGT/PID/PS` 的那一步是**正确**的（实测那 2 行的 FORMAT 因此与 GATK 一致），
   但必须与 INFO 抑制一起落地并上锁，单独一半会把「INFO+FORMAT 不同」变成「仅 INFO 不同」而总数不变。

两次改动均已**回退**（不留未验证、未上锁的代码）；结论与判据留在本节，供下一步一次性做完。

## 第 67 轮：C/E 类的**第三次**尝试同样无效——须先查清这些行的真正来源

按第 66 轮的判据新增 `Record::materialized_reference_only`（在解码阶段的逐坐标展开处打标），
并把抑制器与 `PGT/PID/PS` 清除都挂到该标志上——**整段残差仍为 84，C/E 类一行未变**。

即：C/E 类那 11 条 REF-only 行**既不**是 `materialized_spanning_locus`，**也不是**
`finalized_monomorphic_ref`，**也不是**解码阶段的块逐坐标展开行。聚合路径里还剩两个候选来源：

- `split_reference_blocks_at_variants()`（`:3815` 处 `split.reference_block = true`）产生的分段行；
- 分组阶段调用 `materialize_reference_only()` 的那些行（`:6604` / `:7311`）。

**下一步（先查来源，再谈修复）**：不要继续猜标志。加一个**临时插桩**——在输出前对
`record.pos` 命中这 11 个位点时打印该记录的
`materialized_spanning_locus / materialized_reference_only / finalized_monomorphic_ref /
reference_block / alleles / allele_count`，一次运行即可确定来源；
确认后把抑制挂到正确的判据上，并用整段语料复验（预期 84 → 73）。三次尝试均已回退。

## 第 68 轮：用插桩把 C/E 类的候选标志**逐个排除**

在抑制器入口插桩（`FASTGATK_DEBUG_MAT=1`，打印 `pos / span / mono / block / alleles / alt`），
跑 `-L 20:10008960-10008970` 得到：

```
[MATDBG] pos=10008965 span=0 mono=0 block=0 alleles=1 alt=.
[MATDBG] pos=10008966 span=0 mono=0 block=0 alleles=1 alt=.
[MATDBG] pos=10008967 span=0 mono=0 block=0 alleles=1 alt=.
...
```

即 C/E 类这些 REF-only 行：

- **不是** `materialized_spanning_locus`（span=0）——第 66 轮尝试已否证；
- **不是** `finalized_monomorphic_ref`（mono=0）——第 66 轮尝试已否证；
- **不是** `reference_block`（block=0）——第 67 轮新增的 `materialized_reference_only` 也因此挂不上（它在解码展开处打标，而那些行 block=0）；
- `alleles=1`（ALT='.'），且 10008965 是 `10008952` 那条多等位记录跨度的**最后一个坐标**，
  10008966 起才有输入块记录。

**因此下一轮的正确做法**：插桩要打在**创建点**而不是汇合点——在三个可能产生 REF-only 行的位置
各打一个不同的标记（解码阶段的块展开、`split_reference_blocks_at_variants()` 的分段、
分组阶段的 `materialize_reference_only()` 调用），一次运行即可看出这些行由哪一处产生；
确认后再把抑制挂到正确判据上（预期整段残差 84 → 73）。插桩已回退。

## 第 69 轮：创建点插桩的**尝试失败**（注入的标记在其中一个调用点编译不过），锚点已记下

按第 68 轮的计划在三个创建点注入标记，编译报错：
`genotype_gvcf_tool.cpp:6606: 'record' was not declared in this scope`。
原因：两个 `materialize_reference_only(...)` 调用点的参数名不同——正确锚点为

| 位置 | 实际写法 |
| --- | --- |
| `:6604`（分组阶段，含具体记录变量） | `materialize_reference_only(output_header, record, options.gatk_annotation_compatibility);` |
| `:7311`（另一分组分支） | 形参是 `records[index]`，**不是** `record`（故用 `record.pos` 的标记在此处编译失败） |

其余两个创建点（解码阶段的逐坐标展开 `Record expanded = staged;`、`split_reference_blocks_at_variants()`
里的 `split.reference_block = true;`）注入成功且能编译。

**下一步（锚点已给全，可直接照做）**：在 `:6604` 与 `:7311` 分别用各自的变量名注入标记
（前者用 `record.pos`，后者用 `records[index].pos`），连同另两处一起跑一次即可定位来源；
随后把抑制挂到正确判据上并复验（预期整段残差 84 → 73）。第 68 轮的排除结论
（这些行 span=0/mono=0/block=0）不受影响。

## 第 70 轮：创建点插桩=第二次失败（标记落在变量不可见的作用域），结论与建议

按第 69 轮记录的锚点注入四处标记，编译仍失败，且错误全部落在**注入的标记行**上：

```
6606:75: error: 'record' was not declared in this scope
7315:74: error: no match for 'operator[]' (std::vector<Record> and <unresolved overloaded function type>)
```

说明这两处的「调用首行」锚点实际位于**变量不可见的嵌套作用域/续行**上（两条
`materialize_reference_only(...)` 都是跨行调用，匹配到的行并不是语句起始作用域）。

**建议**：下一次不要靠文本锚点盲插，改为**用 `nm`/`gdb` 之外的两条更稳的路子**任一：
1. 直接读代码把两处 `materialize_reference_only` 调用的**完整作用域**读清楚（两条调用前后各 ±10 行），
   再决定标记插在哪一行；或
2. 反过来做——**不改代码**，用「哪一处被调用」的侧面证据：
   把 `materialize_reference_only()` 的函数体第一行打标记（函数级插桩，作用域一定正确，
   只是无法区分调用者），配合 `-L` 窗口逐步缩小范围（例如把区间缩到只含 10008965-10008970）。

函数级插桩是**零风险**的一步（作用域一定正确），建议下一轮先用它确认「这些行确实由
`materialize_reference_only()` 产生」，再用 ±10 行读取把调用者钉死。

## 第 71 轮：C/E 来源=第四次未命中，给出当前最稳的定位手段

又试了「输出等位基因子集里 `uncovered` 的 REF-only 物化分支」这一假设（在该分支打标 + 抑制），
**整段残差仍为 84，C/E 一行未变**。至此 C/E 的候选来源逐一被排除：

| 假设 | 轮次 | 结果 |
| --- | --- | --- |
| `materialized_spanning_locus`（跨位点物化行） | 66 | ✗ 打标后无变化 |
| `finalized_monomorphic_ref`（含早退出口） | 66 | ✗ |
| 解码阶段块逐坐标展开 | 67 | ✗（且汇合点插桩实测 span=0/mono=0/block=0） |
| `split_reference_blocks_at_variants()` 分段 | 69/70 | 未能安全插桩（锚点作用域问题） |
| 子集内 `uncovered` REF-only 物化分支 | 71 | ✗ |

**当前最稳的定位手段（下一步，作用域一定正确）**：把标记打在**编码器接收处**
（`[&](GenotypeComputed computed)` 或 `[&](GenotypeEncoded encoded)` 分支的第一行）——
那里 `computed.record` 一定在作用域内，打印
`pos / materialized_spanning_locus / materialized_reference_only / finalized_monomorphic_ref /
reference_block / alleles / INFO 键列表`；一次 11 坐标的小窗口运行即可看到这些行在**输出时刻**的
全部状态（本会话已有 `FASTGATK_DEBUG_SUBSET` / `FASTGATK_DEBUG_QD` / `FASTGATK_DEBUG_MAT` 三个
同类插桩先例，都能编译运行）。查询器插桩在本代码库里被反复证明比「猜创建点」可靠。

**注意**：第 68 轮汇合点（抑制器入口）的插桩已给出 `span=0 mono=0 block=0 alleles=1`，
说明这些行在到达抑制器时**三个标志全为 0**；因此**要么**它们从未被打标，**要么**标志在中途被清
（例如 `reference_block = false` 这类赋值）。编码器处的插桩会顺带回答这个问题。

## 第 72 轮：编码器处插桩**成功**——确认三标志全 0 且带 INFO，候选来源收敛到两处

在编码器接收处（`[&](GenotypeComputed computed) -> std::optional<GenotypeEncoded>`，两处遍历各一个）
插桩 `FASTGATK_DEBUG_ENC=1`，跑 `-L 20:10008960-10008970`：

```
[ENC] pos=10008965 span=0 mono=0 block=0 nall=1 a1=. info=yes
[ENC] pos=10008966 span=0 mono=0 block=0 nall=1 a1=. info=yes
...（10008965-10008970 各一条）
```

即这些行在**输出时刻**同样 `span=0 mono=0 block=0`、`allele_count=1`、且**带 INFO**。
结合第 71 轮「在子集 `uncovered` 分支打标无效」，可判定这些行**不经过**该分支。

**剩余的两个候选（下一步只需分别验证一次）**：
1. 解码阶段的块逐坐标展开——但第 67 轮在该处打标**有效设置**后残差未变，
   故除非标志在那之后被清（`reference_block = false` 这类赋值），否则可排除；
2. `split_reference_blocks_at_variants()` 的分段行——它设 `reference_block = true`，
   而输出时刻 `block=0`，说明中途有清标志的赋值；**这是目前最可能的来源**。

**下一步的具体做法（作用域已确认安全）**：在 `:3815` 的 `split.reference_block = true;`
**之后**插一行读回打印（`fprintf(stderr, "[SPLIT] pos=%d\n", split.pos)`），
再在编码器处打印 `block`；若分割行在编码器处 `block=0`，则确认「标志被中途清除」，
随后把抑制挂到**分段行**的判据上（例如给分段行单独加一个不会被清的标志）。

## 第 73 轮：`split_reference_blocks_at_variants()` 也被排除——C/E 来源仍待定

同一次运行里同时插桩两处（分段处 `[SPLIT]` 打印、编码器处 `[ENC]` 打印），
窗口 `-L 20:10008960-10008970`：

```
[ENC] pos=10008965 span=0 mono=0 block=0 nall=1
...（每坐标一条）
（没有任何 [SPLIT] 行）
```

即 **该窗口内 `split_reference_blocks_at_variants()` 根本没有产生分段**，第 72 轮把来源收敛到
split 分段的判断**被证伪**。至此 C/E 的候选来源及其实测结果：

| 候选 | 验证方式 | 结果 |
| --- | --- | --- |
| 跨位点物化行（`materialized_spanning_locus`） | 打标 + 整段复验 | ✗ 残差不变 |
| `finalized_monomorphic_ref` | 打标 + 早退出口 | ✗ 残差不变 |
| 子集 `uncovered` REF-only 分支 | 打标 + 整段复验 | ✗ 残差不变 |
| 解码阶段块逐坐标展开 | 打标 + 整段复验 | ✗ 残差不变 |
| `split_reference_blocks_at_variants()` 分段 | 运行期打印 | ✗ **该窗口内根本没跑** |

**关键事实**：10008966 **有**输入块记录（`10008966 C <NON_REF>`），dense 模式下块会被逐坐标展开，
所以这一行**必然**由某处展开逻辑产生，但输出时刻 `block=0`、`allele_count=1`。
四个候选已被排除 ⇒ 说明还有**第五处** REF-only 生成路径尚未定位（可能就在解码阶段块展开的
另一个分支，或展开后又经 `materialize_reference_only()` 走了 `reference_block = false` 的路径，
而该函数被调用时传的是**分组后的副本**，导致打标对象与输出对象不是同一个）。

**下一步（最直接）**：在 `materialize_reference_only()` **函数体第一行**无条件打标
（函数级插桩，作用域必然正确——这是第 70 轮就提出、至今未试的零风险手段），
跑同一窗口即可确认它是否被调用；若被调用，再在其两处调用点附近±10 行读取完整作用域，
用正确变量名插桩定位调用者。

## 第 74 轮：函数级插桩**命中**——C/E 行确由 `materialize_reference_only()` 产生

按第 73 轮的建议做函数级插桩（在 `materialize_reference_only()` 函数体第一行无条件打印），
窗口 `-L 20:10008960-10008970`：

```
[MRO] pos=10008965 block=1 nall=2
[MRO] pos=10008966 block=1 nall=2
[MRO] pos=10008967 block=1 nall=2
[MRO] pos=10008968 block=1 nall=2
[MRO] pos=10008969 block=1 nall=2
（共 5 次调用）
```

**结论**：这些 C/E 行**确实**由 `materialize_reference_only()` 产生——进入时 `block=1, allele_count=2`
（纯 `<NON_REF>` 参考块），输出时刻变成 `block=0, allele_count=1`（REF-only）✓ 与编码器观测一致。

**因此第 71 轮的做法错在分支**：该函数内部有两条路径，`uncovered` 分支（第 71 轮打标处）与
**正常路径**（紧凑 PL=[0]、保留 DP/AD/GQ 等）；这些行走的是**正常路径**，所以打标没生效。

**下一步（这次有确切位置，一步可成）**：把标志设在 `materialize_reference_only()` 的
**函数入口**（或在正常路径的分支里），再把第 62 轮的抑制器条件扩到该标志；
预期：C 类 9 行与 E 类 2 行的多余 INFO/`PGT:PID:PS` 被清除，整段残差 **84 → 73**。
（另需注意：第 67 轮曾在解码展开处打标却无效，说明该处的展开记录与最终输出记录不是同一对象
路径，故标志要设在**函数内部**而不是调用点。）

## 第 75 轮：标志法**第五次**无效 ⇒ 键是在抑制器**之后**被重新写入的

把标志设在 `materialize_reference_only()` **函数入口**（第 74 轮插桩确认的位置），
并把抑制器条件扩到该标志（含 `PGT/PID/PS` 清除）：**整段残差仍 84，C/E 一行未变。**

至此「标志 + 计算阶段出口抑制」这条路已在**五个位置**失败：
跨位点物化行、`finalized_monomorphic_ref`、子集 `uncovered` 分支、解码展开调用点、
`materialize_reference_only()` 函数入口。而第 62 轮对**跨位点物化行**用同一手法是**成功**的
（109 行减少）——两者唯一的机制差别是：跨位点行由本工具的 pass → merge → 计算阶段产生并输出，
而 C/E 行由**参考块展开**产生。

**这一对比指向一个明确结论**：C/E 行的那四个 INFO 键是在我的抑制器**之后**被写回的
（最可能是编码/文本阶段的 GATK 兼容处理，或参考块路径特有的注释补写），
因此修复点不在计算阶段，而在**编码/文本写出阶段**。

**下一步（一次实验即可判定）**：在编码器接收处（`[&](GenotypeComputed computed)` 首行，
第 72 轮已验证该处插桩可编译运行）打印 `computed.record.value` 的 **INFO 键列表**：
- 若键在编码器处**已不存在**、却在输出文本里出现 ⇒ 确认为文本阶段写回，修复点定在
  `gatk_compatible_record_text()` / `apply_gatk_annotation_compatibility()`；
- 若键在编码器处**仍存在** ⇒ 抑制器根本没跑到这些行（即它们走了某个不经过两个出口的返回路径），
  需再定位该返回路径。

## 第 76 轮：**方法学更正**——窄窗口探针根本没有复现 C/E 类行

在编码器处打印 INFO 键列表，窗口 `-L 20:10008960-10008970`：

```
[ENC] pos=10008965 nall=1 ninfo=4 keys: END DP RCQ RCP
[ENC] pos=10008966 nall=1 ninfo=4 keys: END DP RCQ RCP
...
（输出文本：`10008966  DP=56` 等）
```

即这些行**不含任何 rank sum**（只有 `END/DP/RCQ/RCP`，输出只剩 `DP`），
而类 C/E 的定义正是「native 多写了 `MQ`/`BaseQRankSum`/`MQRankSum`/`ReadPosRankSum`」。

**结论（重要更正）**：`-L 20:10008960-10008970` 这样的**窄窗口**与整段运行（`-L 20:10000000-10099999`）
的**块/覆盖上下文不同**（重叠参考块被区间裁剪），因此窄窗口里那些 REF-only 行**不是**残差地图里的
C/E 行。也就是说：第 68–75 轮所有"窄窗口插桩 → 判定 C/E 行标志"的结论，**观察到的都不是 C/E 行**
（这也解释了为什么五个位置的打标都"无效"——那些行本来就不是要修的行）。

**下一步（必须改正探针条件，一次运行即可）**：把插桩条件从 `pos ∈ 窗口` 改为
`pos == <某个真正的 C/E 位点>`（例如 `10008964` / `10077008`），并**用整段命令运行**
（`-L 20:10000000-10099999`，约 26 s）；在该条件下打印编码器处的 INFO 键列表与各标志，
即可得到 C/E 行的真实状态。第 74 轮的函数级插桩同样要在整段命令下重跑复核。

## 第 77 轮：改正探针后**一次命中**诊断，修复仍差最后一步

按第 76 轮的更正，把插桩条件改成 `pos == 10008964 || pos == 10077008` 并用**整段命令**运行
（`-L 20:10000000-10099999`）：

```
[ENC] pos=10008964 span=1 mono=1 block=0 nall=1 ninfo=8
      keys: BaseQRankSum DP ExcessHet MLEAC MLEAF MQRankSum RAW_MQandDP ReadPosRankSum
[ENC] pos=10077008 span=0 mono=0 block=0 nall=1 ninfo=4 keys: END DP RCQ RCP
```

**结论（C 类，10008964）**：这类行同时是 `materialized_spanning_locus=1` **和**
`finalized_monomorphic_ref=1`，所以在计算阶段走的是**`finalized_monomorphic_ref` 早退分支**——
而第 62 轮加入的抑制器**只挂在两个最终出口**上，早退分支根本不经过它。**这就是前七轮"打标无效"的真正原因**
（不是标志没设上，而是**调用点没覆盖这条路径**）。

**修复尝试（本轮）**：在早退分支加一行 `suppress_materialized_spanning_annotations(...)`，
并把抑制条件扩到 `finalized_monomorphic_ref`（同时清 `PGT/PID/PS`）→ **整段残差仍 84，C 类 11 行未变**。

⇒ 说明**清除动作本身在这些行上没生效**（而非调用点问题）。最可能的原因：这些记录的 INFO
在那一刻尚未 `bcf_unpack(BCF_UN_INFO)`，或 `bcf_update_info_float(..., nullptr, 0)` 对它们返回 -1
（我在实现里一直用 `(void)` 忽略了返回值）。**下一步（很具体）**：
1. 在抑制器里打印每次 `bcf_update_info_float(...)` 的返回值与 `record.value->n_info`（前后对比），
   确认删没删掉；
2. 若返回 -1，改用 `bcf_unpack(record.value, BCF_UN_INFO)` + 按声明类型调用
   `bcf_update_info_float/int32`，或直接在**编码阶段**（`GenotypeEncoded` 分支，文本格式化之前）
   处理——那里的记录同样是 `computed.record.value`，且第 62 轮对跨位点行的成功案例证明
   「在正确路径上清除」是有效的。

## 第 78 轮：抑制失效的**根因**找到——`bcf_update_info_*` 的删除调用返回 0 但**不删除**

按第 77 轮的两步走：在抑制器里先 `bcf_unpack(record.value, BCF_UN_INFO)`，再打印每次清除的返回值与
`n_info`；同时在 `finalized_monomorphic_ref` 早退分支补上抑制调用（该路径第 77 轮已确认必须覆盖）。
整段运行、只看 10008964：

```
[MQrm] tag=MQ             rc=0 ninfo=8
[MQrm] tag=BaseQRankSum   rc=0 ninfo=8
[MQrm] tag=MQRankSum      rc=0 ninfo=8
[MQrm] tag=ReadPosRankSum rc=0 ninfo=8
```

**四次调用全部 `rc=0`（报告成功），而 `n_info` 始终是 8——一条都没删掉。**
这解释了前面所有"打标/条件都对、残差却不动"的现象：**不是判据问题，而是清除动作本身无效**。

（整段残差仍 84，与之一致。）

**下一步（很具体）**：
1. 不要再用 `bcf_update_info_float(hdr, line, tag, NULL, 0)` 这条路去"删除"——
   在此 htslib 版本上它显然不是删除语义（返回 0 但 `n_info` 不变，更像是把值置为 missing）；
2. 改用**文本层过滤**：native 的 GATK 兼容输出有自己的文本组装/格式化路径
   （`format_gatk_float_value()` 与 `gatk_compatible_record_text()`），其中已按 **key** 分派
   （例如 `key == "MQ" || key == "QD"` 决定精度）；在那里对「物化行/REF-only 行」按 key 丢弃
   这四个键**最直接**，且不受 htslib INFO 删除语义影响；
3. 判据可以随行携带（本会话已有的两个标志不足以覆盖，需要在文本层知道该行是物化行——
   可由记录的两个标志组合或再加一个标志传入文本组装函数）。

## 第 79 轮：文本层过滤的**具体改法**（读码得到，下一步照做即可）

已定位到需要改的三处（行号为当前树）：

| 位置 | 现状 | 需要的改动 |
| --- | --- | --- |
| `:4934` | `std::string gatk_compatible_record_text(const std::string& formatted)` | 增加一个默认参数，例如 `bool drop_read_level = false` |
| `:4987` 附近 | 组装 INFO 的循环：`info_text << '=' << format_gatk_float_value(ordered[index].key, ordered[index].value);` | 在该循环内按 key 跳过：`drop_read_level && key ∈ {MQ, BaseQRankSum, MQRankSum, ReadPosRankSum}` → `continue` |
| `:6728`、`:7493` | 两处调用 `gatk_compatible_record_text(<formatted>)` | 传 `computed.record.materialized_spanning_locus \|\| computed.record.finalized_monomorphic_ref`（REF-only 行另需丢 `PGT/PID/PS`） |

原因：该函数**只收到已格式化文本**，拿不到"这一行是物化行"的信息，所以必须由调用点把该状态传进来。
改完后用整段命令复验（预期 84 → 73：C 类 9 行 + E 类 2 行的多余键消失），再跑双后端全量。

**为什么不用 htslib 删除**（第 78 轮实测）：`bcf_update_info_float(hdr, line, tag, NULL, 0)`
四条调用全部 `rc=0` 而 `n_info` 不变（8→8），此版本下它不是删除语义；
文本层按 key 丢弃不受该语义影响，且 native 的 GATK 兼容文本本来就由这段代码自己组装。

## 第 80 轮：文本层过滤**编译通过并运行，但仍无效果**——需要确认这些行的文本到底由谁生成

按第 79 轮的三处改法实施（`gatk_compatible_record_text` 加 `drop_read_level` 默认参数、
INFO 组装循环内按 key 跳过四个键、两处调用点传入 `materialized_spanning_locus ||
finalized_monomorphic_ref`），**编译通过**，整段复验：残差仍 **84**（C 类 9 行、E 类 2 行原样）。

⇒ 说明这些行的 GATK 兼容文本**不是**由 `gatk_compatible_record_text()` 的 INFO 组装分支产出的
（或该函数在这一步拿到的 `computed.record` 已是 moved-from 状态、标志读不到）。
改动已回退。

**下一步（一次插桩即可判定，插桩位置作用域安全）**：在 `gatk_compatible_record_text()` 入口打印
`formatted` 的 POS 与第 8 列（INFO），并在 `format_gatk_float_value()` 里打印被调用的 key；
用整段命令跑、只 grep `pos=10008964`：
- 若该函数**根本没被调用**（无对应输出）⇒ 这些行走的是另一条写出路径（如 `bcf_write` 直写），
  过滤要挂到那条路径上；
- 若被调用但 key 列表里没有那四个键 ⇒ 它们是在**别处**（例如 `apply_gatk_annotation_compatibility()`
  之后、文本组装之前）被写入的。

**本轮结论**：C/E 类的问题已从"来源不明"推进到"**已排除两条修复路径（htslib 删除、文本层 INFO 过滤）**，
并明确了下一次插桩要回答的一个二选一问题"。这是可继续的具体状态，而不是卡住。

## 第 81 轮：文本函数**确实看到**这些键；但两种传标志的写法都无效，剩下的写法已明确

**取证（整段命令 + 只 grep 10008964）**：在 `gatk_compatible_record_text()` 入口打印第 8 列：

```
[TXT] INFO=BaseQRankSum=1.026;DP=63;ExcessHet=0;MLEAC=.;MLEAF=.;MQRankSum=1.844;ReadPosRankSum=0.666
（输出行： BaseQRankSum=1.03;DP=63;ExcessHet=0.0000;MLEAC=.;MLEAF=.;MQRankSum=1.84;ReadPosRankSum=0.666）
```

⇒ 该函数**确实**处理这一行、且 `fields[7]` 里**确实**带着这四个键（所以第 80 轮"文本层过滤"
的思路本身是对的）。

**两种写法实测均无效**（整段残差都维持 84）：

| 写法 | 轮次 | 结果 |
| --- | --- | --- |
| 传入 `computed.record.materialized_spanning_locus \|\| finalized_monomorphic_ref` | 80 | ✗ |
| 传入 `encoded.record.…`（move 之后的对象） | 81 | ✗ |

**结论**：调用点取到的标志在那一刻不可靠（move-from / 标志未随对象传递）。
**唯一还没试过、且证据最充分的写法**：在**编码 lambda 的入口**（`std::move` 之前）先把标志
捕获进一个局部布尔量，再把该局部量传给文本函数——第 77 轮的插桩正是打在那一行，
当时读到 `span=1 mono=1`，**证明该处标志为真**：

```cpp
[&](GenotypeComputed computed) -> std::optional<GenotypeEncoded> {
    GenotypeEncoded encoded;
    const bool drop_read_level = computed.record.materialized_spanning_locus ||
                                 computed.record.finalized_monomorphic_ref;   // ← 新增（move 之前）
    encoded.record = std::move(computed.record);
    ...
    encoded.text = gatk_compatible_record_text(<formatted>, drop_read_level);    // ← 用局部量
```
配合 `gatk_compatible_record_text(..., bool drop_read_level = false)` 与 INFO 组装循环内的
key 跳过（这两处第 80 轮已写好并编译通过，可直接复用）。预期整段残差 **84 → 73**。
（本轮最后一次尝试因锚点已被上一轮改掉而未写入，改动已回退，树保持已验证状态。）

## 第 82 轮：`drop_read_level` 在该行上为假——三种传参写法全部无效，下一步探针已定

按第 81 轮"唯一剩余写法"实施：在**编码 lambda 入口**（`std::move` 之前）捕获
`drop_read_level = materialized_spanning_locus || finalized_monomorphic_ref`，
再传给 `gatk_compatible_record_text(..., drop_read_level)`（函数加默认参数 + INFO 组装循环内按 key 跳过，
两处第 80 轮已编译通过）。编译通过，整段复验：**残差仍 84**（C 类 9 行、E 类 2 行原样）。

三种写法全部无效：

| 写法 | 轮次 | 结果 |
| --- | --- | --- |
| 调用点读 `computed.record.<flags>` | 80 | ✗ |
| 调用点读 `encoded.record.<flags>`（move 后） | 81 | ✗ |
| 编码 lambda 入口捕获局部量（move 前，第 77 轮插桩证实该处为真） | 82 | ✗ |

**推论**：那些 C/E 行的文本**不是**经由本次改动的那个调用点/那条组装循环产出的
（尽管第 81 轮的 TXT 探针在 `-L 20:10000000-10009999` 下确实看到该行带着四个键）。

**下一步（决定性的一次插桩，1 次运行即可）**：把探针放到 `gatk_compatible_record_text()`
**内部**，同时打印
① `drop_read_level` 的取值、② 该行的 POS、③ 组装前 `fields[7]` 与组装后 `info_text.str()`
的键列表；用整段命令跑、只 grep `pos=10008964`。三种可能结果各自直接指向修法：
- `drop_read_level=false` ⇒ 该行走的不是编码 lambda 的那两个调用点；
- `=true` 而组装后键仍在 ⇒ 跳过逻辑没落在实际写键的那一层（可能有两处组装循环）；
- 组装后键已消失而输出仍有 ⇒ 文本在之后被二次处理（回到"谁在最后写这四个键"的问题）。

改动已回退；树保持已验证状态。

## 第 83 轮：**两次关键发现**——文本路径确为写键之处，且"键"与"值"是分开写的

用一个**无条件的**过滤做判别实验（不改签名、不改调用点，只在 INFO 组装循环里按 key 跳过）：

```
row at 10008964 (INFO):
BaseQRankSum;DP=63;ExcessHet=0.0000;MLEAC=.;MLEAF=.;MQRankSum;ReadPosRankSum
残差：84 → 336（变差）
```

**发现 1（路径确认）**：输出**确实**被这次改动影响 ⇒ 这四个键就是在这条 INFO 组装循环里写出的
（此前三次"传标志"无效不是路径问题）。

**发现 2（循环结构，修法的关键）**：被跳过后输出里**键仍在、只是没有值**
（`BaseQRankSum;` 而非 `BaseQRankSum=1.03;`）⇒ 该循环**先把 key 写出去，再单独写 `'=' + 值`**
（我的 `continue` 放在写值那一行之前，所以只丢掉了值）。因此**跳过必须放在循环体最前面**
（写 key 之前），放在写值那一行之前是错的。

**下一步（两步，均已具体化）**：
1. 把跳过逻辑移到循环体**第一行之前**（在任何 `info_text <<` 之前）；
2. 用第 82 轮设计的内部探针（打印 `drop_read_level` 取值 + POS + 组装前后键列表，1 次运行）
   确认标志在该处为真——三次传参写法无效说明**标志在那一刻是假**，
   而第 77 轮在编码 lambda 入口读到的是真，二者矛盾需用该探针一次定论。

本轮改动已回退（无条件过滤会让残差变差到 336，不能保留）。

## 第 84 轮：C/E 类修复的**完整配方**（放置位置已证，键清单已补全）

两个实验（都用 `FASTGATK_DROP_RL` 环境变量做无条件过滤，只为定位机制，不作为实现）：

**实验 1——跳过放在写值行之前（第 83 轮）**：输出变成 `BaseQRankSum;DP=63;...`（键在、值没了），
残差 84 → 336 ⇒ 循环结构是「先写 key（`info_text << ordered[index].key;`），再写 `'=' + 值`」。

**实验 2——改为在循环**之前**过滤 `ordered` 向量（`std::remove_if` 擦除）**：
- 目标行变为 `DP=63;ExcessHet=0.0000;MLEAC=.;MLEAF=.`（四个键**连键带值**干净移除，无分隔符残留）；
- 而 GATK 的同一行是 `DP=63;MLEAC=.;MLEAF=.` ⇒ **还差一个键：`ExcessHet`**（native 写、GATK 不写）。

### 完整配方（下一轮可直接实现）

1. **放置位置**：在 INFO 组装循环**之前**对 `ordered` 做 `std::remove_if` 擦除
   （循环头 `:4983`，键与值由同一项输出，擦除即可同时去掉键与值，且不会留下多余的 `;`）；
2. **键清单（C/E 类需丢弃）**：`MQ`、`BaseQRankSum`、`MQRankSum`、`ReadPosRankSum`、**`ExcessHet`**；
   另有 FORMAT 侧的 `PGT/PID/PS`（E 类）；
3. **条件**：仍需把「该行是物化行 / REF-only 物化行」的状态传进 `gatk_compatible_record_text()`
   ——第 80/81/82 轮三种传参写法实测均无效，说明标志在文本阶段取不到真值，
   故下一轮先用内部探针（打印传入的布尔量与 POS，1 次运行）确认**哪里能读到真值**，
   再按上述配方落地；预期 C/E 共 11 行消失，整段残差 **84 → 73**。

> 说明：无条件过滤会让其他有这些注释的行变差（残差 336），故不能作为实现，仅用于定位机制。

## 第 85 轮：两类物化行的**键集不同**——一个布尔量不够，必须分开

按第 84 轮配方实施（编码 lambda 入口捕获标志 + 文本层擦除五个键，**本次确认构建成功**、
二进制时间戳已更新），整段复验：

| 类 | 修前 | 修后 |
| --- | --- | --- |
| `*`/`*`（跨位点物化行） | 40 | **109**（变差） |
| `.`/`.`（REF-only 物化行，C/E） | 11 | 3+3+2（拆成三小类，键已基本对齐） |
| 总计 | 84 | **150**（变差） |

**结论**：两类物化行的**键集不同**，用一个布尔量（跨位点行 ∨ 单态参考行）导致误伤：

| 行族 | GATK 丢弃的键 | 证据 |
| --- | --- | --- |
| 跨位点物化行（`*` 行） | `MQ` + 三个 rank sum（**保留 `ExcessHet`**） | 第 62 轮已按此修好（40 行只剩 AC/AN 问题）；本轮把 `ExcessHet` 一并丢掉后该类从 40 涨到 109 |
| REF-only 物化行（C 类） | `MQ` + 三个 rank sum **+ `ExcessHet`** | 第 84 轮实验 2：擦除后目标行仅剩 `ExcessHet=0.0000` 与 GATK 不同 |

**最终改法（下一轮，两块布尔量或一个枚举）**：
1. 文本函数接收两个独立开关：`drop_read_level`（两类共有）与 `drop_excess_het`（**仅 REF-only 物化行**）；
2. 擦除逻辑仍放在 INFO 组装循环**之前**（第 84 轮已证：键与值同时移除、无分隔符残留）；
3. 两个开关在**编码 lambda 入口**（`std::move` 之前）分别由
   `materialized_spanning_locus` 与 `finalized_monomorphic_ref` 取值；
4. REF-only 行另需丢 FORMAT 的 `PGT/PID/PS`（E 类）。
预期：`*`/`*` 保持 40（不再误伤）、C/E 的 11 行消失 ⇒ 整段残差 **84 → 73**。

## 第 87 轮：把 `drop_read_level` 收窄到 span 行的尝试**无任何效果**；判别式仍需细化

第 86 轮后余下的 REF-only 残差共 13 行，逐行看清后分三类：

| 类 | 行数 | 例 | GATK | native |
| --- | --- | --- | --- | --- |
| 样本列（GT/GQ） | 3 | 10008964 | `0/0:34:99` | `./.:34:0` |
| QUAL + 样本列 | 3 | 10024301 | `163.67` | `Infinity` |
| INFO（GATK **保留**读级注释） | 4 | 10041698 | `BaseQRankSum=1.89;DP=80;ExcessHet=0.00;…;MQ=58.63;…` | `DP=80;MLEAC=.;MLEAF=.`（被第 86 轮误删） |
| INFO + FORMAT + 样本列 | 2 | 10077008 | `GT:DP:RGQ` + 读级注释 | `GT:DP:RGQ:PGT:PID:PS` + 注释被删 |

据此把 `drop_read_level` 从 `span || mono` 收窄为 `span`（理由：GATK 在"真实记录退化为 REF-only"的行上
**保留**读级注释，如 10041698，而在"被跨接记录覆盖的物化行"上丢弃，如 10008964）。
**实测：行为完全不变**（两种构建逐位点相同：99912/99993 个位点一致；两个探针位点的
INFO/整行匹配结果也一致），残差都维持 81。

**推论**：10041698 那类行的 `materialized_spanning_locus` **也是真**（否则收窄就会恢复其注释）。
最可能的原因：该坐标上**同时**存在「真实记录起点」与「由其它跨接记录物化出来的行」，
两者在分组阶段被并到同一条记录上，标志由物化的那一条带来。因此判别式应当细化到
**"该坐标是否存在输入记录起点"**（而不是"该记录是否被物化"）：
- 有起点 ⇒ 保留读级注释（GATK 用起点记录的注释）；
- 无起点（纯物化）⇒ 丢弃。

**下一步**：在 `materialize_spanning_loci()` 已有的 `occupied` 判据外，额外记一个
`Record::synthetic_locus_without_start`（仅在确实没有起点时置位），并用它驱动文本层开关；
预期把上表第 3、4 类（共 6 行）恢复为与 GATK 一致，残差 **81 → 75**（其余 6 行的 GT/GQ 与
QUAL=Infinity 属另外两类问题）。
（本轮改动已回退，树保持第 86 轮已验证状态。）

## 第 88 轮：与我先前的判断**相反**——那 6 行是第 86 轮**引入**的回归，第 87 轮的"无效果"也需重判

**证据链（本轮实测）**：

| 构建 | 10041698 的 INFO | 说明 |
| --- | --- | --- |
| `native4`（第 66 轮前后，**修复前**） | `BaseQRankSum=1.89;DP=80;ExcessHet=0.0000;MLEAC=.;MLEAF=.;MQRankSum=0.194;ReadPosRankSum=-1.640e-01` | **有**读级注释（与 GATK 一致） |
| `nativeF`（第 82 轮，跳过错放在写值行前） | `BaseQRankSum;DP=80;…;MQRankSum;ReadPosRankSum` | 键在、值没 |
| `nativeJ` / `nativeK`（第 85/86 轮） | `DP=80;MLEAC=.;MLEAF=.` | 键**被删** |

而编码器插桩显示该位点的标志是 `[FL] pos=10041698 span=0 mono=0`：

```
[FL] pos=10008964 span=1 mono=1 nall=1 alt=.
[FL] pos=10041698 span=0 mono=0 nall=1 alt=.
```

⇒ 该行的两个标志都是 0，**按第 86 轮的条件本不该被删**，但它的键确实没了。
**所以第 86 轮在这一点上引入了回归**（净效果仍是 84 → 81：修好 3 行、破坏 6 行中…
需重新逐行核对净账），且第 87 轮"把条件收窄为 span 后行为不变"的结论也需重判
（两次构建逐位点相同，说明**收窄没生效**或**删除并不来自该条件**）。

**下一轮必须先做的一件事**：用 `native4`（修复前）与当前构建**逐位点**统计「谁被修好、谁被弄坏」，
得到第 86 轮的真实净账；若确有回归，则把 `drop_read_level` 的判据换成
**"该坐标无输入记录起点"**（而不是任何现有标志），或直接回退第 86 轮的 read-level 部分、
只保留 `drop_excess_het`（后者已由第 84 轮的实验 2 单独证实）。

## 第 89 轮：第 86 轮的真实净账（逐位点，用修复前/后输出直接统计）

用 `native4`（第 66 轮前后，**修复前**）与 `nativeK`（第 86 轮后）对 GATK 逐位点比对：

```
FIXED by round86: 3   BROKEN by round86: 0
 sample fixed: 10011519, 10090291, 10093570（均为 REF-only 行）
```

⇒ **没有任何位点由"匹配"变成"不匹配"**，第 86 轮净收益 +3 位点。

**但需要更细的判读**：10041698 这类位点在**修复前就已经因别的列（样本列 GT/GQ）不匹配**，
所以它在位点级别上一直是"不匹配"，第 86 轮把它的 INFO 从"与 GATK 一致"改成"缺键"
（field 级别是退步），位点级别看不出来。也就是说：

- **位点级别**：+3 / −0（第 86 轮是净收益 ✓ 应保留）；
- **字段级别**：在那 6 个本来就不匹配的位点上，读级注释被误删（GATK 保留、native 删了）。

**因此下一步的最优做法**：保留第 86 轮的 `drop_excess_het`（第 84 轮实验 2 已独立证实）与
`drop_read_level`，但把 `drop_read_level` 的判据从"`span || mono`"改成
**"该坐标没有输入记录起点"**（第 88 轮的插桩显示 10008964 是 span=1/mono=1 的物化行，
而 10041698 是 span=0/mono=0 —— 说明**现有标志无法区分二者**，需要一个新标志）。
一旦该判据到位，那 6 个位点的 INFO 会恢复与 GATK 一致，届时若 GT/GQ 问题也修好，
它们将从"不匹配"变成"匹配"（REF-only 残差 13 → 7）。

## 第 90 轮：代码核对后矛盾仍在——唯一没做过的那个插桩就是下一步

核对已提交（第 86 轮）的代码，逻辑与设计一致：

```cpp
if (drop_read_level || drop_excess_het) {
    ordered.erase(std::remove_if(..., [&](const auto& entry) {
        if (drop_read_level && (entry.key == "MQ" || entry.key == "BaseQRankSum" ||
                               entry.key == "MQRankSum" || entry.key == "ReadPosRankSum")) return true;
        return drop_excess_het && entry.key == "ExcessHet";
    }), ordered.end());
}
```
两处开关分别在**编码 lambda 入口**（`:6732`/`:7503`）由
`materialized_spanning_locus` 与 `finalized_monomorphic_ref` 取值，并传给文本函数（`:6746`/`:7517`）。

而 10041698 在该 lambda 上的标志实测为 `span=0 mono=0`（第 88 轮插桩）：
**按此代码它不该被删除，但它的键在修复后确实不见了**（第 88 轮的三版对照）。

⇒ 矛盾未解。已知的排除项：
- 不是调用点/传参写法（三种都试过）；
- 不是 erase 的条件写错（本轮逐行核对）；
- 不是该 lambda 取错对象（插桩就打在同一个 lambda 里、`std::move` 之前）。

**唯一从未做过的插桩（下一步就做这一个）**：把探针放进 `gatk_compatible_record_text()`
**内部**，对目标 POS 打印 ① 两个开关的取值、② 擦除**前**的键列表、③ 擦除**后**的键列表、
④ 函数入口收到的第 8 列。一次整段运行即可给出四种可能中的一种：
- 开关为假却键消失 ⇒ 键是在**调用本函数之前**就没了的（回到计算阶段，需查 `ordered` 的来源）；
- 开关为真 ⇒ 说明同一位点上存在**两条**记录（一条 `span=1` 一条 `span=0`），
  编码器插桩因 `sort -u` 或同一 POS 只打印了其中一条 —— 那就按 POS 分别打印每条记录的标志；
- 键在擦除后仍在 ⇒ `remove_if` 的谓词没匹配到（键名有前后缀）；
- 键在入口就已不在 ⇒ 与第一项相同。

第 86 轮的修复**保留**（位点级 +3/−0，无回归）；本轮的核对与结论已归档，改动为 0。

## 第 91 轮：函数内部插桩给出决定性结论——**两个开关对两族都是 1**，现有标志无法区分

在 `gatk_compatible_record_text()` 内部插桩（打印两个开关取值、擦除前后键列表、入口第 8 列）：

```
[TXT2] pos=10008964 drl=1 deh=1 before: BaseQRankSum DP ExcessHet MLEAC MLEAF MQRankSum ReadPosRankSum|INFO=BaseQRankSum=1.026;DP=63;…
[TXT2] after: DP MLEAC MLEAF
[TXT2] pos=10041698 drl=1 deh=1 before: BaseQRankSum DP ExcessHet MLEAC MLEAF MQRankSum ReadPosRankSum|INFO=BaseQRankSum=1.887;DP=80;…
[TXT2] after: DP MLEAC MLEAF
```

**结论**：两族的 `drop_read_level`/`drop_excess_het` **都是 1**，擦除也都生效（键列表都从 7 项变 3 项）。
也就是说：

- 第 86 轮的条件对 10041698 这类行**判成了"应删"**（GATK 实际保留）⇒ 规则的**判据本身**不成立；
- 第 88 轮 `[FL]` 探针在该位点读到 `span=0 mono=0`，与这里 `drl=1`（即 `span||mono` 为真）**矛盾**
  ⇒ 那次的读数不可靠（很可能同一位点上有两条记录，或读到的是另一条实例）；
- 第 87 轮"把条件收窄为 span 后行为不变"也因此**不能作为判据有效的证据**。

**下一步（两条，按顺序）**：
1. 先用一次插桩确认 `[FL]` 矛盾：在**编码 lambda 入口**同时打印
   `pos`、`materialized_spanning_locus`、`finalized_monomorphic_ref`、**`record.value->pos`**、
   `alleles[0..1]`，并**不要用 `sort -u`**（避免去重掩盖同 POS 的多条记录）；
2. 用新判据替换：在 `materialize_spanning_loci()` 内合成记录时置一个**独立**标志
   （例如 `synthetic_no_start`），并把文本层的 `drop_read_level/drop_excess_het` 改为由它驱动；
   预期 10041698 这类行不再被删（INFO 恢复与 GATK 一致），REF-only 残差 13 → 7。

## 第 92 轮：矛盾解开 + 规则**完全确定**（两个开关的判据不同）

**矛盾解开（0-based vs 1-based 的坑）**：改用 0-based 位置重新插桩：

```
[FL2] pos0=10008963 bcpos0=10008963 span=1 mono=1 block=0 nall=1   ← 文本 POS 10008964
[FL2] pos0=10041697 bcpos0=10041697 span=0 mono=1 block=0 nall=1   ← 文本 POS 10041698
```

⇒ 第 88 轮那次 `[FL] span=0 mono=0` 是**匹配到了另一条记录**（`record.pos` 是 0-based，
我却拿 1-based 的文本 POS 去比），那次读数作废。真实情况是：
**10041698 是 `span=0 mono=1`**（真实记录退化为 REF-only），**10008964 是 `span=1 mono=1`**（物化覆盖位点）。

**据此把两个开关的判据完全确定**（与三处实测全部吻合）：

| 观测 | GATK 行为 | 结论 |
| --- | --- | --- |
| 10008964（span=1 mono=1） | INFO 只有 `DP;MLEAC=.;MLEAF=.`（无 rank sum、无 ExcessHet） | 两者都丢 |
| 10041698（span=0 mono=1） | 保留 rank sum/MQ **且保留 `ExcessHet=0.00`** | 两者都留 |
| `*` 跨位点行（span=1 mono=0，第 85 轮） | 保留 `ExcessHet` | 读级注释丢、ExcessHet 留 |

⇒ **`drop_read_level = materialized_spanning_locus`**；
⇒ **`drop_excess_het = materialized_spanning_locus && finalized_monomorphic_ref`**。

**实现要点（下一轮，两行赋值）**：把 `:6732` 与 `:7503` 两处

```cpp
const bool drop_read_level = computed.record.materialized_spanning_locus;
const bool drop_excess_het = computed.record.materialized_spanning_locus &&
                             computed.record.finalized_monomorphic_ref;
```

写进去（第 86 轮的文本函数与擦除逻辑不用动）。预期：`*`/`*` 回到 40（不再误伤）、
REF-only 的 INFO 类（10041698/10098308/10099270/10077008/10077010 五行）INFO 恢复一致，
残差 **81 → 76 左右**（其余为 GT/GQ 与 `QUAL=Infinity` 两类）。

> 本轮两次尝试应用该规则时**都因锚点文本已变而未写入**（第一次把两开关都设成 span，
> 实测 `*` 行涨到 109；第二次的 old 文本不匹配），随后 `git checkout` 回退到第 86 轮已验证状态。

## 第 94 轮：`ExcessHet` 的零值渲染规则**被证否为独立可修项**，并暴露出前置条件

第 93 轮定位到 5 个位点只差 `ExcessHet=0.00`（GATK）vs `0.0000`（native），本轮按
「精确 0 → 2 位小数、否则 4 位」实现并整段复验：

```
differ 81 → 402（其中 252 行是 concrete 行）
positions with matching INFO: 99918 → 99597
```

**证否 1（passthrough 假设）**：那 5 行对应输入记录的 INFO 里 `ExcessHet=0.0000`（不是 `0.00`），
故 GATK 并非原样透传输入文本 ⇒ 该解释作废。

**证否 2（独立可修）**：把「精确 0 → `0.00`」规则全量应用后，**252 个具体变异行**由匹配变为不匹配。
逐条看：这些行的 GATK 值是**微小的非零**（渲染成 `0.0000`），而 native 的值是**精确 0**
（在旧规则下渲染成 `0.0000`，恰好与 GATK 一致；新规则下渲染成 `0.00` ✗）。

⇒ **结论**：`ExcessHet` 的零值渲染规则本身很可能是对的（5 行支持它），但它**不能单独落地**——
前置条件是 **native 的 `ExcessHet` 数值必须与 GATK 一致**：现在 native 在 252 个位点上算出精确 0，
而 GATK 算出一个微小非零值（属"引擎数值层面"的差异，与本会话早前记录的 1e-16 级差异同类）。

**因此下一轮的目标应改为**：先对齐 `ExcessHet` 的**数值**（例：在具体变异行上打印 native 与 GATK 的
原始 double 值，定位是算法/求和顺序还是精度问题），数值对齐之后再启用零值渲染规则，
届时 5 行会被修好且 252 行不受影响。本轮改动已回退（净负收益）。

## 第 96 轮：`PGT/PID/PS` 清除**机制有效**，但判别式需再细化（两种尝试都净负）

目标：第 95 轮清单第 6 项（REF-only 行多带 `PGT/PID/PS`，例 `20:10077008/10077010`）。
在编码阶段对物化行清除这三个 FORMAT 键（`bcf_update_format_int32(..., nullptr, 0)`）：

| 判别式 | 结果 |
| --- | --- |
| `span \|\| mono` | `20:10077008/10077010` 的 FORMAT **确实变成 `GT:DP:RGQ`**（与 GATK 一致 ✓），但整段 differ 81 → **90**（**9 行由匹配变不匹配**） |
| 仅 `span` | differ 81 → **90**（star 行的 FORMAT 也被误清：`('star','star',(8,))` 9 行、以及 20+20 行形状改变） |

**结论**：
1. **机制有效** ✓（FORMAT 清除能在文本之前生效，与第 78 轮 INFO 删除语义的坑不同）；
2. **判别式仍不对**：GATK 在多数 `*` 物化行上**保留** `PGT/PID/PS`，只在「物化且退化为 REF-only」的行上丢弃
   —— 与 `ExcessHet` 的规律一致（第 93 轮），**下一轮应直接用 `span && mono`**（即
   `materialized_spanning_locus && finalized_monomorphic_ref`）再测一次；
3. 需要顺带确认 `10077008` 的两个标志取值（本轮未插桩；若它是 `span=0 mono=1`，则 `span && mono` 不成立，
   判别式还要再拆一层）。

两次尝试均已回退（净负收益）；树保持第 93 轮已验证状态。

## 第 98 轮：把整段残差比对固化成**一条命令**（本会话的测量手段此前都是临时命令）

新增 `fastgatk-native/scripts/measure_dense_residual.py`：对同一输入、同一参考、同一 `-L`
窗口跑 pinned GATK 4.6.2.0 与 native，输出

- 两侧行数/位点数、仅 GATK / 仅 native / 共有但不同的位点数；
- 差异按（GATK ALT 形状、native ALT 形状、不同的前导列）分组计数；
- `--positions` 指定的位点逐列 dump（便于单点排查）；
- 末尾一行 JSON（可直接归档）。

也支持 `--gatk-vcf/--native-vcf` 直接比较已有输出（不必重跑工具），
`--require-identical` 可当作"零残差"断言使用（默认不退出码非零——它是测量工具，不是门禁）。

**已实测复现当前状态**（与本文档此前的数字一致）：

```
# GATK rows=100000 positions=100000
# native rows=99993 positions=99993
# only GATK=7 only native=0 differ=81
GATK shape   native shape differing columns        count
star         star         (7, 9)                   38
star         refonly      (4, 5, 6, 7, 8, 9)       30
refonly      refonly      (7,)                     5
refonly      refonly      (9,)                     3
refonly      refonly      (5, 9)                   3
star         star         (7,)                     2
```

顺带把分类细化了一档：此前文档里把 `*`/`*` 记成"40 行"，现在看清是
**38 行同时差 INFO 与样本列、2 行只差 INFO**；REF-only 类也从"11 行"细化为 5+3+3。

## 第 99 轮：用新工具把两个最小类**逐列看清**（各自只差一两处，且发现 native 缺 MQ）

用 `measure_dense_residual.py` 的分类 + 逐列 dump 定位到两个最小类：

### 类 `refonly/refonly (9,)`——3 行（10008964 / 10008965 / 10076991）

```
GATK   N  .  Infinity  .  DP=63;MLEAC=.;MLEAF=.  GT:DP:RGQ  0/0:34:99
NATIVE N  .  Infinity  .  DP=63;MLEAC=.;MLEAF=.  GT:DP:RGQ  ./.:34:0
```

⇒ INFO 完全相同，只有**样本列**不同：GATK 是**纯合参考调用**（`0/0`、RGQ 99），native 是 no-call
（`./.`、RGQ 0）。修法候选：REF-only 物化路径在「参考块/未覆盖」分支里把 GT 赋成 missing、
GQ/RGQ 清零；GATK 在这些位点保留 hom-ref 调用与参考置信度质量。需要先插桩确认走的是哪个分支
（`materialize_gatk_monomorphic_ref_call` 的 `uncovered` 分支会清 GT/GQ，但它同时会清 INFO 的
`DP`，而这里 `DP=63` 仍在 ⇒ 疑点：不是那条分支）。

### 类 `refonly/refonly (7,)`——3 行（10041698 / 10098308 / 10099270）

```
GATK   ... ExcessHet=0.00   ... ;MQ=58.63; ...
NATIVE ... ExcessHet=0.0000 ... (无 MQ)
```

⇒ 两处差异：
1. `ExcessHet` 的零值渲染（`0.00` vs `0.0000`）——第 94 轮已证**不能单独落地**（前置：数值对齐）；
2. **native 根本没有写 `MQ`**（GATK 从 `RAW_MQandDP` 推导出 `MQ=58.63`/`MQ=60.00`）——
   这是本轮**新发现**的缺口：这些 REF-only 行上 native 缺 `MQ` 注解（不是被误删，
   因为 compute 阶段的抑制器只对 `span=1` 的行生效，而这些行 `span=0`）。

**下一步（两个最小靶子，各自只差一处）**：
- 3 行的 `(9,)`：先插桩确认 REF-only 物化的分支，再把「纯合参考 + RGQ」补齐；
- 3 行的 `(7,)`：补 `MQ`（从 `RAW_MQandDP` 推导，与具体变异行的算法一致）——这一项与
  `ExcessHet` 渲染无关，可独立修好并立即减少 3 行残差。

## 第 100 轮：`MQ` 缺失的**真正原因**（不是被删，而是早退路径不计算它）

按第 99 轮的"最小靶子"补 MQ：在编码阶段对 `allele_count == 1` 的行，从 **bcf 的
`RAW_MQandDP`** 推导 `MQ = sqrt(sumMQ2/DP)`（公式已由两点验证：`275031,80 → 58.63`、
`7200,2 → 60.00`）。整段复验：**残差仍 81、类别不变**；单点核对 `20:10041698` 显示
**native 仍无 `MQ`** ⇒ 该推导没有生效。

**原因**：`RAW_MQandDP` 是**输入侧**的 INFO，native 的**输出 header 并不声明它**，
故 `bcf_get_info_int32(..., "RAW_MQandDP")` 必然取不到 → 我的推导无从下笔（逻辑本身没错，取数取错了层）。

**因此 MQ 缺失的真正机制是**：GATK 的 REF-only 行仍带 `MQ`（`58.63`/`60.00`），而 native 的这些行
走的是 `finalized_monomorphic_ref` **早退分支**——该分支只调用 `apply_gatk_annotation_compatibility()`，
**跳过了 `update_gatk_standard_annotations()`**（具体变异行的 MQ 就是在那里算出来的）。

**修法（下一轮，明确）**：在 `finalized_monomorphic_ref` 早退分支里补上 MQ 的计算——
取值应来自 **Host 侧的原始 MQ 累加量**（与具体变异行同源，参照 `update_gatk_standard_annotations()`
里 MQ 的算法），而不是从输出 bcf 的 INFO 里找 `RAW_MQandDP`。
本轮改动已回退。

## 第 101 轮：`MQ` 缺失的两个门都找到了，但数据本身在 REF-only 行上已丢失

第 100 轮把 MQ 的算法读清了（`update_gatk_standard_annotations()` 里，从 **输出 bcf 的
`RAW_MQandDP`** 取 `[sumMQ2, DP]` 算 `sqrt(sumMQ2/DP)` ✓ 与 GATK 一致）。本轮据此把整段 MQ 代码
**移到 `if (record.allele_count < 2 || record.gt.empty() || record.ploidy <= 0) return;` 之前**
（这就是"REF-only 行没有 MQ"的第一道门），编译通过、整段复验：**`20:10041698` 仍无 MQ**，
残差仍 81。

⇒ 该行在**输出阶段**的 INFO 里**已经不含 `RAW_MQandDP`**（`bcf_get_info_int32` 取不到），
所以第二道门是**数据本身在 REF-only 行上已被丢弃**（`RAW_MQandDP` 是输入侧注解，
REF-only 物化路径没有把它带下去）。

**完整修法（下一轮，三小步，均已定位）**：
1. 在 `Record` 上加两个 int32（例如 `raw_mq_sum_of_squares`、`raw_mq_depth`），
   在**解码阶段**从输入记录的 INFO `RAW_MQandDP` 读入并随记录传递（此时数据还在）；
2. 在 `update_gatk_standard_annotations()` 的 MQ 分支里，优先用这两个 Host 字段
   （而不是输出 bcf 的 INFO）；
3. 保持第 101 轮已验证的"MQ 计算放在 `allele_count < 2` 守卫之前"这一位置（本轮已证该位置正确，
   只是数据缺失）。预期 `(7,)` 类的 3 行恢复 `MQ`，与 GATK 只差 `ExcessHet` 渲染（第 94 轮的前置项）。

本轮改动已回退（无效果），树保持第 97 轮已验证状态。

## 第 102 轮：MQ 缺失修好（第 43 个已修 bug）——在 `RAW_MQandDP` 被移除**之前**推导

第 101 轮定位到「REF-only 行没有 MQ」有两道门：① `update_gatk_standard_annotations()` 的
`allele_count < 2` 守卫；② 数据本身已丢。本轮找到第二道门的确切位置：

```
apply_gatk_annotation_compatibility()（:4710 起）
    ...
    (void)bcf_update_info_int32(output_header, record.value, "RAW_MQandDP", nullptr, 0);   // :4783 把它删掉
```

而 `finalized_monomorphic_ref` 早退路径**只**调用这个函数（不调用 `update_gatk_standard_annotations()`），
所以 REF-only 行既没被算 MQ、`RAW_MQandDP` 也已被删。

**修复**：在 :4783 删除之前**就地推导 MQ**（这是唯一同时覆盖两条路径、且数据仍在的位置）：

```cpp
if (bcf_hdr_id2int(output_header, BCF_DT_ID, "MQ") >= 0) {
    ... bcf_get_info_int32(..., "RAW_MQandDP", ...);
    if (raw_mq[1] > 0 && raw_mq[0] > 0) {
        const float mq = sqrt(raw_mq[0] / raw_mq[1]);
        bcf_update_info_float(output_header, record.value, "MQ", &mq, 1);
    }
}
```

**验证（整段语料）**：`20:10041698` MQ=58.63、`20:10098308` MQ=60.00、`20:10099270` MQ=56.49
**与 GATK 逐一相同** ✓；位点级残差仍 **81**（这些行另有 `ExcessHet` 零值渲染差异，第 94/99 轮已记录）；
五个基因型门禁全部通过。

> 记账方式的教训：本轮我先用"与上一次输出文件比较"的临时脚本得出 `FIXED 9`，
> 但那个基线文件其实是更早（differ=90）的中间构建 ⇒ 数字虚高。**权威口径要以 GATK 为参照、
> 用 `measure_dense_residual.py` 复算**（本轮复算为 81），临时脚本的比较对象必须确认清楚。

## 第 103 轮：`(9,)` 类（REF-only 行的 GT/GQ）候选位置已收窄，下一步一次插桩即可定

该类 3 行（10008964 / 10008965 / 10076991）只差样本列：

```
GATK   ... GT:DP:RGQ  0/0:34:99
NATIVE ... GT:DP:RGQ  ./.:34:0
```

**判据（用来筛选候选代码位置）**：正确位置必须同时满足
① 把 GT 置为 missing（`./.`）、② 把 RGQ 归零、③ **不清除 INFO 的 `DP`**（这里 `DP=63` 仍在，
所以不是第 53 轮读到的 `uncovered` 分支——那条会 `bcf_update_info_int32(..., "DP", nullptr, 0)`）。

在 `src/genotype_gvcf_tool.cpp` 里与 GT/RGQ 相关的位置共约十处（`bcf_gt_missing` 的 10 处赋值、
`"RGQ"` 的 5 处写入）。**主要嫌疑**是 `:2861` 附近的

```cpp
for (const auto* tag : {"DP", "RGQ"}) { ... }
...
if (bcf_update_format_int32(output_header, record.value, "RGQ", rgq.data(), ...) != 0)
```

即 REF-only 物化路径里同时处理 `DP`/`RGQ` 的那段（它符合判据①③，且 `RGQ` 被写成 0 与观测一致）。

**下一步（一次插桩定案）**：在 `:2861` 所在函数入口打印
`pos / allele_count / finalized_monomorphic_ref / materialized_spanning_locus`，
整段跑一次、grep 这 3 个 POS：命中即确认；未命中再用同样方式把其余候选（`:1386`/`:1457`/`:1753-1801`/
`:1860-1927`）逐个二分（每处一个标记，一次运行可见全部）。确认后按 GATK 的
`0/0` + RGQ=99 补齐，预期该类 3 行消失（残差 81 → 78）。

### 第 103 轮续：`(9,)` 类根因已定位到函数与判据（输入侧证据）

`(9,)` 3 行（10008964 / 10008965 / 10076991）的正确位置已确认是
`materialize_gatk_monomorphic_ref_call()`（`src/genotype_gvcf_tool.cpp:2804`），关键代码：

```cpp
std::vector<int32_t> rgq(sample_count, bcf_int32_missing);   // :2821
if (record.source_has_gq)                                     // :2829
    bcf_get_format_int32(output_header, record.value, "GQ", &source_gq, &source_gq_count);
...
const auto gq = source_gq != nullptr && sample < source_gq_count
    ? source_gq[sample] : bcf_int32_missing;                  // :2836
if (depth <= 0 || gq == bcf_int32_missing || gq == vector_end) continue;   // :2838
if (gq > 0)                                                   // :2840  <-- 判据
    for (int copy = 0; copy < ploidy; ++copy) gt[...] = bcf_gt_unphased(0);
rgq[sample] = gq;                                             // :2844
dp[sample] = depth;                                           // :2845
any_rgq = true;                                               // :2846
```

**观测到的输出形状与这段代码自洽且唯一**：`./.:34:0` 要求 `depth > 0`（34）、`gq == 0`
（故 `:2838` 的 `continue` 未命中、`rgq` 被写成 0）、但 `gq > 0` 为假（故 GT 保持
`bcf_gt_missing` → `./.`）。即 **走到了 `:2844-2846`，却在 `:2840` 处不满足**——只有
`gq == 0` 能同时产生这三个观测。

**输入侧（GATK 自己的 `expected.testGVCFMode.gatk4.g.vcf`）在同一坐标的实测**：

| 位置 | 覆盖它的输入记录 | 该记录 sample |
|---|---|---|
| 10008964/10008965 | `POS=10008952 REF=CACACACACACACA ALT=C,CCA,CCACACACACA,CCACACACACACA,<NON_REF>`（REF 跨 10008952–10008965，**下一个记录从 10008966 才开始**） | `GT=2/4, AD=1,3,14,5,11,0, DP=34, GQ=99, SB=1,0,25,8` |
| 10076991 | `POS=10076989 REF=CAA ALT=C,CA,CAAA,<NON_REF>`（REF 跨 10076989–10076991，下一个记录从 10076992 开始） | `GT=2/3, AD=7,4,13,7,0, DP=41, GQ=31, SB=4,3,12,12` |

两条独立结论：
1. 这两个坐标**正是"落在某条记录 REF 跨度内、但不是记录起点"的位置**，与
   `materialize_spanning_loci()` 的物化条件一致，`record.materialized_spanning_locus` 应为真。
2. 10008952 那条记录的 `DP=34`、`GQ=99` 与 GATK 输出该行的 `0/0:34:99` **逐字段相等**——
   说明 GATK 对这类物化行的规则是 **继承来源记录的 DP 与 GQ，并把 GT 写成 hom-ref `0/0`**；
   而 native 的 `DP=34` 也对上了（`:2845` 成功），只有 `GQ` 读成了 **0**、GT 因此没转 hom-ref。

**因此根因只剩两个互斥可能，一次插桩即可判定**：
- (A) `:2830` 的 `bcf_get_format_int32(..., record.value, "GQ", ...)` 读的是**已被前序阶段改写过的
  输出记录**（GQ 此时已是 0），而 `source_has_gq` 为真——于是 `gq == 0`，完全复现观测。
  修法：在 GQ 被清除前缓存来源 GQ（或改读来源向量的 `record.gq`）。
- (B) `source_gq` 读到了正确值但该 sample 的 GQ 本来就是 0；若如此则是 `:2840` 的 `gq > 0`
  判据与 GATK 不一致（`:2824` 的注释已指出 GATK 的依据是 **GQ 的存在性**而非 `> 0`），
  修法即把 `gq > 0` 改为 `gq != missing`。

**插桩方案（一次运行定案）**：在 `:2832` 循环内、`:2838` 之前按 POS 打印
`pos / allele_count / materialized_spanning_locus / finalized_monomorphic_ref / source_has_gq /
source_gq_count / source_gq[sample] / min_dp / dp / gq`，整段跑一次并 grep 这 3 个 POS。
注意 **`record.pos` 是 0-based**（文本 POS 才是 1-based），grep 时别用错。

**每行值不同，需逐行核对**：10008964/10008965 的来源是 `DP=34/GQ=99`，而 10076991 的来源是
`DP=41/GQ=31`；若 GATK 在 3 行上打印的并非同一个 `34/99`，则"继承来源 DP/GQ"这条规则要按行区分，
插桩输出里一并打印 GATK 侧对应行以对拍。

预期收益：该类 3 行消失（残差 81 → 78）；若同一函数也是 `star/refonly (4,5,6,7,8,9)` 那 30 行的
成因，收益更大——插桩时把 `-` 这三行之外再 grep 一条 `star/refonly` 的 POS 一起看。

## 第 104 轮：`(9,)` 类已修（bug 44）——物化跨位点行的样本列读错 GQ 来源

**插桩实测**（`FASTGATK_DEBUG_MONOREF`，`materialize_gatk_monomorphic_ref_call()`）：

```
pos1=10008964 span=1 allele_count=2 source_has_gq=1 source_gq_count=1
             source_gq_val=0 rec_dp=34 depth=34 gq=0 would_homref=0
pos1=10008965 （同上，rec_dp=34）
pos1=10076991 source_gq_val=0 rec_dp=31 depth=31 gq=0 would_homref=0
```

`./.:DP:0` 的成因唯一：`:2838` 未 continue（depth>0）、`rgq` 写成 `gq`(=0)、`:2840` 的
`gq > 0` 为假 ⇒ GT 停在 `bcf_gt_missing`。即**读到的 GQ 是 0**（不是缺失）——假设 (A) 成立：
该函数从**已被计算阶段改写**的 `record.value` 读 GQ，此时它已是 PL 派生值 0。

**GATK 真值逐位点**（整段 dense）：

| 位点 | 来源记录 | 来源样本 | GATK 输出 |
| --- | --- | --- | --- |
| 10008964 / 10008965 | `POS=10008952 REF=CACACACACACACA` | `GT=2/4, DP=34, GQ=99` | `0/0:34:99` |
| 10076991 | `POS=10076989 REF=CAA` | `GT=2/3, DP=31, GQ=41` | `0/0:31:41` |

⇒ 规则：物化行**继承来源 GQ/DP，GT 写 hom-ref**。更正第 103 轮笔记：10076989 的来源是
**DP=31 / GQ=41**（原笔记把 `31:41` 读反了，GATK 输出 `0/0:31:41` 亦可反证）。

**修复**：`Record` 增 `source_gq`，在解码 `read_scalar_format("GQ", destination.gq)` 之后留存
（`source_has_gq` 只保留了有无、没保留值）；物化行处 `gq = record.source_gq[sample]`（仅
`materialized_spanning_locus` 生效）。`gq > 0` 判据**未改**，故非物化路径零影响。

**实测效果**：整段 dense 位点级 **81 → 78**（3 行逐字节一致：`0/0:34:99` / `0/0:34:99` /
`0/0:31:41`），另把原 `(5,9)` 类的第 9 列一并修好（仍差第 5 列）。新残差分类：

| 形态 | 差异列（0-based） | 行数 | 含义 |
| --- | --- | --- | --- |
| `star/star` | (7, 9) | 38 | INFO + 样本列 |
| `star/refonly` | (4,5,6,7,8,9) | 30 | ALT..样本全列（形态分歧） |
| `refonly/refonly` | (7,) | 5 | INFO（ExcessHet 精度） |
| `refonly/refonly` | (5,) | 3 | QUAL |
| `star/star` | (7,) | 2 | INFO |

注意残差地图的列号是 **0-based**：`(7,)`=INFO、`(5,)`=QUAL、`(9,)`=样本列、`(4,)`=ALT。

**新门禁**：`fastgatk-genotype-gvcf-spanning-source-gq-gatk-oracle`（严格；双值对照 99≠41）。

## 第 104 轮的附带收获：7 个 only-GATK 位点已完全刻画（下一轮目标）

用当前树重跑 GATK（`-L 20:10000000-10099999 --include-non-variant-sites`，存
`.diag/round104/gatk-dense.vcf`）后逐位点读出 GATK 的形状——**7 行全部是同一种**：

```
POS=10062936 REF=N ALT=. QUAL=Infinity FILTER=. FORMAT=GT:DP:RGQ SAMPLE=0/0:23:16
   INFO=DP=89;MLEAC=.;MLEAF=.
POS=10087821 REF=N ALT=. QUAL=Infinity FILTER=. FORMAT=GT:DP:RGQ SAMPLE=0/0:53:99
   INFO=DP=103;MLEAC=.;MLEAF=.
POS=10098309 REF=N ALT=. QUAL=Infinity FILTER=. FORMAT=GT:DP:RGQ SAMPLE=0/0:2:6
   INFO=DP=2;MLEAC=.;MLEAF=.
```

三条硬事实：

1. **这 7 个坐标的参考碱基都是 `N`**（实测 `GRCh37.chr20.fa` + `.fai` 逐字节读：
   10062936/37/38、10087821/22、10098309/10 全部 = `N`）。所以 `REF=N` 不是占位符，
   而是真的 N-masked。
2. **它们都有跨接来源记录**（第 103 轮我的快速脚本按 `END=` 起 span，漏掉了
   "无 END 但 REF 很长"的记录，已更正为 `span = max(END, pos+len(REF)-1)`）：

   | 位点 | 来源记录 | 来源样本（GT:AD:DP:GQ:PL:SB） |
   | --- | --- | --- |
   | 10062936-38 | `POS=10062935 REF=CAAA ALT=C,CA,CAA,CAAAA,<NON_REF>`（span 10062935-10062938） | `0/4:...,DP=23,GQ=16` |
   | 10087821-22 | `POS=10087820 REF=CAG ALT=C,CAGAG,<NON_REF>`（span 10087820-10087822） | `0/2:...,DP=53,GQ=99` |
   | 10098309-10 | `POS=10098308 REF=AAT ALT=A,<NON_REF>`（span 10098308-10098310） | `0/0:...,DP=2,GQ=6` |

3. **样本列再次印证第 104 轮刚落地的继承规则**：GATK 的 `23:16` / `53:99` / `2:6`
   分别等于三条来源记录的 `DP:GQ`（site 级 `INFO/DP` 是 89/103/2，**不是**样本 DP）。
   注意这几行是 **monomorphic（`ALT=.`）**，说明"继承来源 DP/GQ"这条规则对 monomorphic
   物化行同样成立。

**形状解读（下一轮的锚点）**：GATK 在这里给的是 monomorphic REF 行——`ALT=.`（不是 `*`）、
`FILTER=.`（未设）、`QUAL=Infinity`、`INFO` 只剩 `DP` 且 `MLEAC=.`/`MLEAF=.`（缺值）。
native 却**整行都没有**。native 已有 `finalized_monomorphic_ref` 这条路径（它正是产出
`ALT=.`、并抑制读级注释的那条），所以下一步不是新写路径，而是查**为什么这 7 行在
N-masked 坐标上走不到/被丢掉**，以及 QUAL 的 `Infinity` 语义：

- 第一嫌疑：`:3939` `base = reference_base(...)` 在 N 位置得到 `"N"`，
  `alleles = "N,*"` ⇒ 合成的是 `*` 行；随后若投影出的基因型不带该删除，`*` 被剪成
  monomorphic 行——需要在剪枝之后确认行是否仍在 `records` 里；
- 第二嫌疑：QUAL 非有限值（`Infinity` / `NaN`）的处理。第 95 轮 backlog 已记下
  3 行 `refonly/refonly (5,)` 是 `163.67` vs `Infinity`（本轮修完样本列后它们**只差 QUAL**，
  见上表 `(5,)` 3 行），与这 7 行的丢行**同根**：native 侧 QUAL 的有限/非有限判定。

**一次插桩定案的方案**：在 `materialize_spanning_loci()` 里对 `position ∈ {10062936,
10087821, 10098309}` 打印 `base / alleles / deletion / projected GT / 是否 push_back`，
再在 compute 出口打印这些坐标的 `allele_count / qual / 是否被 LowQual 丢弃`。

## 第 104 轮附二：78 个残差的**字段级**地图（下一轮的直接靶子）

用当前树跑 dense 双工具并逐字段比对（GATK 100000 位点 / native 99993 位点 / differ 78）：

| 列（0-based） | 含义 | 差异行数 |
| --- | --- | --- |
| 4 | ALT | 30 |
| 5 | QUAL | 33 |
| 6 | FILTER | 30 |
| 7 | INFO | **75** |
| 8 | FORMAT | 30 |
| 9 | 样本列 | 68 |

INFO 内部按键统计：`QD` 70、`AN`/`AF`/`AC` 各 68、`MLEAC`/`MLEAF` 各 63、
`FS`/`SOR`/`DP` 各 44、`ExcessHet` 35。

### 最大一类（38+2 行）：`star/star` 的 `*` 行基因型

同一来源记录 `POS=10008952 REF=CACACACACACACA ALT=C,CCACACACACACA`（REF 跨 10008952-10008965，
输出 GT `1/2`）在它覆盖的 13 个坐标上，GATK/native 的对拍：

| 坐标 | GATK | native |
| --- | --- | --- |
| 10008952（记录自身起点） | `1/2:1,14,11:34:99:1717,539,803,595,0,533` | **逐字节相同** ✓ |
| 10008953-10008963（**11 个**） | `ALT=*`, `QUAL=0`, `LowQual`, `AC=2;AF=1.00;AN=2;MLEAC=2;MLEAF=1.00;QD=-0.00`, 样本 `1/1:1,11:34:62:1184,62,0` | `ALT=*` 但 `AC=0;AF=0.00;AN=0;MLEAC=1;MLEAF=0.500`（**无 QD**），样本 `./.:1,14:34:0:1717,539,803` |
| 10008964-10008965（2 个） | `ALT=.`, `QUAL=Infinity`, `FILTER=.`, `INFO=DP=63;MLEAC=.;MLEAF=.`, 样本 `0/0:34:99` | **已修（第 104 轮）逐字节相同** ✓ |

⇒ 该类的机制是**跨接位点的等位基因投影 + 重新基因分型**：GATK 的
`ReferenceConfidenceVariantContextMerger.replaceWithNoCallsAndDels()` 把"比本记录 REF 短的每个
ALT"都变成 `*`（来源 GT `1/2` 两个等位基因都短于 14bp 的 REF ⇒ 合并后是 `*/*` ⇒ 输出 `1/1`，
PL 由来源 21 个值投影到 `{REFNO_CALL, *}` 得 `1184,62,0`）；native 只按
`called_deletion_allele_index(record)` 取**单个**删除等位基因做投影，来源基因型里的
其它等位基因在投影后的等位基因集合里不存在 ⇒ 退化成 no-call `./.`、
`AC=0;AF=0.00;AN=0`，并且把来源 PL 的对角项 `1717,539,803` 原样带出（不是投影值），
QD 也因此整个缺失。

**未解的边界问题（下一轮先测这个）**：同一跨度内为什么前 11 个坐标是 `*` 行、后 2 个是
monomorphic 行？跨度共 13 个坐标（= len(REF)-1），实测 `*` 只覆盖 11 个
（`REF` 14 - 最长 ALT 12 = 2 恰好是 monomorphic 的个数；但"每个 ALT 的删除覆盖区间取并集"
给的是 13，故现有模型不足以解释）。下一轮建议用同一个来源记录把两个计数都算出来
（`len(REF)`、各 ALT 长度、`called_deletion_allele_index`、以及 GATK 侧 `*` 的实际最大坐标），
确定 `*` 的终止条件是"最长 ALT 长度"还是"某个被保留等位基因的删除长度"。

**第 30 行 `star/refonly (4,5,6,7,8,9)`** 与上面同源（ALT/QUAL/FILTER/INFO/FORMAT/样本全差），
预计同一处修复会一起收掉——即"投影等位基因集合"改为"所有短于本记录 REF 的 ALT 都映射到 `*`"。

### 其余小类
- `refonly/refonly (5,)` 3 行：只差 QUAL（`163.67` vs `Infinity`）——与上文 7 个丢行同根
  （native 对 monomorphic 行的 QUAL 有限/非有限判定）；
- `refonly/refonly (7,)` 5 行：只差 INFO（`ExcessHet` 精度）。

**测量教训（本轮又踩一次）**：`/tmp` 在不同 bash 调用之间**不保留**，跨调用复用的产物必须写到
工作区（本轮改用 `.diag/round104/`）。

## 第 105 轮：跨接位点的 `*` 存活规则已被测出，7 个丢失行补齐（bug 45）

### 先证伪一个假说
第 104 轮末我曾猜"参考碱基 = N 的坐标不给 `*`"。实测证伪：跨接记录覆盖的 **165 个坐标里
162 个参考碱基就是 `N`，其中 145 个 GATK 照样给 `*`**，只有 17 个是 monomorphic。

### 实测规则（19 条跨接记录，17 条完全吻合，另 2 条由"该坐标另有记录起点"解释）
对一条非块记录，设来源样本基因型的等位基因长度为 `L_i`（含 REF 本身），记录 REF 长 `R`：

```
star 覆盖区间 = 该跨度最前面的 max_i(R - L_i) 个碱基
其余被覆盖坐标 = monomorphic 行（ALT=.，QUAL=Infinity，FILTER=.，INFO=DP;MLEAC=.;MLEAF=.）
```

| 来源 POS | R | 来源 GT | 基因型等位基因长度 | 预测 star | 实测 star | 实测 mono |
| --- | --- | --- | --- | --- | --- | --- |
| 10004769 | 11 | 0/1 | 11, 1 | 10 | 10 | 0 |
| 10006819 | 5 | 1/1 | 1, 1 | 4 | 4 | 0 |
| 10008952 | 14 | 2/4 | 3, 13 | **11** | **11** | 2 |
| 10011517 | 3 | 0/2 | 3, 2 | 1 | 1 | 1 |
| 10067049 | 10 | 2/3 | 2, 3 | 8 | 8 | 1 |
| 10076989 | 3 | 2/3 | 2, 4 | 1 | 1 | 1 |
| 10090289 | 3 | 0/2 | 3, 2 | 1 | 1 | 1 |
| 10097436 | 37 | 0/1 | 37, 1 | 36 | 34 | 1 |
| 10098308 | 3 | 0/0 | 3, 3 | 0 | 0 | 2 |

**只有"被调用"的等位基因决定存活范围**：10062935 的 ALTs 里有长度为 1/2/3 的删除等位基因，
但它的 GT 是 `0/4`（被调用的是长度为 5 的插入），实测 star = 0、3 个坐标全是 monomorphic ✓。

两条不完全吻合的记录都另有原因：10068158 与 10097436 的跨度内**另有记录起点**
（10068160、10097437），那些坐标由各自的记录支配（star+mono+被占 = 预测值 ✓）。

### 7 个丢失行的成因与修复
只有 3 条跨接记录的**被调用等位基因全都 ≥ REF 长度**（无 called deletion）：
10062935（GT 0/4）、10087820（GT 0/2）、10098308（GT 0/0）——它们覆盖 3+2+2 = **7 个坐标**，
正是 only-GATK 的那 7 行。native 的 `materialize_spanning_loci()` 原来要求
`called_deletion_allele_index >= 1` 才建立覆盖，于是这 7 个坐标**从未被访问**。

修复：
1. 覆盖遍不再丢弃 `deletion < 0` 的记录（保留 `deletion = -1`）；
2. 合成时投影目标改为"有 called deletion 时用该等位基因，否则用**本记录的 REF/<NON_REF>**"
   ——后者让所有具体 ALT 都成为 NO_CALL，重算出的基因型不携带 `*`，于是计算阶段像 GATK 的
   `filterAllelesForFinalSet()` 一样把 `*` 剪掉，落到既有的
   `materialize_gatk_monomorphic_ref_call()`（QUAL=Infinity、MLEAC/MLEAF 缺失、
   `FORMAT=GT:DP:RGQ`、样本 `0/0:<来源DP>:<来源GQ>`）。

**走过的死路（记录以免重犯）**：直接合成"1 等位基因"的行会撞上
`OUTPUT_CONTRACT_FAILURE: invalid reference-only record`——dense 的分组路径
（`genotype_gvcf_tool.cpp:7388-7391`）对组内每条记录无条件调用
`materialize_reference_only()`，而它要求 `record.reference_block == true`。**必须让记录带
`*` 走进计算阶段再被剪枝**，才能复用第 104 轮已验证的 monomorphic 出口。

### 验证（整段 chr20 dense）
`GATK 100000 位点 / native 100000 位点`（原 99993）、`only-GATK = []`、`only-native = []`、
`differ = 78`（无回归）；7 行**逐字节一致**。门禁
`fastgatk-genotype-gvcf-spanning-source-gq-gatk-oracle` 扩展为 14 个 pin 位点（7 个 `*` 行 +
7 个 monomorphic 行），并对 monomorphic 行额外断言 GATK 自身的 `ALT=.` 与 `QUAL=Infinity`。

### 第 105 轮更正：只有"带 GQ"的分支被落地（第一次全量 309/310 的教训）

本轮先把上面那条修复**无条件**落地，结果第一次全量回归
（`.diag/regression/20260913-101336/`）在两后端都 **309/310 失败**，失败者是既有的
`fastgatk-genotype-gvcf-dense-spanning-loci-gatk-oracle`。测出的原因是 GATK 对
"没有 `*` 贡献的被覆盖坐标"有**两种**形状，判别式是**来源 FORMAT 是否带 GQ**：

| 来源 FORMAT | GATK 实测行 | native（无条件落地时） |
| --- | --- | --- |
| 带 GQ（chr20：16/99/6） | `A . Infinity . DP=…;MLEAC=.;MLEAF=. GT:DP:RGQ 0/0:23:16` | ✓ 逐字节一致 |
| 不带 GQ（dense 夹具 `no-deletion-alt-no-star-rows`） | `A . . . DP=20 GT:AD ./.:0` | ✗ `A . Infinity . DP=20;MLEAC=.;MLEAF=. GT ./.` |
| 不带 GQ（夹具 `pruned-deletion-no-star-rows`） | `A . 192.21 . DP=20;MLEAC=.;MLEAF=. GT ./.` | ✗ `QUAL=Infinity` |

即：`cleanupGenotypeAnnotations()` 只有在被投影的基因型**仍带 GQ** 时才重新装回 hom-ref 调用
（`0/0:DP:RGQ` + QUAL=Infinity）；不带 GQ 时走另一支，发布 `./.` 且**去掉 QUAL**（no-deletion 夹具）
或保留有限 QUAL（剪枝夹具 192.21）。

因此最终落地的是**带 GQ 的分支**：覆盖遍加 `if (deletion < 0 && !record.source_has_gq) continue;`，
其余坐标**保持不发布**（缺口照旧，而不是发出错形状的行）。补上该条件后：
chr20 整段仍是 `100000/100000`、`only-GATK=[]`、`only-native=[]`、`differ=78`、10 个 pin 行逐字节一致；
两个门禁都零违规。

**下一轮的两个明确靶子**（都已定位到分支）：
1. **不带 GQ 的 no-call 分支**（`./.` + `FORMAT=GT:AD` + 无 QUAL / 有限 QUAL）——补齐它才能把
   "没有 `*` 贡献的被覆盖坐标"整类收干净，也能让 dense 门禁的两个 `no-star` 对照组从
   "只断言不含 `*`"升级为逐字节比对；
2. **`star/star` 38 行**的等位基因投影（把"短于本记录 REF 的每个**被调用**等位基因"都映射到 `*`，
   并按 GATK 的 `AlleleSubsettingUtils` 规则把 PL 折到新等位基因集合，实测目标值
   `1184,62,0` vs native 现输出 `1717,539,803`）。

## 第 106 轮：REF-only 行的 ExcessHet 渲染路径（bug 46），残差 78 → 73

### 现象（第 105 轮末的字段级地图又把这一类缩小到 7 行）
`(7,)` 类 7 行（INFO 单列差异）实测分成两种：

| 位点 | GATK | native | 差异键 |
| --- | --- | --- | --- |
| 10041698 / 10077008 / 10077010 / 10098308 / 10099270 | `ExcessHet=0.00` | `ExcessHet=0.0000` | `ExcessHet` |
| 10012573 / 10068159 | `QD=0.00` | `QD=-0.00` | `QD` |

### 机制（**源码级 + 探针级**双重确认，不是拟合）
1. **GATK 源码**：`ExcessHet.java:91` 与 `:272` 都把值**预先格式化成字符串**：
   `Collections.singletonMap(getKeyNames().get(0), (Object) String.format("%.4f", eh))`
   ⇒ htsjdk 原样写出 ⇒ 全语料 361 行都是 `0.0000` ✓（native 现状正确）。
2. **htsjdk 探针**（`javap` 找到 `VCFEncoder.formatVCFDouble(double)`，反射调用，pinned JDK17 + pinned jar）：

   | 输入 | 输出 |
   | --- | --- |
   | `0.0` | `0.00` |
   | `-0.0` | `0.00`（**符号丢失**） |
   | `1e-30` | `0.00` |
   | `1e-20` | `1.000e-20` |
   | `1e-16` | `1.000e-16` |
   | `0.005` | `5.000e-03` |
   | `0.5` | `0.500` |
   | `1.5` | `1.50` |
   | `12.26` | `12.26` |
   | `-3.891` | `-3.891e+00` |

   `formatQualValue`：`0.0→0`、`-0.0→-0`、`0.005→0.01`、`1.0→1`、`12.26→12.26`。
3. 只有**这 5 行 REF-only** 的 ExcessHet 是"裸数值属性"（由 htsjdk 的 `formatVCFDouble` 渲染 ⇒ 精确零 = `0.00`），
   其余 361 行走 GATK 自己 `%.4f` 的字符串路径 ⇒ `0.0000`。

### 修复（最小、可验证部分）
`format_gatk_float_value()` 增 `bool ref_only`；REF-only 且 `|值| < 1e-20` 时直接返回 `"0.00"`，
其余保持 4 位小数。参数由两处编码 lambda 在 `std::move(computed.record)` **之前**捕获
（`computed.record.allele_count == 1`）并透传。

**未实现的部分（同一条 htsjdk 路径的另一半）**：REF-only 行若带**非零** ExcessHet，GATK 会输出
`0.500`/`1.50`/`5.000e-03` 这类形式（探针已给出规则），而 native 仍输出 4 位小数。
chr20 语料里没有这种行，故**不落地未验证行为**，记录在此。

### 验证
整段 chr20 dense：位点级 **78 → 73**，5 行**逐字节一致**，`only-GATK=[]`、`only-native=[]`；
`(7,)` 类从 7 行降到 2 行（剩下的正是 QD 符号零）。门禁
`fastgatk-genotype-gvcf-spanning-source-gq-gatk-oracle` 增加第三组 pin（5 行 REF-only ExcessHet，
断言 GATK 侧确为 `ExcessHet=0.00` 且 native 行逐字节相同）。

### 下一个靶子：`QD` 的零符号（已定位到判别式）
16 行的 `QD` 是 `0.00`（GATK）vs `-0.00`（native），**全部是 `*` 行**：GATK 侧
`QD=-0.00` 的 115 行几乎都是 `GT 0/1`（`MLEAC=1`），`QD=0.00` 的 24 行几乎都是 `GT 1|1`/`1/1`
（`MLEAC=2`/`3`）⇒ 判别式是**位点是否 hom-alt**，而 htsjdk 的 `formatVCFDouble(-0.0)` 会丢符号
⇒ `-0.00` 必然来自另一个保留符号的路径/数值。native 15 行/24 行两侧都是 `-0.00`，
说明 native 的 **QD 数值本身**（或它派生的未取整 QUAL 的符号）在 hom-alt `*` 行上是错的。
下一轮应从"`*`-only 行的未取整 QUAL 符号"入手（AF 计算器 ~1e-16 的舍入方向），而不是格式化层。

## 第 107 轮：`star/refonly` 30 行的根因已定位（块记录的"内部坐标所有权"），修复试过但未达标，已回滚

### 根因（插桩实测）
用既有的 `FASTGATK_DEBUG_SUBSET`（`apply_gatk_output_allele_subset`，全部 100 kb 只有 387 行）
加临时 `[COV]`/`[SPAN]` 打印，对 10097436 那条记录（`REF=CTTTTCTTT…` 37bp、`GT 0/1`、
`ALT=C,<NON_REF>`、`spanEnd1=10097472`）实测：

```
[COV]  srcPos1=10097436 refLen=37 spanEnd1=10097472 deletion=1 gq_ok=1 pushed=1
[SPAN] srcPos1=10097436 … position1=10097438   ← 合成
[SPAN] srcPos1=10097436 … position1=10097439   ← 合成
[SPAN] srcPos1=10097436 … position1=10097440   ← 合成
[SPAN] srcPos1=10097436 … position1=10097441   ← 合成
（10097442 之后一行都没有）
```

跨度明明到 10097472，却只在 4 个坐标合成——原因是 `materialize_spanning_loci()` 的
`occupied` 集合由**所有记录的起点**构成，而 dense 模式在**解码阶段**就把每个参考块展开成
"每坐标一条 REF-only 记录"（`genotype_gvcf_tool.cpp:7225-7280`）⇒ **块的每个内部坐标都成了
某条记录自己的起点**，把跨接变异记录的 `*` 行全部挡掉。

这条规则同时解释了第 105 轮那个"找不到解释的 mid-span monomorphic 行"：
**10097442 正是那个块自己的起点**——GATK 在"有记录始于该坐标"时只保留起始记录
（`GenotypeGVCFsEngine.java:339-354`），于是该坐标是 REF-only 行；其余内部坐标则把块当作
**跨接事件**参与合并（REF→NO_CALL），所以跨接变异的 `*` 照样出现 ✓✓。

### 修复实验（已验证形状变对，但位点级不达标，**已回滚**）
给 Record 加 `expanded_block_interior`（块展开时 `position != start` 置真），并在 `occupied`
里跳过这些内部块记录（只保留真正的起点）。实测整段 chr20 dense：

| 类别 | 修复前 | 修复后 |
| --- | --- | --- |
| `star/refonly (4,5,6,7,8,9)` | 30 | **12** |
| `star/star (7,9)` | 38 | **56** |
| 位点级 differ | 73 | **73（未降）** |

即 **18 行的 ALT 形状恢复成 GATK 的 `*`** ✓，但它们的 INFO/样本列仍与 GATK 不同：

```
POS=10097450  GATK  ALT=* QUAL=0 LowQual  INFO=AC=1;AF=0.500;AN=2;DP=153;…  sample ./.:63,0:63:0:.:.:0,0,557
              NATIV ALT=* QUAL=0 LowQual  INFO=…;DP=90;…                      sample 0|1:57,32:89:99:0|1:10097436_…:1089,0,2243:10097436
```

### 下一个靶子：块 + 跨接变异在同一坐标的**合并语义**
GATK 在这些坐标上把**块也并进来**（`ReferenceConfidenceVariantContextMerger`：跨接变异的等位基因
变 `*`、块的 REF 变 NO_CALL），于是样本变成 `./.`、PL 是合并后的 `0,0,557`、位点 `DP` 也是合并值
153；native 现在只发出跨接变异记录自己的样本（`0|1` + PGT/PID/PS）。要收掉这 12+56 行，需要实现
"同一坐标多条记录（含块）的样本级合并"，而不是只取某一条记录——这与 `star/star (7,9)` 的
PL 投影问题是同一处的两半。

**回滚说明**：该实验不构成回归（`only-GATK=[]`、`only-native=[]`、differ 不变），但也没有把
位点级残差从 73 降下来，按纪律**不留未达标的半成品**；源码已 `git checkout` 还原，两后端已重建，
复测确认仍是 73 / `only-GATK=[]` / `only-native=[]`。

## 第 108 轮：剩下 73 处的**逐类机制**已按 GATK 源码测清（下一步可直接照做）

### 类 1：`star/star (7,9)` 38 行 —— 跨接位点的等位基因投影

对同一来源记录的多个 `*` 行逐位点对拍后，GATK 的映射规则是**可复现的**：

```
输出等位基因 = {位点参考碱基, "*"}
来源基因型第 i 个等位基因 → 若 i == 0 则映射到 REF；若 len(allele_i) < len(来源REF)（或该等位基因本身就是 "*"）
                           则映射到 "*"；否则 → NO_CALL（只把该等位基因置 missing，不整样本置 no-call）
AD[*] = AD[最后一个映射到 "*" 的被调用等位基因]
PL/GT/QUAL = 由 AF 模型**重新基因分型**（不是沿用来源 PL）
```

实测证据（`AD` 全是"最后一个映射等位基因"）：

| 位点 | 来源 GT / AD | GATK 行 | native 行 | native 取错在哪 |
| --- | --- | --- | --- | --- |
| 10008953 | 2/4 / 1,3,14,5,11,0 | `1/1` AD **1,11** PL 1184,62,0 | `./.` AD 1,14 PL 1717,539,803 | 只取**第一个**被调用删除（等位基因 2），等位基因 4 落到目标外 ⇒ `force_no_call_on_dropped_gt` 把整样本置 no-call |
| 10067050 | 2/3 / 2,4,7,9,0 | `1/1` AD **2,9** PL 544,28,0 | `./.` AD 2,7 | 同上 |
| 10068171 | 1\|2 / 0,2,4,0 | `1\|1` AD **0,4** | `./.` AD 0,2 | 同上 |
| 10076990 | 2/3 / 7,4,13,7,0 | `0/1` AD **7,13** PL 211,0,335 | `./.` AD 7,13 PL 252,41,376 | AD 对；但等位基因 3 落到目标外 ⇒ 强制 no-call；GATK 重新基因分型后给 `0/1` |

**缺口清单（实现时缺一不可）**：① 目标等位基因集合改为"所有短于来源 REF 的被调用等位基因 → 单个 `*`"；
② `force_no_call_on_dropped_gt` 在这种"部分等位基因落到 NO_CALL"的情形必须**关掉**；
③ AD 取**最后一个**映射等位基因（现在是第一个）；④ PL 必须由 AF 模型重算（native 现在沿用了来源 PL：
`1717,539,803` = 来源 PL 的 (0,0),(0,2),(2,2) 三个对角项）。

### 类 2：`star/refonly` 30 行 —— 块记录内部坐标的所有权 + **块与跨接变异的样本级合并**
第 107 轮已定位：块在解码阶段被展开成"每坐标一条记录"，其内部坐标被误当作记录起点。
把这一点修好后（实验已回滚）18 行的 ALT 变对，但样本列仍差：GATK **把块也并进同一坐标**
（跨接变异贡献 `*`、块的 REF 变 NO_CALL），于是样本是 `./.`、PL 是合并后的 `0,0,557`、位点 `DP` 是合并值 153；
native 只发跨接变异记录自己的样本（`0|1` + PGT/PID/PS）。
⇒ 需要实现"同一坐标多条记录（含块）的样本级合并"，与类 1 的投影是同一处的两半。

### 类 3：`refonly/refonly (5,)` 3 行 —— QUAL 有限 vs `Infinity`
GATK 源码 `GenotypingEngine.java:158-163` 实测：

```java
final double log10Confidence =
    !outputAlternativeAlleles.siteIsMonomorphic || annotateAllSitesWithPLs
        ? AFresult.log10ProbOnlyRefAlleleExists() + 0.0     // 非单态：补后验
        : AFresult.log10ProbVariantPresent() + 0.0;         // 单态：AF 模型自己的 P(变异存在)
final double phredScaledConfidence = (-10.0 * log10Confidence) + 0.0;   // +0.0 去掉 -0.0
```

而 `AFCalculationResult.java:118-119`：
`log10ProbVariantPresent() == MathUtils.log10OneMinusPow10(log10PosteriorOfNoVariant)`
——**与 native 的公式完全相同** ⇒ 差异不在公式，而在 **`log10PosteriorOfNoVariant` 的数值本身**：
native 得 `0.0`（精确 1）⇒ 补后验 `-inf` ⇒ QUAL=Infinity；GATK 得 `≈ -1.7e-17`
⇒ QUAL = 163.67。这是 AF 计算器 log10-sum-exp 的 ~1e-16 级浮点差（与第 51 轮记下的符号零同族），
**不可能靠改公式修掉**，需要对齐 AF 状态的求和实现。

### 类 4：`star/star (7,)` 2 行 —— QD 的零符号
GATK 侧统计：`QD=-0.00` 115 行（几乎都是 `GT 0/1`、`MLEAC=1`），`QD=0.00` 24 行（几乎都是
`GT 1|1`/`1\|1`、`MLEAC=2/3`）⇒ 判别式是**位点是否 hom-alt**。htsjdk 探针显示
`formatVCFDouble(-0.0) == "0.00"`（丢符号）而 GATK 大量输出 `-0.00`，说明 QD 走的是**保留符号**的
浮点路径（`%.2f`）⇒ 差异在**数值**：native 在两类行上都得到 ≈ -0.0。
注意 `GenotypingEngine` 里 `phredScaledConfidence` 加了 `+0.0` 去负零，而写进 VC 的是
`builder.log10PError(log10Confidence)`（没有那次 `+0.0`）——**QD 看的是后者**，这就是"QUAL 打印 `0`
但 QD 打印 `-0.00`"的来源。要修需要在 native 侧把"写进记录的 QUAL"与"QD 用的未取整 QUAL"分开，
并按 hom-alt/het 复现符号。

**本轮结论**：四类都测清了机制与缺口，但每一类都需要成体系的实现（投影重写 / 样本级合并 /
AF 求和实现 / QUAL 双值拆分），不适合在一轮里半成品落地。本轮**未改生产代码**，树保持
第 106 轮基线（位点级 73、双后端 310/310）。

## 第 109 轮：换一个**独立语料**跑，立刻暴露两个新的文本层缺陷（bug 47 / 48）

已调优的那份 chr20 语料（`expected.testGVCFMode.gatk4.g.vcf`，NA12878）已经把门禁跑绿，
于是换一份**native 从未针对过的真实 gVCF**：`CombineGVCFs/YRIoffspring.chr20snippet.g.vcf`
（GIAB 风格、同一条 chr20、另一个样本；它的 INFO 集合里有 **`ClippingRankSum`**，而调优语料没有）。
同一份 `-R testdata/chr20/reference/GRCh37.chr20.fa`，两种模式各跑 GATK / native 逐字节比对：

| 模式 | 修复前 differ | 修复后 differ |
| --- | --- | --- |
| **default（无 dense 标志）** | **132 / 300** | **0 / 300 逐字节一致** ✓ |
| dense（`--include-non-variant-sites`） | **285 / 100001** | **81 / 100001** |

### bug 47：`ClippingRankSum` 没有被当作"秩和族"注释
两份 GATK 行的差别只有两处，都指向同一个遗漏：

```
GATK  …;BaseQRankSum=-2.640e+00;ClippingRankSum=0.00;DP=35;…;SOR=0.446
NATIV …;BaseQRankSum=-2.640e+00;DP=35;…;SOR=0.446;ClippingRankSum=0
```

1. **键序**：native 把未知键追加到 INFO 末尾，而 GATK/hstsjdk 把它放在 `BaseQRankSum` 之后
   （`standard_order` 与头部 `info_order` 两张表都缺这一项）；
2. **数值格式**：`rank_sum` 判定只列了 BaseQRankSum/MQRankSum/ReadPosRankSum ⇒ `ClippingRankSum`
   落到"原样输出"分支，于是 `0` 而不是 `0.00`（也不做科学计数）；
3. **物化行抑制**：dense 下 materialized spanning 行必须丢掉读级秩和注释，
   `suppress_materialized_spanning_annotations()` 的名单同样缺 `ClippingRankSum`
   ⇒ native 在 GATK 不写该键的行上多写了一个 `ClippingRankSum=0.00`（修复后 dense 138 → 81）。

三处都是"这个键没被认成秩和族"的同一个根，本轮一起修掉。

### bug 48：REF-only 行的 `ExcessHet` 非零渲染（第 106 轮规则的补全）
第 106 轮只落地了"REF-only 且值为 0 ⇒ `0.00`"这一半（当时没有非零样本）。这份语料给出了非零证据：

```
POS=10002166  GATK …;ExcessHet=3.01;…        NATIVE …;ExcessHet=3.0103;…
```

`3.0103 → "3.01"` 正是 htsjdk `formatVCFDouble` 对 ≥1 值的两位小数规则（探针已验），
于是把 REF-only 分支补成完整的 `formatVCFDouble` 规则：`|v|<1e-20 → "0.00"`、
`|v|<0.01 或负数 → %.3e`、`0.01≤|v|<1 → 三位`、`|v|≥1 → 两位`
（0.0→0.00 ✓ chr20、3.0103→3.01 ✓ 本语料、0.5→0.500、0.005→5.000e-03、-3.891→-3.891e+00 均由探针验证）。

### 独立语料同时**证实**了第 108 轮的 QD 符号规则
dense 剩下的 81 行：`star/star (7,)` 51 行 + `star/star (7,9)` 30 行，
**列只用到 7(INFO) 与 9(样本)**，没有任何形状列差异；INFO 键差异里 `QD` 出现在 75 行。
其中 QD 差异全部是**同一位点是否 hom-alt**：

```
POS=10006820  GATK QD=0.00    NATIVE QD=-0.00     （行：AC=2;AF=1.00;AN=2 → hom-alt）
```

⇒ 与第 108 轮在 chr20 上得到的判别式一致（hom-alt → `0.00`，het → `-0.00`），
两个语料合计 ~97 行。结合 `GenotypingEngine.java:158-163` 的 `+ 0.0` 归一与
`builder.log10PError(log10Confidence)` 写值路径，可以判定：**QD 符号与 `QUAL=163.67 vs Infinity`
（类 3）同根**——都是 native 的 AF 后验在 ~1e-17 量级上"饱和成精确 0/1"，而 GATK 保留了那个极小的
符号。⇒ 下一轮只需攻**一处**（AF 状态 log10-sum-exp 的求和实现），即可同时收掉类 3 与类 4。

### 新门禁
`fastgatk-genotype-gvcf-independent-corpus-gatk-oracle`
（`scripts/verify_genotype_gvcf_independent_corpus_gatk_oracle.py`）：
default 模式断言**整行逐字节一致**（含行数对照 300）；dense 模式断言
①行集合一致（不缺失/不多出）、②**形状列 REF..FILTER 必须一致**、③INFO 差异键必须落在已知集合内
——即把当前残差"钉形状"，任何新形状类或新注释键都会立刻报警。
