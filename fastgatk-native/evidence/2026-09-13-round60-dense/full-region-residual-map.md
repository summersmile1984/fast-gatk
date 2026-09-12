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
