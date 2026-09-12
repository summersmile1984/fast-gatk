# 第 53–54 轮：dense 记录供给的 `-L` 裁剪 与 块起始坐标的 REF 来源 —— 真实语料与最小夹具证据

## 最小夹具（合成）

输入：一条纯参考块记录 `chr1 2 . A <NON_REF> . . END=11;DP=40 GT:DP:GQ:MIN_DP 0/0:40:99:40`
（100 bp 全 `A` 的 chr1 参考），两侧同一 `-L chr1:5-7`、`--include-non-variant-sites`。

| | 行数 | 位点 |
| --- | --- | --- |
| GATK 4.6.2.0 | 3 | 5, 6, 7 |
| native（修前） | 10 | 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 |
| native（修后） | 3 | 5, 6, 7 |

行文本（两侧逐字节一致）：`chr1\t<pos>\t.\tA\t.\t.\t.\tDP=40\tGT:DP:RGQ\t0/0:40:99`

该夹具已成为注册门禁 `fastgatk-genotype-gvcf-dense-spanning-loci-gatk-oracle` 的用例
`dense-reference-block-clipped-to-intervals`（先确认其在修前失败：`record count differs: GATK=3 native=10`）。

## 真实语料（GATK 自带 chr20 HaplotypeCaller gVCF，1291 条）

参考 `testdata/chr20/reference/GRCh37.chr20.fa`；dense 模式；两侧同一窗口。

| 窗口 | GATK 行 / 位点 | native 修前 | native 修后 |
| --- | --- | --- | --- |
| `20:10000000-10003000` | 3001 / 3001 | 3021 / 3021（位点 10000000-10003020） | **3001 / 3001**（位点 10000000-10003000） |
| `20:10050000-10051000` | 1001 / 1001 | 4809 / 4809（位点 10046639-10051447） | **1001 / 1001**（位点 10050000-10051000） |

修后逐位点比较：

| 窗口 | 仅 GATK 有的位点 | 仅 native 有的位点 | 共同但行不同的位点 |
| --- | --- | --- | --- |
| `20:10000000-10003000` | 0 | 0 | 82 |
| `20:10050000-10051000` | 0 | 0 | 3 |

剩余 85 个差异**全部是 REF-only 行且只有 REF 字母不同**（下节），本轮的裁剪修复已把
「行数/位点集合」这一层做平。

## 剩余差异的精确定位（下一轮目标，仍未修）

**块展开行的 REF 取错来源。** 例：

```
POS 10000000   GATK   : 20	10000000	.	T	.	.	.	DP=64	GT:DP:RGQ	0/0:64:99
               NATIVE : 20	10000000	.	N	.	.	.	DP=64	GT:DP:RGQ	0/0:64:99
```

该位点的 gVCF 记录是 `20 10000000 . T <NON_REF>`，而 `GRCh37.chr20.fa` 在该坐标为 `N`
（语料与该参考不一致），**GATK 用记录自身的 REF（`T`），native 用 FASTA 碱基（`N`）**。
窗口 2 的 3 个同类位点同理（块内出现 `A` vs `N`）。

判据：**块自身的起始坐标应保留记录的 REF**（那个位点就是记录本身），
只有块内部坐标才用 FASTA 碱基。

## 复现命令

```bash
JAR=third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar
R=testdata/chr20/reference/GRCh37.chr20.fa
V=gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/haplotypecaller/expected.testGVCFMode.gatk4.g.vcf
third_party/jdk17/bin/java -Xmx2g -jar $JAR IndexFeatureFile -I $V          # 一次

third_party/jdk17/bin/java -Xmx2g -jar $JAR GenotypeGVCFs -R $R -V $V \
  --include-non-variant-sites -L 20:10050000-10051000 -O gatk.vcf --create-output-variant-index false

fastgatk-native/build/fastgatk-genotype-gvcf -R $R -V $V \
  --include-non-variant-sites --gatk-compatible-annotations -L 20:10050000-10051000 -O native.vcf
```

比较：按 POS 建索引后比较「行集合」，并统计仅 GATK / 仅 native / 共同但不同的位点数。

---

# 第 54 轮：块起始坐标必须保留记录自身的 REF（第 36 个已修 bug）

## 规则与最小夹具

locus 的合并 REF 取自**在该坐标起始的记录**；只有在没有记录起始的坐标才回落到参考碱基
（merger 收到 `ref.getBase()`，`ReferenceConfidenceVariantContextMerger.merge()`）。

夹具：全 `A` 的 chr1 参考 + 一条块记录 `chr1 2 . C <NON_REF> . . END=5;DP=40 GT:DP:GQ:MIN_DP 0/0:40:99:40`
（记录的 REF 与参考碱基故意不一致），dense 模式：

| | POS 2 | POS 3 | POS 4 | POS 5 |
| --- | --- | --- | --- | --- |
| GATK | **C** | A | A | A |
| native（修前） | A ✗ | A | A | A |
| native（修后） | **C** | A | A | A |

该夹具已成为注册门禁用例 `dense-block-start-keeps-the-record-ref`
（修前确认失败：`row 0 ... GATK='chr1\t2\t.\tC...' NATIVE='chr1\t2\t.\tA...'`）。

修复落在两处：解码阶段的逐坐标展开（仅当 `position != start` 时才用 FASTA 碱基），
以及 `split_reference_blocks_at_variants()` 中起点等于块起点的分段（保留 `original.alleles.front()`）。

## 真实语料复验

| 窗口 | 修前差异位点 | 修后差异位点 | 其中 |
| --- | --- | --- | --- |
| `20:10000000-10003000` | 82 | **13** | 9 个只有 QD 不同；4 个 native 发 REF-only 而 GATK 发变异行（**既有**） |
| `20:10050000-10051000` | 3 | **0** | 两个窗口的「仅 GATK 有 / 仅 native 有」位点仍为 0 |

窗口 2 现在**逐字节完全一致**（1001/1001 行）。

## 剩余两类差异（已定位，未修）

1. **QD 抖动（9 个位点）**：两侧等位基因相同、只有 `QD` 不同。
   例：`20:10001298` GATK `QD=27.24`、native `QD=30.97`；而 `20:10000758` GATK 恰为 `30.97`。
   即 `QualByDepth.fixTooHighQD()`（`QD >= 35` 时映射为 `30 + gaussian*3`）在两边的
   **随机数消耗顺序/次数不同**，native 目前按自己的顺序抽取。
2. **变异行退化为 REF-only 行（4 个位点，既有）**：例 `20:10000758`
   输入是 `T A,<NON_REF>`（同型变异，MLEAC=2），GATK 输出 `T A ... AC=2;AF=1.00;...`，
   native 输出 `T . 3853.06 . DP=97;...;MLEAC=.;MLEAF=.`——变异 ALT 被剪掉、
   走了 REF-only 物化路径。**与第 54 轮的改动无关**：该行的修前/修后输出逐字节相同。
