# 第 53 轮：dense 记录供给裁剪到 `-L` 区间 —— 真实语料与最小夹具证据

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
