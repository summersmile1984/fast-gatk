# 第 52 轮：dense 跨位点物化在**真实语料**上的验证（chr20）

本轮修复（`materialize_spanning_loci()`）用 GATK 自带的 HaplotypeCaller gVCF 语料验证，
而不是只用合成夹具。

## 输入

| 项 | 值 |
| --- | --- |
| 输入 gVCF | `gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/haplotypecaller/expected.testGVCFMode.gatk4.g.vcf`（1291 条记录，contig `20`） |
| 参考 | `testdata/chr20/reference/GRCh37.chr20.fa`（`.fai` 名为 `20`，与 FASTA 头 `>20` 一致） |
| oracle | pinned GATK 4.6.2.0（`third_party/jdk17` + `third_party/gatk-package/gatk-4.6.2.0/...-local.jar`） |
| native | `fastgatk-native/build/fastgatk-genotype-gvcf --gatk-compatible-annotations` |
| 模式 | `--include-non-variant-sites`（dense），两侧同一 `-L` 窗口 |

## 窗口 20:10000000-10003000（含跨接删除记录）

| 指标 | 结果 |
| --- | --- |
| GATK 行数 / 位点数 | 3001 / 3001 |
| native 行数 / 位点数 | 3021 / 3021 |
| **GATK 有而 native 没有的位点** | **0** ✅（修复前这一整类坐标会缺失） |
| **`*` 行逐字节不一致的位点** | **0** ✅（本轮新增的行与 GATK 完全一致） |
| 其余不一致位点 | 69 个 REF-only 行（见下，**既有**分歧）、9 个 other |
| native 多出的位点 | 20 个，全部是窗口末端之后的 REF-only 行（既有分歧，见下） |

分类统计（`(kind, has-N)` → 条数）：`('refonly', False): 73`、`('other', False): 9`；
其中 `star` 类为 **0**。

## 同一语料暴露的两处**既有**分歧（与本轮改动无关，未修）

### (1) dense 的参考块展开不遵守 `-L` 区间

窗口 `20:10050000-10051000`（1000 bp，且该窗口内没有删除等位基因，故本轮的 pass 完全不触发）：

| 工具 | 行数 | 界内行 | 越界行（POS > 10051000） |
| --- | --- | --- | --- |
| GATK | 1001 | 1001 | 0 |
| native | 4809 | 4362 | 447 |

GATK 恰好每个坐标一行；native 在 1000 bp 窗口内产出 4362 行（**同一位点多行**），
并在区间末端之外又产出 447 行。窗口 1 里 native 多出的 20 行（`20:10003001-10003020`，
全部 `ALT='.'`、`REF=N`）是同一机制。

### (2) 参考块展开行的 REF 取错来源

窗口 1 中 69 个位点上两侧都是 REF-only 行（`ALT='.'`），但 REF 字母不同：
**GATK 写该位点自己记录的 REF（`T`/`C`/`A`），native 写 FASTA 的碱基（`N`）**。
实测：`20:10000000` 在 `GRCh37.chr20.fa` 中确实是 `N`，而该位点的 gVCF 记录是
`20 10000000 . T <NON_REF>`；GATK 输出 `T`，native 输出 `N`。
即：位点若有自己的记录，GATK 用记录的 REF；native 的块展开路径回落到 FASTA 碱基。

## 复现

```bash
# 索引语料（一次）
third_party/jdk17/bin/java -Xmx2g -jar third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
  IndexFeatureFile -I <语料 gVCF>

# 两侧同窗口
third_party/jdk17/bin/java -Xmx2g -jar <jar> GenotypeGVCFs -R testdata/chr20/reference/GRCh37.chr20.fa \
  -V <语料 gVCF> --include-non-variant-sites -L 20:10000000-10003000 -O gatk.vcf --create-output-variant-index false
fastgatk-native/build/fastgatk-genotype-gvcf -R testdata/chr20/reference/GRCh37.chr20.fa \
  -V <语料 gVCF> --include-non-variant-sites --gatk-compatible-annotations -L 20:10000000-10003000 -O native.vcf
```

比较方式：按 POS 建索引，逐位点比较「行集合」，并把差异按 `*` 行 / REF-only 行分类。
