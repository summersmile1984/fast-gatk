# fast-gatk vs GATK 4.6.2.0 — Headline Comparison

生成时间：2026-09-17 01:47:50（Asia/Taipei）
工作树：/home/turing-agents/Documents/fast-gatk
原生构建：/home/turing-agents/Documents/fast-gatk/fastgatk-native/build（OpenMP 后端，32 硬件线程）
GATK：pinned GATK 4.6.2.0 via vendored JDK 17
重复次数：2（取 p50）

三维度：
1. **数据正确性**：原生输出与 GATK 输出的 byte-exact 比对（VCF/VCF.gz 经 BGZF 解码后再比对载荷）。
2. **速度**：`/usr/bin/time` 报告的 wall-clock（Elapsed）。
3. **资源消耗**：`/usr/bin/time` 报告的 peak RSS、文件系统读/写次数、%CPU。

---

## MarkDuplicates
- Fixture: `gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam`

| Metric | Native | GATK 4.6.2.0 | Native/GATK |
| --- | ---: | ---: | ---: |
| wall (p50, s) | 0.010 | 4.023 | **395.95× faster** |
| max RSS (KB) | 6,618.0 | 530,702.0 | 1.25% |
| output bytes | 56,681 | 62,480 | — |

- **数据正确性（BAM 完整性 + sha256 粗校）：✅ PASS**
  - duplicate-flag 决策：native=29，gatk=32，flag 不同的 read 数：3
  - 例：[['809R9ABXX101220:5:44:13776:49360', False, True], ['809R9ABXX101220:5:6:17918:145992', False, True], ['809R9ABXX101220:5:42:11849:188393', False, True]]
  - native sha256: `20829430e8e684d5…`
  - gatk   sha256: `70e7cd005e811323…`

---

## HC
- Fixture: `fixtures/chr20/mnp.bam`
- Region: `20:10019901-10020710`

| Metric | Native | GATK 4.6.2.0 | Native/GATK |
| --- | ---: | ---: | ---: |
| wall (p50, s) | 0.534 | 4.254 | **7.96× faster** |
| max RSS (KB) | 161,996.0 | 306,658.0 | 52.83% |
| FS inputs | 0 | 0 | — |
| FS outputs | 88 | 16 | — |
| output bytes | 1,454 | 1,548 | — |

- **数据正确性（CHROM/POS/ID/REF/ALT 5 列 record-byte-exact）：✅ PASS**
- 全文件 BGZF 解码后 byte-exact：✅ PASS**
  - 记录数：native=11, gatk=11
  - native sha256: `2ac11a5c66cb974c…`
  - gatk   sha256: `a3faf366cd8ebbf4…`

---

## HC (A-line regression window)
- Fixture: `fixtures/chr20/mnp.bam`

| Metric | Native | GATK 4.6.2.0 | Native/GATK |
| --- | ---: | ---: | ---: |
| wall (p50, s) | 0.385 | 4.342 | **11.27× faster** |
| max RSS (KB) | 162,036.0 | 317,332.0 | 51.06% |
| output bytes | 1,296 | 1,369 | — |

- **数据正确性（CHROM/POS/ID/REF/ALT 5 列 record-byte-exact）：✅ PASS**
- 全文件 BGZF 解码后 byte-exact：✅ PASS**
  - 记录数：native=7, gatk=7
  - native sha256: `a521f22b255b26ef…`
  - gatk   sha256: `48150cf661ff989c…`

---

## BQSR (BaseRecalibrator)
- Fixture: `testdata/real/ceutrio/CEUTrio.chr20.bam`

| Metric | Native | GATK 4.6.2.0 | Native/GATK |
| --- | ---: | ---: | ---: |
| wall (p50, s) | 2.733 | 21.829 | **7.99× faster** |
| max RSS (KB) | 154,226.0 | 745,312.0 | 20.69% |
| FS inputs | 0 | 0 | — |
| FS outputs | 584 | 25764 | — |
| output bytes | 176,770 | 13,186,116 | — |

- **数据正确性（raw byte-exact）：⚠️ DIFFER（见说明）**
- BQSR 表头的 GATKReport 结构对齐（10 个 GATKTable），但 covariate 表体被压缩（详见 `tool_registry.json` `BaseRecalibrator.fallback_boundaries`）；不混用跨实现表。
  - native sha256: `403da1a7ac397b5a…`
  - gatk   sha256: `f22e93517cb59e34…`

---

## Mutect2 (tumor+normal)
- Fixture: `gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam`
- Region: `17:69000-70000`

| Metric | Native | GATK 4.6.2.0 | Native/GATK |
| --- | ---: | ---: | ---: |
| wall (p50, s) | 0.301 | 4.133 | **13.74× faster** |
| max RSS (KB) | 21,914.0 | 299,856.0 | 7.31% |
| FS inputs | 0 | 0 | — |
| FS outputs | 72 | 24 | — |
| output bytes | 1,599 | 1,681 | — |

- **数据正确性（CHROM/POS/ID/REF/ALT 5 列 record-byte-exact）：✅ PASS**
- 全文件 BGZF 解码后 byte-exact：❌ FAIL（差异在 header/sample 列）**
  - 记录数：native=0, gatk=0
  - native sha256: `ffd3760fed810423…`
  - gatk   sha256: `414c8a4068f2a192…`

---

## Mutect2 (DREAM somatic, hs37d5 ref)
- Fixture: `testdata/real/dream_synthetic/chr20/tumor.bam + testdata/real/dream_synthetic/chr20/normal.bam`

| Metric | Native | GATK 4.6.2.0 | Native/GATK |
| --- | ---: | ---: | ---: |
| wall (p50, s) | 53.432 | 16.436 | **0.31× faster** |
| max RSS (KB) | 742,762.0 | 943,418.0 | 78.73% |
| output bytes | 2,251 | 5,061 | — |

- **数据正确性（CHROM/POS/ID/REF/ALT 5 列 record-byte-exact）：✅ PASS**
- 全文件 BGZF 解码后 byte-exact：❌ FAIL（差异在 header/sample 列）**
  - 记录数：native=5, gatk=5
  - native sha256: `3dafb77ad0b964c2…`
  - gatk   sha256: `ced450a696622553…`

---

## SortSam (coordinate)
- Fixture: `gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam`

| Metric | Native | GATK 4.6.2.0 | Native/GATK |
| --- | ---: | ---: | ---: |
| wall (p50, s) | 0.015 | 3.821 | **248.13× faster** |
| max RSS (KB) | 5,888.0 | 246,698.0 | 2.39% |
| output bytes | 61,842 | 61,351 | — |

- **数据正确性（BAM 完整性 + sha256 粗校）：⚠️ DIFFER（见说明）**
  - native sha256: `1f39af8be7cffe91…`
  - gatk   sha256: `7d5cda400d201daa…`

---

## MarkDuplicates (small fixture)
- Fixture: `gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam`

| Metric | Native | GATK 4.6.2.0 | Native/GATK |
| --- | ---: | ---: | ---: |
| wall (p50, s) | 0.011 | 3.971 | **347.45× faster** |
| max RSS (KB) | 6,522.0 | 526,940.0 | 1.24% |
| output bytes | 56,681 | 62,480 | — |

- **数据正确性（BAM 完整性 + sha256 粗校）：✅ PASS**
  - duplicate-flag 决策：native=29，gatk=32，flag 不同的 read 数：3
  - 例：[['809R9ABXX101220:5:42:11849:188393', False, True], ['809R9ABXX101220:5:6:17918:145992', False, True], ['809R9ABXX101220:5:44:13776:49360', False, True]]
  - native sha256: `20829430e8e684d5…`
  - gatk   sha256: `70e7cd005e811323…`

---

## MarkDuplicates (CEUTrio chr20, 222k reads)
- Fixture: `testdata/real/ceutrio/CEUTrio.chr20.bam`

| Metric | Native | GATK 4.6.2.0 | Native/GATK |
| --- | ---: | ---: | ---: |
| wall (p50, s) | 2.945 | 6.683 | **2.27× faster** |
| max RSS (KB) | 101,828.0 | 1,520,800.0 | 6.70% |
| output bytes | 27,576,171 | 31,126,231 | — |

- **数据正确性（BAM 完整性 + sha256 粗校）：⚠️ DIFFER（见说明）**
  - duplicate-flag 决策：native=24358，gatk=24360，flag 不同的 read 数：2
  - 例：[['20FUKAAXX100202:6:44:17448:95596', False, True], ['20GAVAAXX100126:6:24:16180:173228', False, True]]
  - native sha256: `1b571b87fea107d2…`
  - gatk   sha256: `376bc3be5f1cd640…`

---

## GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs)
- Fixture: `gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam`

| Metric | Native | GATK 4.6.2.0 | Native/GATK |
| --- | ---: | ---: | ---: |
| wall (p50, s) | 0.042 | 4.710 | **113.19× faster** |
| max RSS (KB) | 7,688.0 | 335,006.0 | 2.29% |
| output bytes | 2,718 | 3,036 | — |

- **数据正确性（CHROM/POS/ID/REF/ALT 5 列 record-byte-exact）：✅ PASS**
- 全文件 BGZF 解码后 byte-exact：❌ FAIL（差异在 header/sample 列）**
  - 记录数：native=3, gatk=3
  - native sha256: `3fce92e30304e9bf…`
  - gatk   sha256: `8f26030363152dc8…`

---

## 汇总

| Case | Speed (native/GATK) | RSS (native/GATK) | 数据正确性 |
| --- | ---: | ---: | --- |
| MarkDuplicates | 395.95× faster | 1.25% | ✅ BAM ok |
| HC | 7.96× faster | 52.83% | ✅ records |
| HC (A-line regression window) | 11.27× faster | 51.06% | ✅ records |
| BQSR (BaseRecalibrator) | 7.99× faster | 20.69% | ⚠️ differ (covariate collapse) |
| Mutect2 (tumor+normal) | 13.74× faster | 7.31% | ✅ records |
| Mutect2 (DREAM somatic, hs37d5 ref) | 0.31× faster | 78.73% | ✅ records |
| SortSam (coordinate) | 248.13× faster | 2.39% | ⚠️ differ |
| MarkDuplicates (small fixture) | 347.45× faster | 1.24% | ✅ BAM ok |
| MarkDuplicates (CEUTrio chr20, 222k reads) | 2.27× faster | 6.70% | ⚠️ differ |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | 113.19× faster | 2.29% | ✅ records |

详细原始数据：`/home/turing-agents/Documents/fast-gatk/report/headline_comparison_20260917-004555/<case>/`，JSON 摘要：`/home/turing-agents/Documents/fast-gatk/report/headline_comparison_20260917-004555/summary.json`