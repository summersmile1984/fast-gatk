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

## 最近验证状态

| 验证项 | OpenMP | Serial | 备注 |
| --- | --- | --- | --- |
| 四倍体多 ALT / max ALT / max genotype oracle | 通过 | 本轮改动前通过 | 最新 Host 注入改动后仍需复跑 Serial。 |
| chr20 MNP + gVCF `*` + max ALT oracle | 通过 | 本轮改动前通过 | 新增的重叠 given-allele 用例已在 OpenMP 通过。 |
| 通用 HC `--alleles` oracle | 通过 | 待本轮复跑 | 含多 ALT、indel、filtered feature、gVCF。 |
| HC/Mutect2 完整回归 | 最近一次改动前 72/72 通过 | 未作为本轮证据 | 最新 Host 注入、gVCF 修改后必须重跑。 |

> 注意：表中“本轮改动前通过”不是最新实现的完成证据，只说明相同方向的历史基线。最终结论以最新二进制的完整回归为准。

## 未完成事项

### 近期必须完成

1. 完成最新改动后的 Serial 编译与下列 oracle：

   ```bash
   python3 fastgatk-native/scripts/verify_hc_alleles_gatk_oracle.py
   python3 fastgatk-native/scripts/verify_hc_chr20_max_mnp_gvcf_gatk_oracle.py \
     --native fastgatk-native/build-serial/fastgatk-hc-call
   python3 fastgatk-native/scripts/verify_hc_complex_multiallelic_oracle.py \
     --native fastgatk-native/build-serial/fastgatk-hc-call
   ```

2. 运行最新二进制的完整 HC/Mutect2 回归，不能使用“修改前 72/72”代替：

   ```bash
   third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest \
     --test-dir fastgatk-native/build --output-on-failure \
     -R 'fastgatk-(hc|mutect2)'
   ```

3. 对 gVCF 流式/分区路径扩大 oracle 覆盖，特别是多个重叠 indel、符号等位基因、注释 read-overlap margin 与重分块的组合。

### 仍需审计的核心一致性边界

- `--alleles` 的更多复杂注入情形：重复序列中的 indel、跨 AssemblyRegion 边界、多个相邻或重叠 feature event、base-haplotype 排名并列。
- 复杂 EventMap 与图路径的完整性：Native 当前以已物化的 K-best 路径和 Host CIGAR 作为受限注入母体；需要用更多 GATK fixture 证明其在高复杂区域等价。
- Mutect2：需要在最新 HC 共享 calling pipeline 改动后重跑完整 oracle，确认正常样本重放、联合 AssemblyResultSet 和 somatic gVCF 不受影响。
- BQSR 与其余核心路径的扩大真实数据/极端参数验证；当前已有实现和覆盖，但尚不足以声称全参数空间 1:1。

### 暂不优先的范围

- 长尾小工具和不影响 HC / Mutect2 / BQSR 加速库核心价值的工具面。
- 这些工具不是“已完成”；只是依据项目目标被明确排在核心一致性之后。

## 下一步执行顺序

1. 结束并确认 Serial 构建，完成 Serial oracle。
2. 跑完整 HC/Mutect2 回归；若有失败，先以 GATK 输出和源码定位到具体边界。
3. 将通过或失败的命令、fixture、GATK 版本和差异行追加到本文件。
4. 只有当每个核心选项都有匹配范围的当前证据时，才讨论“1:1 完成”。

