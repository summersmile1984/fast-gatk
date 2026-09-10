# fast-gatk 模块实施蓝图：C++ Host + Kokkos Kernel

> **ApplyBQSR report-owned ContextCovariate（2026-09-01）**：修正 consumer 端把
> mismatch context 固定为 2-mer 的模型漂移；native 现在从 GATKReport `Arguments`
> 解析并校验 `mismatches_context_size`（1..13），再以相同宽度构造 ContextCovariate，
> flat quality/delta 仍由 Kokkos RangePolicy 变换。新增 size-3 高覆盖 fixture，确保
> context posterior 实际改变 QUAL；pinned GATK 4.6.2.0 的五张 report table 与
> 200/200 Java-decoded ApplyBQSR records 在 OpenMP/Serial 均严格比较。其他自定义
> covariate/plugin 与云端模型语义仍保持 explicit fallback。

> **全量回归快照（2026-09-01）**：OpenMP 与 Serial 当前各 209 项 CTest；此前 194 项
> 已连续 194/194 通过（1198.16 s / 1184.40 s），新增 DepthOfCoverage deletion-site
> 与 MarkDuplicates pair-key oracles 在双后端专测通过。该结果只表示配置级回归稳定；完整 Java
> 算法等价、bit-identical、GPU 实机与多节点生产语义仍按各模块 fallback 说明执行。

> **AnnotateIntervals feature-query boundary（2026-09-01）**：native 补齐 GATK
> `--feature-query-lookahead`（含 0/负值），并在 OutputManifest/telemetry 记录查询缓存
> 控制；pinned Java/native interval oracle 双后端通过，annotation 数值保持不变。云端
> FeatureManager 查询仍为 explicit fallback。

> **HC/Mutect2 CIGAR indel activity ownership（2026-09-01）**：共享 Kokkos
> ActivityProfile 输入新增按 locus 对齐的显式 CIGAR I/D 支持计数；Host 从 ReadBatch
> 统一抽取 indel anchor，kernel 将其作为 active seed。这样仅有 `60M2D58M` 等 CIGAR
> indel、没有 mismatch 的 reads 也会进入 GATK 等价 AssemblyRegion/PairHMM 路径。pinned
> GATK 4.6.2.0 synthetic oracle 在 OpenMP/Serial 均得到 `chr1:250 GCA>G` 及相同
> genotype/likelihood payload，未承诺完整 graph/assembly parity。

> **BaseRecalibrator long-read CycleCovariate（2026-09-01）**：native 现在与 GATK
> 默认 `maximum-cycle-value=500` 一致拒绝超限 cycle，并支持 `--maximum-cycle-value` /
> `--max-cycle` 扩展确定性 cycle 域；checkpoint、report Arguments 和 manifest 保留该
> 维度。BQSR Host projection 同时将 CIGAR insertion read bases 作为 GATK `isSNP=0`
> 的 substitution observations 纳入质量/读组/cycle 表。1,500M + `50M2I48M` pinned
> GATK 4.6.2.0 oracle 在 OpenMP/Serial 的四张 report table 均 bit-identical；完整云端
> 分布式 BQSR 仍是 explicit fallback。

> **PairHMM no-tag indel quality（2026-09-01）**：修正 GATK 默认 `BI/BD` 缺失时的
> raw Q45 规则；Conservative PCR 以原始 read 坐标只调整非终位，clipping 后的末位若
> 原本是内部重复位点可为 Q20/Q40，`NONE` 模式保持 Q45。HC/Mutect2 共享 Host
> preparation，706-row pinned Java PairHMM 输入 oracle 双后端通过，显式 BI/BD 路径不变。

> **ModelSegments 多样本 kernel（2026-09-01）**：native 现在保留完整
> per-sample copy-ratio/MAF 向量，线性 kernel 使用拼接向量，Gaussian kernel 对各样本
> anchor kernel 求和；这与 GATK 4.6.2.0 的 MultisampleMultidimensionalKernelSegmenter
> 一致，避免反相关样本事件在均值归约中被错误抵消。双后端 synthetic oracle exact；
> 完整多样本 posterior/SVD raw-bit parity 仍为 explicit fallback。

> **VariantRecalibrator retry/iteration boundary（2026-09-01）**：pinned GATK
> 4.6.2.0 将 `--max-attempts`（模型构建重试，默认 1）与 `--max-iterations`
> （VBEM 迭代上限，默认 150）分开；native 已分别传递到模型构建/迭代路径并写入
> OutputManifest。OpenMP/Serial help/runtime oracle 均通过；完整 VBEM/raw-bit parity
> 仍保持 explicit fallback。

> **VariantRecalibrator VBEM raw-bit 对齐（2026-09-03）**：native 现已在
> `fastgatk-native/src/vqsr_fastmath.hpp`（commons-math3 3.5 FastMath.log 全表
> 反射转储 + Gamma.digamma 精确移植，探针 17/17 位级一致）与
> `fastgatk-native/src/vqsr_jama_lu.hpp`（GATK 捆绑 JAMA 1.0.3 left-looking LU
> decompose/det/solve 移植，随机 SPD 探针 300/300 位级一致）中固化 Java 的
> 数值路径；M-step 的 sumProb 1e-10 基数序、per-datum regCovar 结合序与
> final `timesEquals(1/sumProb)` 均已按 Java 关联顺序改写，`log_determinants`
> 语义改为 log(LU 积)，E-step/score kernel 改为 Java log10 域组合（cachedDenom
> 每项除以 log10、两阶段 crossProdTmp 二次型、`10^x` 归一、log10sumLog10、
> eval denom = log10(pow(2π,-N/2)) + log10(pow(det,-0.5))），VQSLOD 按 GATK
> `String.format("%.4f", lod)` 文本写出（不再 float 编码）。**pinned fixture 上
> VQSLOD 文本 raw-bit（delta 0.0，oracle 收紧为 token 相等）**，双后端 7/7 VR
> oracle 通过。已知边界：模型 report（%.16E）正模型迭代路径残余 ~5e-9
> （≤1e-6 gate），raw-bit 收益低未继续追。

> **VariantRecalibrator annotation ordering（2026-09-01）**：native 在归一化后按
> GATK `VariantDataManager` 的 training/non-training 均值偏移重排 annotation dimensions，
> 使用稳定排序处理 ties，再进入 Kokkos VBEM；该顺序会影响随机协方差初始化。pinned
> GATK 4.6.2.0 oracle 双后端验证 MQ,QD 顺序与 VQSLOD 语义一致；完整多高斯 raw-bit
> VBEM/resource calibration 仍为 explicit fallback。

> **LearnReadOrientationModel utility CLI（2026-09-01）**：native 支持 GATK
> `--QUIET/--quiet`、`--tmp-dir`、`--verbosity/--VERBOSITY` 与 JDK codec Boolean
> aliases，并写入 manifest telemetry；两样本 orientation-prior tar writer pinned
> oracle 在 OpenMP/Serial 均通过，gzip/member bytes 仍明确非 bit-identical。另已对齐
> `--max-depth`：保留输入直方图的高深度 bins，仅限制 EM 观测循环；双样本 max-depth=2
> pinned oracle 与 GATK 最大概率差约 1.9e-12。

> **CollectReadCounts CNV utility CLI（2026-09-01）**：native 补齐 GATK
> `--QUIET`、`--read-validation-stringency/-VS`、`--disable-bam-index-caching/-DBIC`、
> `--tmp-dir` 和 JDK codec Boolean/value controls；既有 pinned Java/native TSV
> interval/exclusion oracle 在 OpenMP/Serial 均通过，choices 写入 telemetry。

> **GatherBQSRReports utility CLI（2026-09-01）**：native 补齐 GATK
> `--QUIET/--quiet`、`--tmp-dir`、`--verbosity` 与 JDK deflater/inflater Boolean
> forms，并将选择写入 manifest telemetry；既有 pinned GATK 4.6.2.0 五表 oracle
> 以这些 controls 运行仍 bit-identical，空报告/不兼容维度继续 fail-closed。

> **SortSam Picard CLI boundary（2026-09-01）**：native 补齐大写
> `--INPUT/--OUTPUT/--REFERENCE_SEQUENCE`、`-SO/--SORT_ORDER`、
> `--MAX_RECORDS_IN_RAM/--TMP_DIR`、`--CREATE_INDEX`、`--QUIET`、codec 开关、
> `--COMPRESSION_LEVEL`、`--VALIDATION_STRINGENCY` 与 `--VERBOSITY`；pinned
> GATK 4.6.2.0 oracle 在 OpenMP/Serial 比较排序 header/records 和 manifest writer
> telemetry。`duplicate` order 现按 HTSJDK `SAMRecordDuplicateComparator` 的 library、
> unclipped read/MC mate coordinate、orientation、mapped-end、reference-length score、
> name/end keys 做 external-memory merge，并由独立 pinned oracle 在双 backend 比较
> 完整记录顺序；cloud/MD5/byte identity 仍是显式 fallback。

> **GatherVcfs CLI boundary（2026-09-01）**：native 补齐 Picard/GATK
> `--INPUT/--OUTPUT/--REFERENCE_SEQUENCE/--QUIET/--VERBOSITY` 大写别名，并由 pinned
> GATK 4.6.2.0 oracle 验证 comment、compression level=0、CREATE_INDEX=false 与
> 两 shard 数据行一致；完整 cloud/TMP/MD5 仍是显式 fallback。

> **LeftAlignAndTrimVariants CLI/interval boundary（2026-09-01）**：native 新增
> `-XL/--exclude-intervals`、`-ip`、`-ixp`，按 END span 实现 include padding 与
> exclusion padding，并支持 `-isr/-imr/-OVI` aliases；`dont-trim/split/keep-original`
> 现接受 separated Boolean。pinned GATK 4.6.2.0 CLI oracle 在 OpenMP/Serial 均通过，
> manifest 记录 selector/padding/skip/index 状态。复杂 symbolic/cloud/bit-identical
> normalization 仍为 explicit fallback。

> **AnnotateIntervals interval validation（2026-09-01）**：native 对齐 GATK
> 4.6.2.0 的 CNV 参数校验：`OVERLAPPING_ONLY` merging、`UNION` interval-set，
> 且 interval/exclusion padding 必须为 0；include/exclude 选择及非零 padding 的
> fail-closed Java 诊断在双 backend oracle 中覆盖。完整 cloud 与 release-specific
> CNV 语义仍保持显式 fallback。

> **ApplyVQSR exclusion interval boundary（2026-09-01）**：native 支持重复
> `-XL/--exclude-intervals` 与 `-ixp/--interval-exclusion-padding`，在 `-L` include
> 之后以记录 span 做 exclusion union；pinned GATK 4.6.2.0 oracle 覆盖 include+exclude、
> exclusion-only、padding 三种路径，VQSLOD/FILTER/site rows 与 manifest 统计在 OpenMP/
> Serial 均通过。完整 VQSR model/recalibration-table provenance 仍为 explicit fallback。

> **HaplotypeCaller sites-only writer（2026-09-01）**：native 新增 GATK
> `--sites-only-vcf-output` optional Boolean。VCF/GVCF 的 FORMAT/sample 证据在
> genotyping 与 reference-confidence 阶段保持完整，最终 writer 才发布 8 列
> site-only 记录；pinned GATK 4.6.2.0 true/false/bare oracle 在 OpenMP/Serial
> 均通过。GVCF block 分段差异不被此 writer slice 隐藏，concrete candidate rows
> 逐字段一致。

> **HaplotypeCaller GVCF floor-blocks（2026-09-01）**：native 新增 GATK optional
> Boolean `--floor-blocks`。普通 `REF/<NON_REF>` reference block 在 RCM/genotyping
> 后将 GQ 向下取到当前 `-GQB` 分箱下界，并发布 `GT:DP:GQ`（不含 `MIN_DP/PL`）；
> candidate rows 和 BP_RESOLUTION 不受影响。扩展的 pinned 4.6.2.0 GQ-band oracle
> 比较 Java/native rows、header 与 manifest，OpenMP/Serial 均通过；完整 RCM 仍为
> explicit fallback。

> **Mutect2 normal-lod boundary（2026-09-01）**：native 支持 GATK
> `--normal-lod`（默认 2.2），在 somatic likelihood 物化后按 matched-normal ALT
> log10-odds 控制候选发射；pinned GATK 4.6.2.0 low-threshold/zero-threshold
> oracle 在 OpenMP/Serial 均通过，阈值写入 stats/manifest。完整 joint posterior
> calibration 仍是 fallback。

> **GenotypeGVCFs exclusion interval boundary（2026-09-01）**：native 现在支持 GATK
> `-XL/--exclude-intervals`，在 include interval set 之后以记录 span 对 exclusion union
> 做过滤；aggregate 与 `--stream-by-locus` 共用 Host predicate，排除-only 调用仍完整
> 扫描输入。pinned GATK 4.6.2.0 concrete-variant oracle 比较 Java/native 数据行和
> manifest telemetry，OpenMP/Serial 均通过；完整 cohort posterior calibration 仍是
> fallback。

> **GenotypeGVCFs STARTS_IN writer boundary（2026-09-01）**：native 新增 GATK
> `--only-output-calls-starting-in-intervals` optional Boolean。`-L` 的输入 traversal
> 仍按 span overlap 保留跨边界 deletion 参与 joint genotyping，最终 writer 才按
> VariantContext 起点过滤；aggregate 与 `--stream-by-locus` 共用同一 Host predicate。
> 自包含 pinned GATK 4.6.2.0 oracle 覆盖 interval 外起点、interval 内 SNP、bare/true
> 语义及 manifest 计数，并在 OpenMP/Serial 比较 Java/native 数据行。

> **ValidateVariants symbolic-ALT boundary（2026-09-01）**：native 对齐 GATK
> `VariantContext.validateAlternateAlleles()`，只对 concrete ALT 做 called-genotype
> usage 检查；`<DEL>`/`<CNV>` 等 symbolic ALT 在 0/0 或 no-call 记录中不再误报。
> pinned GATK 4.6.2.0 mixed symbolic/concrete oracle 与双 backend CTest 已加入，具体
> ALT 的严格失败语义保持不变。

> **VariantFiltration genotype no-call boundary（2026-09-01）**：native 新增
> `--set-filtered-genotype-to-no-call`，将带非-PASS FORMAT/FT 的 called genotype 改为
> unphased `./.`，仅在发生变化时重算 INFO/AC、INFO/AN、INFO/AF，并保留既有 FT 历史。
> pinned GATK 4.6.2.0 oracle 覆盖 bare/true/false、existing FT 与 embedded `=`
> fail-closed 边界，双 backend CTest 已接入。

> **LeftAlignAndTrimVariants sites-only writer（2026-09-01）**：native 新增
> `--sites-only-vcf-output` optional Boolean，并在 normalization/split/PL-AD 投影
> 完成后才用独立 writer header + `bcf_subset_format` 移除 FORMAT/sample，生成 GATK
> 4.6.2.0 一致的 8 列 VCF。pinned Java/native oracle 覆盖 bare、separated true/false
> 形式，OpenMP/Serial 均通过；完整复杂 symbolic normalization 仍是 fallback。

> **SelectVariants sites-only writer（2026-09-01）**：native 新增 GATK 4.6.2.0
> `--sites-only-vcf-output` optional Boolean；processing header 保留完整样本和 FORMAT
> payload，最终 writer header 才执行 `bcf_subset_format`，输出 8 列 site-only VCF。
> pinned `fastgatk-select-variants-sites-only-gatk-oracle` 覆盖 bare、分离 true/false
> 及 embedded `=` fail-closed，OpenMP/Serial 均通过；完整 JEXL/pedigree/random/cloud
> 语义仍保持 explicit fallback。

> **BQSR covariate aliases（2026-09-01）**：BaseRecalibrator 现在接受 GATK
> `RecalibrationArgumentCollection` 的 `-ics`/`-mcs` 短别名，并与长参数共享同一
> C++/Kokkos report 路径；pinned BI/BD indel oracle 已覆盖 `-ics 4` 的四张表。

> **ApplyBQSR report alias（2026-09-01）**：ApplyBQSR 现在接受 GATK
> 4.6.2.0 的 `-bqsr` 短别名，并与 `--bqsr-recal-file` 共享同一 report 选择路径。
> `fastgatk-apply-bqsr-alias-gatk-oracle` 比较 Java/native 的长短拼写和 493 条
> Java 解码记录，OpenMP/Serial 均通过；仅声明 CLI 直接替换边界，完整 BQSR 模型与
> cloud/MD5 仍保持 explicit fallback。

> **CallCopyRatioSegments z-score CLI boundary（2026-09-01）**：native 已接受 GATK
> 4.6.2.0 的完整参数名 `--outlier-neutral-segment-copy-ratio-z-score-threshold`
> 与 `--calling-copy-ratio-z-score-threshold`，同时保留原短别名。pinned
> Java/native 输出与 `.igv.seg` sidecar 在 OpenMP/Serial 均通过。

> **CallCopyRatioSegments non-finite input boundary（2026-09-01）**：GATK 的
> `CopyRatioSegment` TSV decoder 接受 `NaN`、`Infinity`、`-Infinity` 的
> `MEAN_LOG2_COPY_RATIO`；native 现让这些 IEEE-754 值穿过 Kokkos 变换/统计路径，
> 输出 neutral call，并在 called/`.igv.seg` 表中保留 Java 的文本表示。新增 pinned
> `fastgatk-call-copy-ratio-segments-nonfinite-gatk-oracle`，三种值的输出 bytes、
> sidecar 与 NaN 统计 manifest 均在 OpenMP/Serial 对齐；完整 ModelSegments 概率
> 分段仍是 fallback。

> **CallCopyRatioSegments compensated-sum boundary（2026-09-01）**：两次长度加权
> 统计现与 Java 17 `DoubleStream.sum()` 一样，在 Kahan 高位累加后减去其负补偿字，
> 并保留同号 infinity fallback。此前丢弃低位补偿会让一个精确 calling-z 边界由
> neutral 翻为 amplification；新增 pinned byte-exact oracle 锁定 called 与 `.igv.seg`
> 输出，OpenMP/Serial 均通过。该切片不宣称完整 ModelSegments posterior 等价。

> **VariantFiltration missing-value Boolean boundary（2026-09-01）**：
> `--missing-values-evaluate-as-failing` 现在遵循 GATK 的 bare/separated Boolean
> 解析并拒绝 embedded `=`, 缺失 INFO 与 `== null` JEXL 语义由 pinned Java/native
> oracle 覆盖，OpenMP/Serial 均通过。

> **DenoiseReadCounts 双输出 writer 边界（2026-09-01）**：GATK 4.6.2.0 要求
> 同时提供 `--denoised-copy-ratios`（或 `-O/--output`）和
> `--standardized-copy-ratios`。native 已拒绝单输出调用，避免漏发标准化副本；
> pinned Java/native missing-output guard 与 OpenMP/Serial 回归已通过。

> **DenoiseReadCounts HDF5 metadata boundary（2026-09-01）**：SimpleCountCollection 的
> HDF5 输入把完整 SAM sequence dictionary 存在 `/locatable_metadata/sequence_dictionary`，
> GATK 会将全部 `@SQ`（含 AS/M5/UR/SP）复制到 standardized 与 denoised TSV。native
> 现已在 HDF5 解码后保留该 header；pinned GATK 4.6.2.0 oracle 比较两份输出的 header
> 与数据行，OpenMP/Serial 均通过。HDF5 其余 release-specific metadata、Spark/cloud
> 语义仍按 registry 保持 explicit fallback。

> **MarkDuplicates CREATE_INDEX 默认边界（2026-09-01）**：GATK/Picard 4.6.2.0
> 默认 `CREATE_INDEX=false`。native 已修正此前默认生成 `.bam.bai` 的行为；显式
> `--create-output-bam-index=true` 才发布索引。默认/显式索引的 pinned
> Java/native 断言与 contract/oracle 在 OpenMP/Serial 均通过。

> **MarkDuplicates OpticalOnly tagging boundary（2026-09-01）**：GATK/Picard
> `TAGGING_POLICY=OpticalOnly` 只为 optical duplicate 写 `DT:Z:SQ`，普通 library
> duplicate 不写 `DT`。native 已将该策略接入既有 optical classification 与 writer
> 路径（不改变 duplicate flags/metrics）；新增 pinned Java/native 14-record
> oracle 与独立 CTest，OpenMP/Serial 均通过。`DontTag|All` 继续保持原语义，UMI/flow
> 和 cloud 语义仍 explicit fallback。

> **MarkDuplicates paired ReadEnds key（2026-09-01）**：native 首遍按 read name
> 配对 mapped mates，并使用两条实际 unclipped 5′ 坐标构造统一 duplicate key，修复
> 无 `MC` 标签或 terminal S/H clipping 时同一 pair 被拆开的差异；必要 key metadata
> 已纳入 spill checkpoint/resume。pinned GATK/Picard 4.6.2.0 4-record oracle 与
> OpenMP/Serial CTest 均通过。

> **GatherBQSRReports empty-report boundary（2026-09-01）**：native 现在复现 GATK
> 4.6.2.0 `RecalibrationReport.gatherReports()` 的可用性判定：所有输入若仅含
> Arguments/Quantized、RecalTable0/1/2 均无数据，则以 `there is no usable data in
> any input file` fail-closed，且不生成 primary/sidecar。pinned Java/native oracle
> 与现有 CTest 在 OpenMP/Serial 均通过。AnalyzeCovariates 当前 pinned CSV 已
> bit-identical；完整 PDF/R 绘图仍不宣称等价。

> **ReblockGVCF overlap trim/split boundary（2026-09-01）**：native 现在在 GQ-band
> merge 前复现 Java `ReblockingGVCFBlockCombiner` 的 buffered writer 边界：concrete
> variant/deletion 覆盖 hom-ref block 时裁剪或拆分 block，拆分 tail 用 indexed
> reference 重写 REF/POS，并在 OutputManifest 记录 trim/split/drop telemetry。新增
> `verify_reblock_overlap_gatk_oracle.py` 的 pinned GATK 4.6.2.0 fixture 在 OpenMP/Serial
> 均 exact（69000-69004 / 69005 / 69006-69010），CTest 和 `verify_all.sh` 已接入；
> 完整 annotation/posterior/deletion 语义仍是 fallback。

> **GatherPileupSummaries all-empty writer boundary（2026-09-01）**：native 已与
> GATK 4.6.2.0 对齐：先移除空 shard，再只对非空输入校验 SAMPLE；所有输入为空
> 且 SAMPLE 冲突时成功输出仅含列名的 table header，不写 SAMPLE metadata。新增
> pinned Java/native oracle 在 OpenMP/Serial 均逐字节通过；cloud 与完整 Java
> dictionary edge semantics 仍保持 fallback。

> **Collect/Gather Barclay short-name aliases（2026-09-01）**：GATK 4.6.2.0 的
> `CollectF1R2Counts` 与 `GatherPileupSummaries` 实际接受大写长别名 `--I`/`--O`；
> native parser 与 dispatcher registry 已补齐，同时保留 `-I`/`-O` 和 lower-case
> aliases。Collect archive 与 Gather table contract/oracle 在 OpenMP/Serial 均复核
> alias 输出，确保 Nextflow 生成的 GATK 参数可直接替换。

> **VariantRecalibrator sampling provenance（2026-09-01）**：native 已对齐 GATK
> 隐藏参数 `--sample-every-Nth-variant` 的零基 traversal phase（首条记录保留），
> 且 writer 只物化 sampled recalibration-table records；新增 pinned GATK
> `verify_variant_recalibrator_sample_every_gatk_oracle.py` 在 OpenMP/Serial 固定
> `[1,3,5]` 坐标、score/label/8-column shape 与 manifest count。此边界只收敛
> sampling/scatter provenance；full-covariance 多高斯 VBEM、完整模型收敛与 raw-bit
> 等价仍保持 fallback。

> **BaseRecalibrator mismatch-context boundary（2026-09-01）**：native 现在接受
> GATK 4.6.2.0 `--mismatches-context-size` 的 1..13 范围，并在 substitution
> ContextCovariate、GATKReport Arguments、checkpoint 和 manifest 中保持一致。
> pinned Java/native read-filter oracle 的 size=3 场景在 OpenMP/Serial 均逐表通过；
> full BQSR model、large corpus 与 cloud I/O 仍保持 explicit fallback。

> **VariantRecalibrator full-covariance default（2026-09-01）**：GATK 4.6.2.0 的
> `MultivariateGaussian` 始终使用完整协方差矩阵；native 已将默认模型从 diagonal
> 修正为 full covariance，旧 `--full-covariance` 作为兼容拼写保留。pinned VBEM
> score、Java model-artifact load 和 broad contract 在 OpenMP/Serial 均通过；完整
> 多高斯 raw-bit VBEM、AS model provenance 与大样本收敛不因此宣称完成。

本文把已有的 PairHMM、四类代表模块验证和 GATK 源码调研，收敛成可以逐步编码、
逐步替换 GATK 的工程方案。目标不是把 Java 对象逐个翻译成 C++，而是保持
GATK 的调用/文件/结果契约，同时把适合批量计算的部分重构为 Kokkos 数据流。

> **LearnReadOrientationModel multi-sample writer（2026-09-01）**：标准
> CollectF1R2Counts tar 可包含多个样本；native 现在按样本隔离 histogram/alt table
> 聚合并输出每个 `<sample>.orientation_priors` 成员。pinned GATK 4.6.2.0 双样本
> oracle、OpenMP/Serial build 均通过；仅 tar 成员路径前缀/顺序与 gzip 字节不保证一致。

> **第三波边界（2026-09-01）**：PairHMM fragment-first 聚合使用 Kokkos 有序
> fragment×haplotype 求和并保留 `-Infinity`，六 cell pinned GATK oracle 在两个
> backend raw-bit exact。ModelSegments manifest 分离 GATK total/burn-in/retained
> 与 native conditional-point estimate，并标记完整 Gibbs/MCMC/RNG parity=false。
> GenomicsDBImport native workspace 明确为 `fastgatk-portable-sparse-v1`，不是
> TileDB、不可由 GenomicsDB query；真实 TileDB 仍由 external adapter/legacy bridge
> 提供。相关双后端边界测试均通过。

> **当前状态快照（2026-09-01）**：本文件包含按日期追加的历史增量；当前可复核的
> registry 为 48 项（31 `prototype`、16 `contract-compatible`、1 `adapter`），主计划
> 与 [progress_score.json](</home/turing-agents/Documents/fast-gatk/progress_score.json>)
> 的加权快照为 `global_score=0.5385`（约 54%）。OpenMP/Serial 当前 CTest 清单各为
> **209** 项；最新新增 HaplotypeCaller sites-only/floor-blocks/CIGAR indel-activity writer、GenotypeGVCFs/ApplyVQSR interval exclusion、orphan spanning-deletion 与 STARTS_IN writer boundary、ModelSegments default KernelSegmenter/multisample vector kernel/first-oriented-ALT segmentation/count-backed AF likelihood、VariantRecalibrator annotation-order/zero-variance/retry-iteration guards、CallCopyRatioSegments non-finite IEEE-754/compensated-sum boundary、MarkDuplicates OpticalOnly/pair-key tagging boundaries、DepthOfCoverage default read-filter/deletion-site boundaries、LeftAlignAndTrimVariants sites-only/CLI boundary、Mutect2 force-active/normal-lod activity/emission/PCR indel-quality 与 NLOD boundaries、ValidateVariants symbolic-ALT validation、GatherVcfs CLI aliases、CollectReadCounts/CollectF1R2Counts exclusion、CollectAllelicCounts read-filter、SortSam CLI aliases/duplicate comparator、BaseRecalibrator/ApplyBQSR aliases/context-size covariate、CallCopyRatioSegments full z-score aliases、VariantFiltration missing-value/set-no-call、ApplyVQSR sites-only writer、DenoiseReadCounts HDF5 metadata、LearnReadOrientationModel multi-sample tar writer、GatherBQSRReports CLI/empty-report、FilterMutectCalls sample-aware、SelectVariants sites-only writer、AnalyzeCovariates -bqsr alias、VariantEval keep-ac0 aggregate/ValidationReport boundaries、FastaAlternateReferenceMaker homozygous-IUPAC selection、GenomicsDBImport native interval-span workspace、PairHMM fragment aggregation、GenomicsDB native-storage capability、workflow-local Nextflow/SLURM 闭环
> contamination-table/GermlineFilter、CombineGVCFs interval/reference-block 与
> VariantRecalibrator `--sample-every-Nth-variant` oracle。各项只收敛其明确的
> writer/过滤/interval/sampling 边界，不改变完整 VQSR 模型 fallback 边界。
> **此前增量（2026-09-01）**：CTest 清单曾为 159 项；在此前 140/140 连续基线之外，PoN v7 sample-metadata、HC RCM PL-range、GenotypeGVCFs GP/PG、ModelSegments allele-fraction initializer、ApplyVQSR default cutoff、GenomicsDBImport sample-map、CallCopyRatioSegments interval validation、CombineGVCFs interval/reference-block、LearnReadOrientationModel multi-sample writer 与 FilterMutectCalls sample-aware contamination-table/GermlineFilter 边界已在 OpenMP/Serial 双 backend 通过。新增边界共 20 项，均独立计证，不将配置数量误报为连续 159/159；当前清单已由本文件顶部的 165 项快照覆盖。
> 当前快照补充：CTest 配置已增至 152 项；在此前 140/140 连续基线之外，CombineGVCFs PL-less、FilterMutectCalls orientation-joint/normal-artifact、CreateReadCountPanelOfNormals 退化 SVD 与 PoN v7 sample-metadata、GenotypeGVCFs maximum-ALT、ModelSegments smoothing、Mutect2 normal-evidence、HC AssemblyRegion boundary/RCM PL-range、Mutect2 global-mismapping-rate 与 DenoiseReadCounts integer-input 共十二项 pinned oracle 均已双 backend 定向通过。修复 GenotypeGVCFs dense-reference INFO/DP writer 边界后，按同一版本分段执行的当前完整清单合计 140/140：#1–#59 各 59/59、#60–#67 各 8/8、#68–#140 各 73/73（OpenMP 总 2213.42 s、Serial 总 2266.39 s）；新增边界测试均独立计证，不把 152 项配置误报为连续 152/152。此前 134/134 是修复前的历史整套基线；本次新增的 ModelSegments supplied-partition、CopyRatio/allele-fraction 条件模型与 FilterMutectCalls `-OVI`/contamination-joint 测试均在双 backend 定向通过；本次新增的 GenotypeGVCFs legacy
> AF-calculator oracle 双 backend 定向均 1/1；新增 `fastgatk-resource-limits-contract` 与
> `fastgatk-slurm-resource-wrapper-contract` 已在两个 backend 各定向通过，
> `fastgatk-progress-score-contract`、`fastgatk-fixture-digest-contract` 和
> CountBasesInReference contract 均已定向通过。历史段落中的 79/79、81/81、82/82、
> 83/83 等数字只代表各自日期的测试
> 快照，不应覆盖当前清单或被相加。

> **GenomicsDBImport sample-name-map 增量（2026-09-01）**：native sparse workspace
> 现在按 GATK 约定校验 sample-map 的两列格式和唯一 sample 名；对于本地单样本
> VCF/GVCF，只在 `native-inputs/` 私有物化副本中重写 header，保留源文件不变，
> 使后续 native GenotypeGVCFs 的 sample column 与 GATK `callset.json` 一致。
> `verify_genomicsdb_import_sample_map_gatk_oracle.py` 已在 OpenMP/Serial 通过；
> 该边界不宣称 sparse workspace 是 TileDB。

> **GenomicsDBImport 增量 workspace（2026-09-01）**：外部 adapter 现在透传
> `--genomicsdb-update-workspace-path`，要求已有 workspace，保留完整输入索引并验证
> GATK 4.6.2.0 的初始 interval set 在 update 时继续生效；native sparse workspace
> 对本地输入实现了 deterministic atomic input-index rebuild，保留旧/新输入顺序与
> sample-map，manifest 标记 non-TileDB。新增 update-workspace pinned oracle，
> OpenMP/Serial 均通过；完整 TileDB append/storage 仍由 GATK backend 负责。

> **GenomicsDBImport native interval-span 增量（2026-09-01）**：`--fastgatk-native-workspace`
> 现在按 GATK 4.6.2.0 的记录 span 语义处理 `-L/--intervals`（含 interval-list/BED
> selector）：记录与 selector 相交即保留完整记录，REF confidence block 不在 writer
> 边界被截断；重复 selector 支持 UNION/INTERSECTION。过滤后的物化输入不会复用原始
> VCF/CSI index，避免 stale chunk 造成 query 漏数；native checkpoint 同时固定
> interval selector/set-rule，增量 update 复用初始 interval 集。新增
> `verify_genomicsdb_import_native_interval_gatk_oracle.py`，用 pinned GATK workspace
> bridge 比较 block/variant 两个边界，OpenMP/Serial 均通过。完整 TileDB storage/query
> 仍由 GATK backend/legacy bridge 负责。

> **ReblockGVCF RGQ 小数阈值（2026-09-01）**：`--rgq-threshold`/`--rgq-threshold-to-no-call`
> 按 GATK Java double 语义解析，`10.5` 已由 pinned Java/native oracle 在双 backend
> 验证；完整 GVCF reblocking 仍是 prototype。

> **ReblockGVCF 当前边界补充（2026-08-31）**：新增真实 GATK 4.6.2.0 shard
> oracle，验证两个同一样本、非重叠 `.idx` shard 的 block/ALT/PL/AD/GT/DP/GQ
> 核心字段在 OpenMP/Serial 均 `1/1` exact，并固定 GATK 对多样本输入的拒绝语义。
> native 同时修正 `--rgq-threshold-to-no-call VALUE` 长别名的 separated/inline
> 解析，以及显式 RGQ/PL[0] 阈值不依赖 `--drop-low-quals` 的 low-quality deletion
> 转 GQ0 reference block 行为，并补齐 high-quality variant 的 `RAW_MQandDP`/
> legacy `RAW_MQ`→`MQ_DP` annotation parity；相关 contract 与 shard oracle 均通过。完整
> annotation/posterior/deletion-overlap 语义仍未达到直接替换门槛；同时补齐 GATK 正式
> CLI 的 `--format-annotations-to-remove` 与 `--do-qual-score-approximation` 别名，
> 并使 `--floor-blocks` 按 GATK 移除 reference-block 的 `PL`/`MIN_DP`；
> `--drop-low-quals` 现在在合并前丢弃 GQ0/低于 RGQ 阈值的 block，避免跨空洞合并；native 多样本
> sample-major 路径仍仅为 prototype。concrete variant 的 drop-mode 低置信 ALT 丢弃与高置信
> hom-ref→REF/<NON_REF> PL/GQ 投影也已接入独立 Java oracle，OpenMP/Serial 均 exact；
> `ctest -R fastgatk-reblock-gvcf` 两个 backend 均 4/4。

> **FilterMutectCalls MMQ 默认值补充（2026-08-31）**：GATK 4.6.2.0 的
> `M2FiltersArgumentCollection` 原始字段虽以 `-1` 表示未显式设置，但 getter 在
> 正常模式解析为有效阈值 30、microbial mode 解析为 20；native 保持该有效默认，
> 不把 help 中的 sentinel 当作最终过滤值。正确 `MMQ` `Number=R`（`REF,ALT=60,29`）
> 的 Java/native contamination oracle 与 OpenMP/Serial contract 均通过；完整
> ErrorProbabilities joint model 仍是 fallback。

> **FilterMutectCalls orientation ROQ writer 增量（2026-09-01）**：新增
> pinned GATK 4.6.2.0 oracle，固定 `QualityUtils.errorProbToQual()` 的
> `[1,93]` phred 边界。native orientation posterior 为 0/1 时不再写出
> `ROQ=1000/0`，而是复用 bounded quality helper；canonical/reverse-complement
> prior、`F1R2/F2R1 Number=R` 与 MNV fixture 的 GATK/native `ROQ=[93,93,1]`
> 和 orientation FILTER pattern exact，OpenMP/Serial 均通过。完整
> FilterMutectCalls joint posterior 仍是 fallback。

> **FilterMutectCalls `-OVI` CLI 别名增量（2026-09-01）**：GATK 4.6.2.0
> 将 `--create-output-variant-index` 注册为 `--create-output-variant-index,-OVI
> <Boolean>`；native parser、dispatcher registry 和 help 均已补齐 `-OVI`
> separated Boolean。pinned Java/native oracle 对 `-OVI false`/`true` 的 `.tbi`
> 缺失/存在边界在 OpenMP/Serial 均通过。GATK 不接受 inline `-OVI=true`，因此
> 不将该扩展形式列入 direct-replacement 语法契约。

> **FilterMutectCalls contamination joint-learning 增量（2026-09-01）**：经验
> `OPTIMAL_F_SCORE` pre-pass 现在通过统一 Kokkos AD/POPAF batch helper 计算
> `ContaminationFilter` 后验，并在 Java `ErrorProbabilities` 的 `NON_SOMATIC`
> 类型内与 GermlineFilter 取 max，再学习联合阈值。新增 pinned GATK 4.6.2.0
> `fastgatk-filter-mutect-calls-contamination-joint-oracle`，固定 AD 梯度上的
> contamination FILTER 集合、阈值范围和 posterior telemetry，OpenMP/Serial 均
> 通过；其余 release-specific somatic calibration 仍保持 fallback。

> **FilterMutectCalls orientation joint-learning 增量（2026-09-01）**：启用
> `--orientation-bias-artifact-priors` 时，Kokkos F1R2/F2R1 weighted-median
> `ReadOrientationFilter` 后验现在作为 `ARTIFACT` contributor 同时进入经验
> `OPTIMAL_F_SCORE` pre-pass 和最终 ALT 联合后验。新增
> `fastgatk-filter-mutect-calls-orientation-joint-oracle`，以 pinned GATK 4.6.2.0
> 验证 orientation FILTER pattern、阈值 0 边界及无 priors baseline；OpenMP/Serial
> 均通过。完整 release calibration 仍是 fallback。

> **FilterMutectCalls contamination-table sample fallback（2026-09-01）**：审计
> GATK 4.6.2.0 `ContaminationFilter` 确认 `--contamination-table` 是按 tumor
> sample 名称查找的 map；不匹配样本使用 `--contamination-estimate`（默认 0.0），
> 不能把其它样本的最大值当作全局污染率。native 已保留 per-sample lookup，表默认
> 回退为 0.0，并在显式 estimate 存在时仍加载表使匹配行覆盖 fallback。现有
> `fastgatk-filter-mutect-calls-contamination-oracle` 扩展 OTHER/TUMOR 不匹配
> fixture，Java/native contamination FILTER pattern 在 OpenMP/Serial 均通过；
> 完整 release-specific joint calibration 仍是 fallback。

> **GenotypeGVCFs posterior-QUAL 边界补充（2026-08-31）**：真实 Java/native
> oracle 现在同时运行 `--use-posteriors-to-calculate-qual`（`--gp-qual`）。对
> HaplotypeCaller 生成的仅含 PL/GQ、没有 FORMAT/GP 的 gVCF，三条 variant row
> 的 QUAL 与整行文本在 OpenMP/Serial 均 exact，native posterior kernel calls
> 与 samples 均为 0；只有显式 `USE_POSTERIOR_PROBABILITIES` assignment 才进入
> Kokkos posterior kernel，保持 GATK `hasPosteriors` 条件。完整 GP 输入、joint
> posterior 与 release-specific calibration 仍为 fallback。

> **Launcher contract 增量（2026-08-31）**：`fastgatk` registry 现在冻结
> `--java-options`/`-java-options`、`--gatk-config-file`、递归 `@args`（深度 16、
> `@@` 转义）、前置 `--` 分隔符和 unknown-option fail-closed 语义。JVM 参数在
> native 路径只作为 launcher metadata，不会进入 native binary；GATK properties
> 文件会自动转入 Java fallback，避免静默忽略 codec/cloud/tool 默认值。
> post-tool `--` 参数组在显式 fallback 时完整保留。新增
> `fastgatk-launcher-contract` 双 backend CTest，覆盖 native argv、fallback argv、
> 嵌套参数文件和未知参数阻断。

> **当前运行时/文件边界增量（2026-08-31）**：本轮 BQSR Java oracle 已把四张
> GATKReport 表与 493 条 Java-`PrintReads` 稳定记录（含 `OQ`/
> `--use-original-qualities`）纳入双 backend 回归；HC reverse-strand indel
> oracle 固定 reverse-complement、anchored CIGAR 与 left-normalization。runtime
> contract 现在同时验证 cgroup/SLURM hard/free host budget、visible-device
> memory、scratch/local-SSD、remote-input、byte-bounded backpressure 与
> `RESOURCE_EXHAUSTED` fail-closed 行为。远端输入在 native 前做原子 staging；
> 普通 VCF/Tribble `.idx` 由共享 writer 在 primary output 校验后原子发布。
> SelectVariants 的 Java 4.6.2.0 semantic oracle 覆盖 type/sample/JEXL/concordance/
> allele-subset，并校验稳定的 GT/AD/PL/GQ/AC/AN/AF payload；
> `fastgatk-runtime-output-smoke` 覆盖完整 primary/index/manifest 发布、缺索引、
> 临时残留 fail-closed 与失败提交回滚。上述兼容性边界聚焦回归在 OpenMP/Serial
> 均 `54/54`，不改变尚未完成的完整算法、GPU 实机与多节点发布门槛。

> **BQSR I/D event-table 增量（2026-09-01）**：新增
> `fastgatk-bqsr-indel-gatk-oracle`，用 pinned GATK 4.6.2.0 对含确定性
> per-base `BI`/`BD` Z tag 的同一批 reads 逐表比较 `Quantized`、`RecalTable0/1/2`。
> BaseRecalibrator 现在按 `(read-group, reported-quality, event)` 保留 I/D
> quality key，使用 Java RecalDatum 的 empirical/prior 规则生成 report；缺失
> tag 仍回退 Q45/显式 default，非空长度错误 fail-closed。OpenMP/Serial 均
> `1/1` 通过，同时重编译 `fastgatk-bqsr` 与 `fastgatk-apply-bqsr`；该 bounded
> slice 不宣称完整大样本、cloud 或全量 BQSR 模型等价。

> **ApplyBQSR preserve-threshold 增量（2026-09-01）**：新增
> `fastgatk-bqsr-preserve-gatk-oracle`，以 pinned GATK 4.6.2.0 report 和混合
> Q0..Q40 输入质量验证非默认 `--preserve-qscores-less-than 10`。native 与 Java
> `PrintReads` 后的 493 条记录 exact，且低于 Q10 的质量逐碱基保持不变；C++
> decode/encode 共用阈值谓词，确保低质量碱基绕过 Bayesian model 与 quantizer。
> OpenMP/Serial 各 `1/1` 通过；不宣称完整 BQSR 模型等价。

> **BaseRecalibrator read-filter 增量（2026-09-01）**：native 现在真实解析
> GATK 4.6.2.0 的七个默认 read filter（MappingQualityNotZero、
> MappingQualityAvailable、Mapped、NotSecondaryAlignment、NotDuplicate、
> PassesVendorQualityCheck、Wellformed），支持 `-DF/--disable-read-filter`、
> `-RF/--read-filter` 与 optional-boolean `--disable-tool-default-read-filters`。
> 解析顺序固定为“默认集合（或清空）→禁用默认→显式追加”，不支持的 filter 名称及
> inverted filter fail-closed。新增 `fastgatk-bqsr-read-filter-gatk-oracle`，用含
> MAPQ 0/255、secondary、duplicate、vendor-fail、supplementary 的合成 SAM，四种
> 配置与 pinned GATK 4.6.2.0 的 `Quantized`/`RecalTable0/1/2` 均逐表 exact，
> OpenMP/Serial 各 `1/1` 通过；该 bounded slice 仍不扩展为完整 BQSR 模型或 cloud
> 语义等价。

> **HaplotypeCaller k-mer-list 增量（2026-09-01）**：ReadThreadingAssembler 的
> `--kmer-size` 现在按 GATK 4.6.2.0 语义支持重复参数；默认列表为 `[10,25]`，Host
> 先排序去重后交给共享 Kokkos graph retry，旧 `graph_kmer_size` 标量仍记录最后一次
> 显式值。新增 `verify_hc_kmer_list_gatk_oracle.py` 与 CTest，在 OpenMP/Serial 对默认
> 和乱序 `[25,10]` 进行了 GATK Java VCF 字节级比较并检查归一化 telemetry。该 slice
> 只收敛 k-mer 参数/尝试顺序，不宣称完整 graph/path 或复杂 assembly 等价。

## 1. 设计结论

### 1.1 两个明确的边界

每一个 native tool 都拆成两层：

```text
C++ Host / Tool orchestration
  参数、对象语义、HTSlib/GenomicsDB、region 状态机、traceback、排序、writer
             |
             |  flat SoA batch + Kokkos::View + explicit plan
             v
Kokkos Kernel
  parallel_for/reduce/scan、TeamPolicy、ThreadVectorRange、Kokkos SIMD、Kokkos Kernels
```

Host 层拥有生命周期和兼容性；kernel 层只消费固定布局的数值 batch。kernel 不
调用 HTSlib、文件系统、网络、SLURM、Java 或 C++ 异常，也不保存指向 Host 对象
图的指针。

### 1.2 由实测数据得到的硬约束

| 已有证据 | 对方案的约束 |
|---|---|
| PairHMM AVX-512 steady-state kernel 为 GKL 实际矩阵吞吐的 1.09–1.28×，但 16 核摊入 prepare 后为 0.94× | 所有高频 kernel 必须有持久化 `Plan`、workspace pool、cached transition/packing；不能每次调用重新分配和 deep-copy |
| Smith-Waterman 规则 DP 从约 355M 到 4.94G cell updates/s（1→16 核） | 规则 DP 用 batch + `RangePolicy`/跨 work-unit SIMD，优先验证强扩展 |
| BQSR 细粒度 bin 原子在高核数退化 | 使用 per-thread/per-team local table，再做确定性整数 reduction；全局 atomic 只作为小批量 fallback |
| k-mer 图度数统计高核数不稳定 | graph ownership、hash/CSR/path traversal 保留 Host；Kokkos 只做编码、排序、计数、边评分和 compact |
| region 小批量的计算吞吐随核数下降，I/O 与 kernel 可分离 | reader/writer 独立 CPU pool，按 byte capacity 做 bounded queue，只有达到 batch 阈值才送 Kokkos |
| scalar/AVX2/AVX-512 三个变体结果一致；GPU 尚未实测 | 用同一 Kokkos 源码构建多产物；严格数值模式先固定 CPU oracle，GPU 先做 correctness backend 再优化 |

当前同一源码已通过 OpenMP 与 Serial Kokkos backend 的完整构建、kernel API
smoke 和 benchmark；CUDA/HIP/SYCL 仅在对应编译器和设备存在时启用，缺失时由
backend matrix verifier 明确报告 skip，不把 host 结果伪装成设备验证。

### 1.3 兼容模式和快速模式

```cpp
enum class DeterminismMode { Strict, Fast };
enum class Backend { Auto, OpenMP, Threads, Cuda, Hip, Sycl, GatkFallback };

struct ToolContext {
  ResourceProbe resources;
  TraversalPlan traversal;
  IoContext io;
  DeterminismMode determinism;
  Backend backend;
  TelemetrySink telemetry;
};
```

- `Strict`：固定输入顺序、bucket 顺序、归约树、排序 tie-break、随机种子和数学
  函数；PairHMM/SW 通过 GATK oracle 做 raw-bit 或逐 CIGAR 比较。
- `Fast`：允许 backend 选择更适合的 reduction/float，但每个工具定义 tolerance
  和生物学结果门槛，不能默默变成 Strict。
- `GatkFallback`：保持原始 GATK 命令、JVM 选项和输出；native 失败时由 launcher
  显式记录 fallback 原因。

## 2. 总体代码结构

建议的 C++ 仓库结构：

```text
cpp/
  compat/       gatk launcher、argument file、registry、fallback
  runtime/      ResourceProbe、AdaptiveController、queues、telemetry
  traversal/    interval/shard/halo、AssemblyRegion、deterministic order
  io/           HTSlib/BGZF/CRAM/VCF、reference、cloud range、spill
  model/        flat read/variant/haplotype/feature schemas、arena
  kernels/
    pairhmm/    regular/flow、CPU SIMD、GPU wavefront
    sw/         affine score、trace checkpoints、CIGAR host assembly
    kmer/       encode/sort/unique/edge primitives
    pileup/     activity、histogram、BAQ、filters
    bqsr/       covariates、local tables、merge
    genotype/   likelihood、genotype math、reference confidence
    dense/      Kokkos Kernels GEMM/SVD/GMM primitives
  tools/        HC、Mutect2、BQSR、GenotypeGVCFs、CNV、VQSR adapters
  formats/      output manifest、VCF/BAM metrics、schema/version
  tests/        oracle、format、tool、pipeline、resource tests
```

每个 kernel 模块统一提供四个对象：

```cpp
template<class ExecSpace>
class KernelPlan;  // 长寿命：bucket、View、workspace、scratch 规格

struct HostBatch;  // std::vector/arena 所有权，只在 Host
struct DeviceBatch; // Kokkos::View 所有权和 mirror，不能暴露内部指针

template<class ExecSpace>
PrepareResult prepare(KernelPlan<ExecSpace>&, const HostBatch&, ExecSpace);

template<class ExecSpace>
ComputeResult execute(KernelPlan<ExecSpace>&, const DeviceBatch&, ExecSpace,
                      DeterminismMode);

CollectResult collect(const ComputeResult&, HostBatch&, DeterminismMode);
```

`prepare` 和 `execute` 必须分开计时。生产 pipeline 重复使用 `KernelPlan`；
短 batch 可由 Host 直接执行 reference kernel，不强行发起设备拷贝。

## 3. 基础模块方案

### 3.1 `compat`：GATK 入口和工具注册

**C++ Host**

- 实现同名 `gatk` dispatcher，支持 `--list`、`--help`、`--version`、`--dry-run`、
  `--java-options`、`--gatk-config-file`、`@args` 和 `--` 分隔符。
- 解析 launcher 选项，但不吞掉 tool 参数；保留原始 argv 供 fallback。
- `tool_registry.json` 登记 native/fallback、支持参数、输出/索引/sidecar、
  backend、determinism 和 golden corpus。
- 选择顺序：参数/输出契约检查 → native backend 探测 → 资源估算 → native；
  backend 不可用或契约失败时 explicit fallback 或 fail-closed。
- 远程输入在 native dispatch 前通过有界并行、无 shell 的 argv-based curl staging 落到
  content-addressed cache，原子提交并按工具类型获取 FASTA/HTSlib sidecar；
  远程输出默认 fail-closed，显式 `FASTGATK_REMOTE_OUTPUT_MODE=upload` 时按
  sidecar→主文件顺序提交 HTTP(S)/presigned URI，fallback 保留原始 URI。

**Kokkos kernel**

- 无 kernel。该模块只是 Host 的兼容和安全边界。

**验收**

- 与 GATK launcher 的 exit code、帮助文本关键字段、参数错误和 `--dry-run`
  行为一致；Nextflow process 不改命令格式即可运行。

### 3.2 `runtime/traversal`：资源、interval、shard、halo

**C++ Host**

- `ResourceProbe` 读取 cgroup v2、SLURM CPU/memory/GRES、local scratch/inode、进程文件描述符上限、
  远端 I/O 延迟；不能用物理机总内存替代 job limit。
- `AdaptiveController` 根据 host/device/scratch/inflight bytes 调整 batch；
  只在 assembly region、read group、VCF/GVCF block 等安全边界缩放。
- `TraversalPlan` 保留 `core_interval`、`halo_interval`、`output_ordinal`、
  `shard_id`；halo 只计算不输出。
- Host 维护 contig 顺序、active-region 状态机、downsampling 和 deterministic
  read order。

**Kokkos kernel**

- 固定窗口的 activity/pileup、filter mask、base count、候选计数使用
  `RangePolicy` 或 `TeamPolicy`。
- 跨窗口的 profile 状态转移仍在 Host 顺序执行；kernel 不能通过锁等待其他 team。

**内存/调度**

- 三队列：decoded reads → compute batches → encoded results；以 byte capacity
  背压，不以 item 数量背压。
- 默认 host inflight 目标为可用 budget 的 60%，设备 60%，scratch 80%；连续压力
  窗口降低 batch，不能在 OOM 后无限重试。

补充：`ReadBatch::bytes()` 统一统计 flat payload，`HtsReader::set_batch_records()`
只在 batch 边界生效；CollectReadCounts/GetPileupSummaries 已接入
`AdaptiveController::next()`，在观测到 host byte pressure 后缩小后续 decode
批次，并把 initial/effective/reduction telemetry 写入 OutputManifest。这样不
改变输入顺序，也不把“按条数”误当成内存上限。

运行时现在还提供共享的 `ThreeStagePipeline<Decoded,Computed,Encoded>`：三个
字节有界队列连接独立的 decode/compute/encode Host stage，最终 sink 在调用线程
串行 drain，确保输出 ordinal 不因并发改变；任一 stage 异常会关闭全部队列、等待
线程退出并重新抛出，manifest 可记录每段 item/byte/peak occupancy。CountReads、
FlagStat、CollectReadCounts、GetPileupSummaries、CollectAllelicCounts、DepthOfCoverage
与 CollectF1R2Counts 已接入该执行器，作为所有读级工具复用的最小生产
示例；Kokkos
kernel 仍只在 compute callback 中执行，避免把 HTSlib/文件对象暴露给 device。
CollectF1R2Counts 也复用该执行器：Host 通过 HtsReader 完成 BAM/CRAM、CIGAR
和 F1R2 orientation observation 解码，再以 bounded observation slices 进入 Kokkos
`RangePolicy` 原子计数，最终 sink 合并到稳定 locus arrays；因此 Mutect2
的 F1R2 输入路径不会因全量 observation 物化而绕过统一背压契约。

ApplyBQSR 现同样复用该执行器：decode stage 保持 HTSlib BAM/CRAM 对象和
covariate/CIGAR projection，compute stage 只把扁平 quality/delta 数组交给 Kokkos
`RangePolicy`，encode stage 承担 preserve-Q/quantization，sink 在调用线程按输入
顺序写回并更新 checkpoint。BAM 对象通过 move-only RAII 跨队列传递，三段队列均
按 byte capacity 背压；manifest/benchmark 记录 item、byte、peak occupancy。
OpenMP/Serial 的 BQSR report、ApplyBQSR SAM 核心字段和 checkpoint oracle 均通过，
长读/大样本及 object-store 语义仍是后续工作项。

HC decode 现在也使用同一 `ReadBatch::bytes()` 和 batch-boundary setter。由于
AssemblyRegion/PairHMM 仍要求稳定的 aggregate HostBatch，当前只自适应下一批
读取，并对累计 staging 做饱和加法；`requested_batch_records`、
`effective_batch_records`、`adaptive_batch_reductions` 写入 summary/manifest。它
是内存压力控制的真实增量。另有显式 `--stream-by-contig` 路径：要求输入按
coordinate 排序，在 contig 边界释放 aggregate HostBatch，增量写出 VCF/GVCF
body，并记录 `streamed_contigs`/`streamed_peak_host_bytes`；该路径已通过双
contig normal-vs-stream exact contract，但仍不冒充按 AssemblyRegion 的全基因组
固定内存上限，后者仍是后续工作项。
本轮进一步把 `--stream-by-contig` 接入 `ThreeStagePipeline`：decode 线程只负责
坐标检查和 contig-boundary `ReadBatch`，compute 线程调用共享 assembly/PairHMM
Kokkos 路径，sink 在唯一写线程渲染并合并 `Result`。独立只读 header reader 避免
HTSlib 状态与并行 decode 竞争；manifest 记录三段 item/byte/peak telemetry，
OpenMP/Serial 的 contig 与 region contract 均通过。
本轮再补 `--stream-by-region N`：对有 `.bai/.crai` 的 BAM/CRAM 使用 HTSlib
indexed iterator，按固定 core tile 加 activity/genotyping halo 解码；`locus_intervals`
把 halo 保留给 activity/graph/PairHMM，而 `intervals` 仍只控制 core 输出，避免小
tile 边界漏掉变异。tile 完成后释放 reads/matrix，只将 core-domain VCF body 增量写出；manifest 记录
`streamed_regions`、`stream_indexed` 和峰值 Host bytes。单 tile 与 aggregate
路径逐字节一致，多 tile 已验证 VCF locus identity；多 tile GVCF 现在保留跨 tile
pending HomRefBlock，并按 GQ band、逐位点 depth、MIN_DP 和 PL 规则 stitch 后输出。
跨 tile 的完整 AssemblyRegion graph context、paired-read state 和 GVCF block 的
更广泛 GATK oracle 仍需继续补齐，不能宣称已经达到全基因组 bit-identical。
`benchmark_hc_broad.py` 已将该 region pipeline 纳入同一 fixture benchmark，并对
region-streaming 的重复运行做 SHA-256 稳定性检查；region 与 aggregate 的 SHA-256
差异作为结果字段报告，不把独立 tile 的局部 activity/graph 证据误报为 exact。benchmark
结果只用于当前执行空间/资源配置的回归，不代表完整 genome speedup。

**验收**

- scatter/gather 后 core 区域无重复；padding 首尾不越界；cgroup 限制、scratch
  不足和失败重试均有可读错误。

### 3.3 `io/model`：HTS 输入、flat record 和输出

**C++ Host**

- HTSlib/BGZF/CRAM/VCF/BCF/index、reference range 和 cloud URI 负责解码/编码；
  `ReadRecord`、`HaplotypeRecord`、`VariantRecord` 使用 flat arena + offsets，
  不把 Java 对象图搬到 device。
- CIGAR 解析、reference projection、header/program record、压缩等级和 index
  由 Host 保持 GATK/htslib 语义。
- reader/writer 在独立 CPU pool；每个 chunk 携带 CRC、input offset、output
  ordinal 和 schema version，支持 spill/resume。

**Kokkos kernel**

- 在已解码 batch 上做 base/quality transform、filter mask、CIGAR op 分类、
  BAQ/pileup primitive；`View` 用 SoA，byte/uint16 对齐。
- 不把压缩 block、文件指针、字符串或可变长 CIGAR 直接传入 device。

**验收**

- samtools/htslib read-back、header/contig/index、BGZF/CRAM writer、远端 range
  读取和 partial output 清理全部通过。

## 4. 计算模块逐项方案

### 4.1 PairHMM（P0，最高优先级）

**C++ Host**

- `PairHmmPlan` 长期持有 length buckets、SoA `HostBatch`、quality/transition
  table、三行 workspace pool、prefix-cache metadata。
- 先按 `(error_model, read_len, hap_len, lane_width)` 分桶；连续 haplotype
  保留 `hap_start_index` 和 permutation，完成后恢复 GATK matrix 顺序。
- Strict 使用 GATK quality table、(A+B)+C 运算顺序和 fdlibm log10；Fast 可选
  backend 优化但必须有 tolerance oracle。
- regular/flow PairHMM 使用不同 `KernelPlan` 和测试向量，不能用 regular kernel
  模拟 flow。

目前 regular 路径继续使用三行 M/I/D + Kokkos SIMD；flow 路径已独立为
`compute_kokkos_flow`：Host 按 `(read_flow_len, hap_flow_len)` 分桶，传入
FlowBasedRead 校准后的 `[read_flow][256]` probability table，device 使用
四步长 flow recurrence 和五行环形 scratch，输出稳定 request order。API smoke
和 kernel benchmark 已分别验证 `error_model=flow`、same-flow 优于 mismatch-flow
以及重复运行 checksum；原始 BAM `tp`/flow-order 解码、GATK FlowBasedRead
概率校准现由 Host `fastgatk::io::decode_flow_read` 按 `tp`/`t0`/boundary-flow
及 `clipProbs` 规则提供；`encode_flow_key` 已暴露与 GATK
`FlowBasedHaplotype` 一致的 key-to-base、reverse-key 和左右 clipping 映射。
HtsReader/ReadBatch 已直接保留原始 `tp`/`t0`/FO/mc aux 元数据，HC/Mutect2
在完整且一致的 flow batch 上自动选择该模型；GATK `sample.bam`/`sample.t0.bam`
key/matrix corpus 已作为 2 个 CTest oracle，另有 4 个固定 read/haplotype case
对 GATK `FlowBasedPairHMM` 做 Strict raw-bit likelihood oracle。基于 Host
read/reference span 与局部 haplotype-vs-reference SW CIGAR 的 flow haplotype
clipping 已接入 request-local 窗口并记录 clipping/fallback；完整 GATK
haplotype CIGAR 坐标映射、hmer/uncertainty 以及 alignment-derived flow clipping、GKL flow raw-bit
corpus 和 GPU anti-diagonal wavefront tiling 仍是后续完整语义；当前 flow recurrence 已统一通过 Kokkos `TeamPolicy` 启动，并在
Serial/OpenMP API smoke、benchmark 和 backend matrix 中记录 `TeamPolicy`。

regular PairHMM 的 Host 准备阶段已按 GATK 默认路径执行 MAPQ cap、Q18 以下
base-quality 归一到 Q6，并使用 `ReadUtils.DEFAULT_INSERTION_DELETION_QUAL=45`
和 GCP=10；因此不同 Kokkos execution space 收到相同的离散输入。SAM/BAM 的
BI/BD per-base indel-quality Z tag 已由 HTSlib Host 解码为 offset/payload，经过
过滤批次和 HC/Mutect2 batch aggregation 后进入 PairHMM；缺失 tag 仍按 GATK
缺失 BI/BD tag 的原始默认始终为 Q45；Conservative/Hostile/Aggressive 的 PCR
曲线只调整非末位位置（repeat length 0 的有效值为 Q40），末位仍为 Q45；
NONE 保持全数组 Q45，长度错误则 fail-closed。PairHMM base-quality floor 可通过
`--base-quality-score-threshold`（默认 Q18，范围 Q6..Q255）调整，
`--disable-cap-base-qualities-to-map-quality` 可关闭 MAPQ cap；HC/Mutect2
共享这两个 GATK 兼容控制项。
`--pcr-indel-model NONE|HOSTILE|AGGRESSIVE|CONSERVATIVE` 也在两者间共享：Host
按 GATK `ReadLikelihoodCalculationEngine.findTandemRepeatUnits` 的 8bp unit/20bp
repeat 上限和 rate-factor 曲线调整 BI/BD 后再进入 Kokkos；缺失 tag 先使用
Q45，默认 Conservative 的有效非末位起点为 Q40，末位保持 Q45，NONE 全程
Q45。model、rate factor、实际调整位置数均写入
summary/manifest/stats，避免把不同模型的结果混为同一输入。

**Kokkos kernel**

- CPU：`RangePolicy<ExecSpace>` 每个 work item 处理一个 SIMD pair group，
  `Kokkos::Experimental::simd<double>` 跨独立 pair；三行 M/I/D 原位更新，
  不改变单 lane recurrence。
- GPU：`TeamPolicy` 一个 team 处理一个 pair 或同长度小组；scratch 保存 tile，
  anti-diagonal/tiled wavefront 通过 team barrier，不做跨 team spin-lock。
- 输出 fixed-size likelihood/scaled-sum；Host 负责 permutation、log10 strict
  和最终对象映射（或按 backend 支持在 kernel 内执行相同 strict 函数）。

**内存/调度**

- prepare/deep-copy 不能在每次 invocation 重复；region/batch 生命周期内复用 plan。
- device 只在 pair 数量大、长度分布窄、拷贝可摊销时启用；小 batch 走 CPU。

**验收**

- GATK EXACT/ORIGINAL/LOGLESS/GKL 对照；4096 matrix raw-bit strict；empty/invalid
  输入；CPU scalar/AVX2/AVX-512 和 CUDA/HIP/SYCL correctness。
- 性能同时报告 kernel-only、prepare、amortized、端到端；不能只报告热 kernel。

### 4.2 Smith-Waterman（P0）

**C++ Host**

- `AlignmentRequest` 保存 ref/alt span、SW 参数、overhang、stable id；Host 维护
  traceback、CIGAR builder、start/end 和 GATK tie-break。
- batch 前按长度/参数/overhang 分桶；短序列可直接 CPU，长批量再设备化。
- 不把“最大 score 正确”当作完整兼容；输出必须包含 CIGAR、start/end。

**Kokkos kernel**

- score 阶段维护 affine `H/E/F` 滚动行；CPU 用跨 alignment SIMD 或
  `RangePolicy`，GPU 用 tile/anti-diagonal `TeamPolicy`。
- 当前 CPU 路径已用 `Kokkos::Experimental::simd<int>` 对同长度 alignment
  batch 做跨请求向量化；异构长度按稳定 chunk 使用相同 Kokkos scalar kernel，
  `simd_width/simd_groups` 进入 API telemetry 和 benchmark schema。
- kernel 输出 score、argmax、traceback checkpoint/压缩方向；变长 CIGAR 由 Host
  compaction 生成，避免 device 动态分配。

**验收**

- GATK Java/GKL 对照随机序列和真实 assembly region：score、CIGAR、overhang、
  tie-break、空序列/边界；然后比较 PairHMM 前置 assembly 流程。

### 4.3 k-mer/read-threading/local assembly（P0/P1，Host-first）

**C++ Host**

- 每个 assembly region 建立 NUMA-local flat arena；canonical k-mer node id、
  edge CSR、support、read/haplotype offsets 都由 Host 拥有。
- graph traversal、prune、reference-connected path、异常恢复、候选 haplotype
  生成和 debug dump 先留 Host，使用 stable sort 和明确 tie-break。
- 先实现 read error correction、k-mer size/threshold、reference padding 和
  homopolymer/flow 参数的 GATK 语义，再考虑设备。

**Kokkos kernel**

- `RangePolicy`/SIMD：2-bit rolling encode、reverse-complement、canonical key。
- `Kokkos::sort`/Kokkos Kernels：key sort、unique/run-length count、edge compact、
  support/coverage 计算；每 team local bins，避免实验中观察到的全局原子瓶颈。
- GPU 只传固定 edge arrays、support 和候选 score；复杂 path traversal 仍 Host。

当前 graph API 还接入了 GATK-shaped `min_dangling_branch_length`、
`min_dangling_matching_bases` 与
`recover_all_dangling_branches`：不回到 reference scaffold 的短分支在
bounded traversal endpoint 被确定性抑制；显式 recovery 只保留 `tid=-1`
的 unanchored Host path，并记录 branch path/base telemetry；现在对保留的
disconnected path 复用共享 SW 做 bounded reference recovery，达到匹配下限
（`min_dangling_matching_bases`，默认 3）后恢复 `(tid,start,end)`，其余噪声仍禁止直接进入坐标型 variant call。新增
SeqGraph 风格的非分支链压缩计数、read-edge multiplicity 支持和
support-descending/stable tie-break，并以 GATK 的 linear-chain、SNP
rethreading、双 ALT 分支及反向链 read corpus 做回归。现在还把物化
haplotype path 送入实际的 sequence-level SeqGraph rewrite，按 GATK 顺序执行
diamond 的共享 prefix/suffix split/merge、共同 suffix tail merge 和
linear-chain zip，重建并返回等价 path language；rewrite 节点/边和各类变换计数
通过 API/HC/Mutect2 telemetry 暴露。参考重复 k-mer 默认按 GATK 语义拒绝当前
k 并由 Host 递增重试，低层 API 仍可显式允许。全部 likelihood-based pruning
参数和真实大 assembly-region corpus 仍需继续收敛。

Assembly 前置还支持显式 `error_correct_reads`：共享 Kokkos
read-error-correction kernel 对 exact k-mer count 建立稀疏 k-mer 到最近
solid k-mer（最多两个 Hamming mismatch）的 correction map，再按重叠窗口做
strict consensus；这与 GATK NearbyKmerErrorCorrector 的核心语义一致，避免
真实低频 ALT 被单个窗口改写。修正后的 read 仅进入 assembly/graph，原始
filtered read 保留给 PairHMM/RCM；默认关闭，参数和 corrected/solid/
uncorrectable telemetry 写入 HC/Mutect2 manifest。低质量 tail clipping 与
long-homopolymer policy 仍需单独 oracle。另一个 GATK correction selector
`pileup_error_correction_log_odds` 已接入：Host 按 CIGAR/locus 构造
plurality pileup，复用 Mutect2 的 flat-Beta log-likelihood ratio 和
三处 mismatch indel guard，Kokkos 只执行最终 edit mask；有限阈值优先于
Nearby 路径，默认 `-inf` 关闭，并在 manifest 记录 mode/locus/skip telemetry。

Graph orchestration 还实现了 GATK ReadThreadingAssembler 的周期策略：对
read-supported directed cycle 按 `k+10` 递增重试，最大 `k=31` 且不超过本地
read/reference 长度；reference 重复 k-mer 也按 GATK 默认拒绝当前 k 并重试，
`dont_increase_kmer_sizes_for_cycles` 可显式关闭。每次尝试仍调用相同 Kokkos
graph/count API，Host 只负责有界重试和确定性选择，并暴露 selected k、iteration
count、cycle/rejection telemetry。graph 现在把物化 path 送入实际的 sequence-level
SeqGraph rewrite：linear-chain zip、diamond shared prefix/suffix split/merge、
共同 suffix tail merge，并重建等价 path language；`num_pruning_samples` 和
可选 adaptive pruning（initial error、LOD/seeding LOD、max-unpruned-variants）
也已接入 API/HC/Mutect2。linear-chain、diamond/tail、SNP rethreading、
MultiSampleEdge sample-gate 和 edge-support ordering corpus 已进入 API smoke；
linked-de-Bruijn raw-topology path 也已接入，并实现 uncovered pivotal-edge
artificial-haplotype recovery 与 disable toggle；AdaptiveChainPruner 的
constant-error likelihood 现在复用 Commons Math/GATK 的 digamma cutoff、
递归展开顺序和 binomialCoefficientLog 分支，并以 GATK 4.6.2.0 数值 oracle
覆盖小计数、digamma cutoff 和高深度路径；真实大 assembly-region corpus
仍需补齐。另已按 GATK
`AdaptiveChainPruner.getMaxWeightChain` 修正 adaptive seed：先比较链内最大
单 edge multiplicity，再比较链长；sum-vs-max 区分 corpus 已加入 API smoke。

**验收**

- GATK ReadThreadingGraph/SeqGraph test corpus；graph edge/path/CIGAR/haplotype
  顺序；同一 region 多次运行结果一致；host-first 与 device primitive 结果一致。

### 4.4 HaplotypeCaller orchestration（P0/P1）

#### 已落地的 HC fragment overlap 语义

在 read filter 之后增加 Host-side copy-on-write 阶段，严格要求 read name、
read-group、paired/mate flags 与 reciprocal mate coordinates 均可验证，再按
CIGAR 将两端投影到共同 reference coordinate。默认复现 GATK
`FragmentUtils.adjustQualsOfOverlappingPairedFragments`：相同碱基两端质量
封顶到 Phred 20，冲突碱基两端质量置零；结果质量数组被 assembly、graph、
SW/PairHMM 共同消费；reference-confidence 保留 GATK 的 filtered-region
evidence map 质量语义。`--do-not-correct-overlapping-
quality`（另有 descriptive alias）为显式禁用开关；缺失 fragment metadata 时 fail-closed，并
在 manifest/telemetry 中报告 metadata unavailable，而不是按坐标猜配。静态
双端 SAM 回归覆盖 concordant/conflicting overlap 计数和禁用开关。

PairHMM 进入 genotyping 前现在复刻 GATK
`ReadLikelihoodCalculationEngine.filterPoorlyModeledEvidence`：Host 按 trimmed
AssemblyRegion 的 HMM read length 计算
`min(2, ceil(length * expectedErrorRatePerBase)) * -4.0`，将低于阈值的 read
从 genotype likelihood/AD/DP 证据中移出，但保留给 assembly 与 RCM 的原始
filtered-region evidence map。默认 `expectedErrorRatePerBase=0.02`，兼容参数
`--expected-mismatch-rate-for-read-disqualification` 已接入 HC/Mutect2；阈值、
移出计数和 trimmed read length 进入 manifest。Host 同时按
AssemblyRegionTrimmer 的 SNP/INDEL genotyping padding 生成局部 PairHMM 窗口，
避免整个 activity halo 改变 disqualification 语义。VCF writer 的
`GT/GQ/DP/AD` 现在直接使用保留下来的 PairHMM evidence，PL 保持 prior-free，
genotype prior 只用于 GT 选择；bundled chr17 sentinel 已与 GATK 的
`1/1, GQ=3, DP=1, AD=0,1` 对齐，QUAL/PL 数值仍显式不是 bit-identical。

GVCF reference-confidence 分箱也已接入同一 Host/Kokkos 选项边界：重复的
`-GQB/--GVCF-GQ-bands` 采用 GATK 的“上界列表”语义（隐含下界 0，自动补齐
终止上界 100），并由 reference-block state machine 和 writer header 共享，避免
只改 header 而实际 block 仍按旧分箱合并。`10,20,50` 的
`0-10/10-20/20-50/50-100` header、记录和非法/乱序输入校验均已对 bundled
GATK 4.6.2.0 做 OpenMP/Serial oracle 回归。

**C++ Host**

- 复刻 `AssemblyRegionWalker`：读取 shard、建 active/inactive region、halo、
  read filters/downsampling、参考和 feature cache。
- `assembleReads` 负责调用 assembly Host；候选 haplotype 按 stable id 排序。
- 对每个 region 组装后调用 SW → PairHMM → genotyping；只有 core interval 写出。
- GVCF 模式保持 reference-confidence block 边界、header、index、PG record。
- Primitive CIGAR insertion 统一消费 BAM reverse flag：反向链的 SEQ 片段先
  reverse-complement 到 reference orientation，再进入 anchored ALT 与
  left-normalization；`fastgatk-hc-reverse-strand-indel-oracle` 在 OpenMP/Serial
  backend 独立验证该契约。

**Kokkos kernel**

- 复用 k-mer/SW/PairHMM/genotype kernels；不在 HC Host 内重复构造设备 batch。
- region 级 `PipelinePlan` 预估 reads × haplotypes、scratch 和输出大小，触发
  CPU/GPU backend。

**验收**

- HaplotypeCaller 4.6.2.0 integration tests、GVCF golden、VCF/header/index、
  GIAB truth/ confident region；native/fallback 逐参数对照。

### 4.5 Pileup/activity/BAQ/downsampling（P1）

**C++ Host**

- 维护按 coordinate 的 state machine、read filter 顺序、downsampling seed、
  BAQ 边界和 active-region profile；复杂 predicate 保持 Host。

**Kokkos kernel**

- 对固定窗口的 decoded reads 用 `RangePolicy`/`TeamPolicy` 做 base/quality/mapQ
  histogram、soft-clip mask、BAQ primitive 和 activity probability。
- per-team local histogram 后合并，避免低覆盖小窗口启动开销。

`GetPileupSummaries` 的 `GoodCigarReadFilter` 必须是 per-read predicate：Host
`ReadBatch::cigar_record_layout_valid(record)` 只检查当前 CIGAR span 与 read
consumption；`reference_end` 和 `project_read_offset` 也遵循同一边界。批级
`cigar_layout_valid()` 仍保留为严格 API/诊断检查，但不能作为 walker 的逐 read
过滤条件。zero-length CIGAR 与正常 read 混合在一个 batch 的 OpenMP/Serial
contract 覆盖该语义。

当 AF-bearing population records 存在但没有记录通过 AF/interval/biallelic-SNP
selection 时，`GetPileupSummaries` 按 GATK 成功写出 metadata + 六列表头的空表；
仅 AF header 缺失或全部记录缺 AF 才报告 `BAD_INPUT`。

**验收**

- pileup/count/filter/downsampling 随机和真实 shard；同 seed 严格一致；窗口边界
  与 halo 不产生重复或漏计。

当前 `compute_activity_profile_kokkos` 在有参考投影时接收每个位点的
`reference_bases`，在 Kokkos kernel 中计算 HC 的 non-reference activity；Mutect2
额外传递按 locus 分组的 base-quality 列表，在同一 Kokkos kernel 中复现
`Mutect2Engine.logLikelihoodRatio`（GATK digamma、flat-Beta entropy、Phred error
model 和多替换质量修正），再与 `initial-tumor-lod` 生成二值活动状态。无质量输入
继续使用显式 count/diversity fallback，避免 smoke/adapter 输入改变语义。

ActivityProfile/AssemblyRegion 控制现在也已统一接入 HC/Mutect2：
`--active-probability-threshold`、`--assembly-region-padding`、
`--max-assembly-region-size` 和 `--max-prob-propagation-distance` 分别控制
activity 门限、halo、active-span 硬上限和概率传播距离。Host state machine 不再把
最大区域长度当作 gap 阈值。Mutect2 按 GATK 保留独立的 raw AssemblyRegion，即使
padding/halo 重叠也不合并；HC 仍在共享 API 层使用有界 union。参数和最终 region
数量写入 manifest/telemetry，Kokkos API smoke 覆盖硬上限与传播距离分离。这样
graph/SW/PairHMM 的 scratch 预算可以由 AssemblyRegion 上限直接约束。

HC 的 bounded AssemblyRegion 兼容切片现已补齐：`pipeline.hpp` 默认
`assembly_region_padding=100`、`max_probability_propagation_distance=50` 与 GATK
4.6.2.0 一致；`calling_pipeline.cpp` 的 Host trim gate 计算原始
`unclippedReadLength`（减去 CIGAR soft-clip），仅要求裁剪结果非空，不再把短的
trimmed edge fragment 当作短读段丢弃。Pinned
`verify_hc_assembly_region_boundary_gatk_oracle.py` 在 max-size=50 下对照 Java
GATK，锁定 chr17 17:69368 的 `AD=20,22 DP=42` 和完整 VCF byte identity；
`fastgatk-hc-assembly-region-boundary-gatk-oracle` 在 OpenMP/Serial 均通过。
该边界改善跨 AssemblyRegion/shard 的 read ownership，但完整 HC graph、hmer 和
release-specific assembly 仍保持显式未完成。

### 4.6 BQSR / ApplyBQSR（P1）

**C++ Host**

- 解码 read、生成 covariate descriptor（read group、quality、cycle、context、
  event type）；管理 recalibration report/table schema、collapse/round 和输出
  quality 编码。
- BaseRecalibrator 在 Host decode 后固定执行 GATK 标准 BQSR read-filter 组合
  （MAPQ 非 0/255、mapped、非 secondary/duplicate/QC、Wellformed），并在
  covariate 计算前移除 CIGAR soft-clip；低质量尾端只在 ContextCovariate 中写成
  `N`，内部低质量碱基保持可用。
- `--known-sites` 按 GATK repeatable 参数收集多个 VCF/BCF，HTSlib Host 统一读取、
  按 BAM dictionary 映射，并按 VariantContext 的 REF/INFO-END 区间排序合并后做确定性
  union；checkpoint 签名与 manifest 记录完整输入集合。
- 每个 shard 产生整数 local table；按 GATK 固定顺序 merge，避免浮点 reduction。
- ApplyBQSR 使用 Host decode/encode/sink + Kokkos `RangePolicy` 的 bounded
  三阶段流式 transform；只有扁平 quality/delta batch 进入 device，BAM 对象不跨
  C ABI/设备边界。
- ApplyBQSR 的 `--use-original-qualities` 在 Host 解码 `OQ:Z` FASTQ 质量用于
  covariate lookup，保留 OQ tag 并对 malformed/长度不符输入 fail-closed；不把
  SAM tag 指针带入 Kokkos kernel。
- `--emit-original-quals` 在 Host encode 前只对缺失 OQ 的记录写入原始 QUAL，
  已存在的 OQ 不覆盖，并纳入 checkpoint/manifest 参数签名。
- ApplyBQSR 的静态质量分箱由 Host 按 GATK probability-space nearest/lower-bin
  规则生成固定映射；`--static-quantized-quals` 可重复指定，
  `--round-down-quantized` 只与静态分箱连用，`--allow-missing-read-group` 明确
  允许未知 RG 旁路 covariate model 并走该映射。
- 正数 `--quantize-quals N` 由报告 `Quantized` histogram 在 Host 侧重建 GATK
  `QualQuantizer` 的 N-level greedy adjacent merge；`0` 保持质量，`-1` 使用报告
  原有 map，三种模式在 ApplyBQSR encode 边界显式区分。

**Kokkos kernel**

- `count_bqsr_quality_kokkos` 已用 `TeamPolicy` 为每 team 建 188-bin local
  quality/mismatch counter，再按 league 顺序在 Host 做 deterministic integer
  merge；不再把所有 work-items 汇入一个全局 bin atomic。
- `apply_bqsr_quality_kokkos` 使用 `RangePolicy` 做独立质量变换；CIGAR、reference、
  known-sites 与 covariate-key 仍由 Host 投影，两个 execution policy 都写入 manifest。

BaseRecalibrator 与 ApplyBQSR 均由 runtime `ThreeStagePipeline` 驱动，队列按
byte capacity 背压，sink 保持 report/table/read order 的确定性；OpenMP/Serial
共享同一源代码和 benchmark driver。

ApplyBQSR 的 `--create-output-bam-index` 采用 GATK optional-boolean 语义：默认
生成 `.bai`/`.crai`，显式 `false`/`0` 时不发布索引；该选择进入 checkpoint
签名和 OutputManifest，避免恢复作业改变主输出之外的 sidecar 集合。

**验收**

- BaseRecalibrator 已输出 GATKReport v1.1 Arguments/Quantized/RecalTable0/1/2
  schema，并由 GatherBQSRReports 确定性合并；`--checkpoint` 在完整 batch
  边界原子写入 bins/covariates 和 input/options signature，`--resume-checkpoint`
  会校验 signature、跳过已完成 batch，并在 chr17 fixture 上逐字节重建 report
  与 sidecar。ApplyBQSR 现已在校验签名后复用已完成 BAM/CRAM 前缀，原子发布并
  重建 index；read qualities、GATK metrics、分片 merge 和完整 GATK recalibration
  model 仍需继续对照。
  GatherBQSRReports 会保留首个输入的 `Arguments`/`quantizing_levels`，按该
  levels 重建 Java `QualQuantizer`，并在合并前检查 covariate 与 indel-context
  table dimensions；非默认 quantization 及结构不一致的 fail-closed 行为已由
  独立 GATK 4.6.2.0 oracle 在 OpenMP/Serial 双 backend 验证。
  AnalyzeCovariates 已接入 Java/native report 输入、before/after/BQSR role
  选择、GATK intermediate CSV 生成和无外部 R 依赖的确定性 PDF；CSV 在 pinned
  chr17 report 上与 GATK 4.6.2.0 逐字节一致，accuracy 批量变换统一走
  `HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect`/`RangePolicy`
  并写入 execution-space/policy/records telemetry；多 report 合并前会校验共享 recalibration
  arguments，不一致时 fail-closed，PDF 绘图样式仍不宣称与 R 脚本逐字节一致。
  GATK 4.6.2.0 的 `AnalyzeCovariates` 短别名是小写 `-bqsr`；native 已补齐该
  alias，同时保留既有 `-BQSR` 兼容拼写。专用 pinned Java/native alias oracle
  固定 lower/long CSV 与 native lower/upper CSV 等价，并验证 `-bqsr=...` 的
  Barclay embedded-equals 形式 fail-closed；OpenMP/Serial 均纳入 CTest。
  `--compute-indel-bqsr-tables` 已覆盖 I/D default-quality、CIGAR event error、
  indel context 和四-base cycle cushion，并在 bundled chr17 fixture 上逐行
  对齐 Java GATK；同一 fixture 的 ApplyBQSR 输出再经 Java `PrintReads` 转成
  SAM 后，493/493 个稳定 read record（flag/CIGAR/sequence/逐碱基质量）与 GATK
  bit-identical，不再只比较 quality sum；长读/大样本与 cloud/restart corpus
  仍需扩展。计数必须整数一致，输出 read order/header 不变。
  `--use-original-qualities` 的 OQ reset/保留 tag 语义也已由 Java PrintReads 对
  493 条完整 record 逐行 oracle 验证；大样本与异常 OQ corpus 仍需扩展。
  BQSR 的 HTSlib decode 现在也接入 runtime `AdaptiveController`：首个 batch
  先按 SLURM/cgroup 保守 Host 预算限流，普通运行在安全 batch 边界观察
  `ReadBatch::bytes()` 后可继续减小下一批；checkpoint/resume 保持原有批次边界
  不变。`effective_batch_records`、`adaptive_batch_reductions` 和资源快照写入
  manifest，并由 1 MiB 调度分配回归验证。

### AnalyzeCovariates report-argument compatibility (2026-09-01)

AnalyzeCovariates before/after report 合并现按 GATK 4.6.2.0
`RecalibrationArgumentCollection.compareReportArguments()` 的 14 个固定字段比较；
`indels_context_size`、`covariate` 和 report 路径等 Java 明确忽略项不再误报，旧
report 缺失字段按 Java 默认值补齐。OpenMP/Serial 的 CSV oracle 逐字节一致，PDF
保持无 R 依赖的确定性语义输出，R 绘图布局仍显式 fallback。

### BQSR native-report producer/consumer gate

`BaseRecalibrator` 的 native writer 现会在完成 report 后重建 GATK v1.1
fixed-width table layout，并为 `Quantized` 写出 integer column types。此前
虽然四张数据表的值可与 GATK 对齐，但 TAB/短列头会让 HTSJDK 的
`GATKReportTable` 按列头位置错切，导致 Java `ApplyBQSR` 无法直接消费 native
report。新增 `verify_bqsr_report_roundtrip_gatk_oracle.py` 和
`fastgatk-bqsr-report-roundtrip-gatk-oracle`：pinned GATK 4.6.2.0 的
`ApplyBQSR` 直接读取 native report，再与 native `ApplyBQSR` 经 Java
`PrintReads` 比较 493 条 record exact；OpenMP/Serial 均通过。该 gate 证明
report serialization/direct replacement 边界，不宣称完整 BQSR 模型或任意旧版
GATKReport 变体。

### BQSR CRAM bounded gate

`ApplyBQSR` 的 CRAM writer 固定使用 HTSlib CRAM 3.0，并启用 reference-MD5
兼容模式：HTSlib 默认 CRAM 3.1 会被 pinned GATK 4.6.2.0/HTSJDK 拒绝，且
fixture 的 sequence-dictionary M5 与测试 reference 可能不同。新增
`verify_bqsr_cram_gatk_oracle.py` 使用真实 `dictFix.cram`/CRAI、known-sites，
比较 BaseRecalibrator 四张 GATKReport 表，并让 Java HTSJDK 重开 Java/native
ApplyBQSR 的 CRAM；OpenMP/Serial 均验证 493 条 SAM record exact。该设置只保证
CRAM 容器/reference/index 的直接替换边界，不宣称完整 BQSR 模型或 CRAM 字节相同。

### 4.7 Haplotype likelihood/genotyping/reference confidence/GVCF（P1/P2）

**C++ Host**

- 构造 allele/read likelihood matrix、genotype model、annotation、reference
  confidence block；管理 allele order、sample order、VCF/GVCF writer。
- GVCF block 是顺序 state machine，边界和合并不能由设备任意重排。
- 处理 low-depth/special allele/NaN/underflow 和 Strict/Fast numeric mode。

**Kokkos kernel**

- PairHMM 输出矩阵后的小型 likelihood transform、log-sum-exp、genotype
  likelihood、allele count 和 reference-confidence 数值 batch。
- CPU 用 SIMD/vector range；GPU 用小矩阵 TeamPolicy；归约树在 Strict 模式固定。
- 使用 Kokkos Kernels dense primitives，但不让 vendor BLAS 自行改变排序/NaN 语义。

**验收**

- GenotypeGVCFs/HC `-ERC GVCF` golden；VCF/GVCF record、block boundary、AD/DP/GQ、
  annotation、header/index 和 bio-identical truth。

当前 HC 已完成同一 filtered/projected observation stream 上的 REF-vs-any SNP
PL/GQ/AD、GQ-band block state、`MIN_DP` 和 GVCFBlock header；RCM 现在保留
observation→read 映射，按 GATK 默认 `indelSizeToEliminateInRefModel=10` 对
read/reference suffix 做 insertion/deletion plausibility 检查，且参数写入
Manifest。GATK 4.6.2.0 fixture 已覆盖区间、非重叠、PL shape 和 sentinel structural
oracle；终端 soft-clip 已按参考方向转换为左锚定 insertion candidate，并有正向/
reverse-strand regression；完整 unseen-haplotype posterior 和 bio-identical PL/GQ
仍是后续门禁。对于包含 concrete ALT 的候选位点，native 已按 GATK
`AlleleLikelihoods.updateNonRefAlleleLikelihoods` 的 qualified-concrete median
规则生成逐 read `<NON_REF>` likelihood，并在完整三角 PL 中参与 genotype
combination；graph path 的 bounded SW/CIGAR 反投影已覆盖 exact soft-clip 支持
的内部插入子集，graph deletion branch 也进入候选契约；CIGAR/soft-clip indel
在候选 map 前按 reference repeat 做确定性 left-normalization，并由 mixed
`M/=/X/I/D` regression 覆盖；短距离 compound `I/D` 组合也生成 anchored
ref/alt 候选；同一 bounded CIGAR span 内的 `I/D`+内部 soft-clip 现在也生成
anchored ref/alt，并有独立支持阈值；尚未覆盖的是 GATK 完整 assembly
haplotype 集合和跨 region 的复杂 compound soft-clip 证据。
终端 clip anchor read 在 RCM indel-informative 判断中按 GATK realignment-aware
的 max-indel-size suffix fail-closed，避免将 clipped evidence 当作 clean suffix；
insertion/deletion anchor 与 deletion element 也按 GATK pileup 排除条件
fail-closed，避免同一 read 同时被计为 indel event 和 clean-suffix evidence；
bundled chr17 的 terminal POS→END boundary 已由 Java oracle 对齐。

RCM Number=G 归约现在使用稳定 locus-segment staging：Host 只负责 stable sort 与
offset materialization，Kokkos `RangePolicy` 对每个 `(locus, genotype)` 仅遍历该
locus 的观测段；段内顺序与旧实现一致，保持 Strict PL/GQ bit semantics，同时避免
`locus_count × observation_count` 的无效扫描。polyploid kernel 的 benchmark
prepare/kernel p50、checksum 与 HC oracle 均纳入验证。三倍体
`17:69000-69100` 的 reference block（`PL=0,5,14,135`）以及候选位点的
GT/AD/DP/GQ/PL、MLEAC/MLEAF 已与 GATK 4.6.2.0 逐字段一致，OpenMP/Serial
均通过；这项 sentinel parity 不等于完整 assembly/全工具 bit-identical。

reference-confidence 的 Host indel-informativeness 判断也复用稳定的 observation→locus
CSR 索引：每个 locus 只遍历自己的 observation segment，并保留输入顺序与 40-read
上限，因此把原来的 O(loci×observations) 扫描降为 O(loci+observations)，不改变
PL/GQ 归约或 GVCF block 边界。该路径已通过 OpenMP/Serial HC 与 GenotypeGVCFs
oracle；reference-block prior、跨样本 posterior 和完整 Java 数值 parity 仍按
后续门禁推进。

2026-09-01 RCM PL range 增量：GATK `ReferenceConfidenceModel` 的
Number=G PL 不使用 PairHMM candidate envelope 的 999 截断；共享
`calculate_reference_confidence_genotypes_kokkos` 现在仅受 signed-int 输出范围约束，
Host indel RCM 分支也遵循相同边界。pinned
`fastgatk-hc-rcm-pl-range-gatk-oracle` 在 `--indel-size-to-eliminate-in-ref-model 1`
的高深度 BP_RESOLUTION 窗口中固定 PL=1080 与 Java 4.6.2.0 逐字段一致，OpenMP/Serial
均通过。该项只收敛 RCM PL 表示范围，完整 assembly、跨 region 和 GVCF 全文件
bit-identical 仍是后续边界。

PairHMM allele max-marginalization 的共享接口现在接受显式 concrete-allele
count，并以 row-major `[read-row][REF,ALT...]` 形式保留任意数量的 ALT；原有
两参数入口继续固定为 biallelic 兼容层。这样后续 HC/Mutect2 的完整 multi-ALT
read×allele posterior 可以直接复用同一 Kokkos `RangePolicy` kernel，而不会在
API 边界把证据悄悄截成 REF/ALT。HC 现已按稳定的
`(group,locus,source-read)` 行构造完整 concrete-ALT 索引（全局矩阵 stride 只用于
容纳不同位点的最大 ALT 数，未使用格保持 `-inf`），再把每个候选的具体 ALT
似然写回共享矩阵；多 ALT/多倍体 oracle 覆盖了该路径。Mutect2 现已把同一
`(contig,position,REF)` 下的 concrete sibling ALT 作为一个 locus-level VariantBlock
写出：VCF 使用逗号分隔 ALT，Number=A 的 TLOD/NLOD/AF/posterior/F1R2/R1F2
向量与 REF+ALT 的 AD/FAD，normal ALT 子集也按 row 的显式 alternate 索引映射，
避免把同位点 sibling 当作 REF。跨 AssemblyRegion 的完整 somatic/haplotype
posterior、bit-identical Java 数值和 joint filtering 仍按后续门禁推进。

HC 的 indexed `--stream-by-region` 现在使用 allocation-aware interval queue：当
当前 halo tile 的 HostBatch 超过 runtime 的保守内存预算时，Host 丢弃尚未输出的
tile 并把 core 稳定二分后重新查询；左右 child 按原坐标顺序消费，halo 每次重新
计算，输出只覆盖 core，因此不会牺牲跨边界证据或产生重复记录。安全预算额外保留
decoded-batch/HTS 元数据余量；单碱基 tile 若仍落在硬 SLURM/cgroup 配额内则允许
入队，超出硬配额才 fail-closed 为 `RESOURCE_EXHAUSTED`。`streamed_region_splits`
与峰值 staged bytes 写入 summary/manifest，作为内存门禁的可审计依据。

### 4.8 Mutect2 和 somatic sidecar（P1/P2）

**C++ Host**

- 复用 HC 的 AssemblyRegion/assembly/SW/PairHMM；Tumor/Normal sample mapping、
  GATK-shaped parameter validation、filter、VCF/stats/F1R2/somatic GVCF 保持
  Host，posterior 数值项交给共享 Kokkos kernel。
- sidecar 写入采用 `OutputManifest`，关闭 writer 前验证全部文件和 interval。
- FilterMutectCalls 的 native 硬过滤子集包括 TLOD、germline、污染 AF floor、
  minimum allele fraction、minimum reads per strand、unique ALT-read count、
  concrete ALT count、median MBQ/MMQ/MPOS thresholds 和 F1R2/R1F2 balance；
  过滤名、顺序、stats 计数和
OutputManifest 必须稳定。
`--tumor-sample` 已接通，AD-based hard filter 按 header 中的样本名选择，默认优先
`TUMOR`、否则使用首个样本，并把实际选择写入 stats/OutputManifest。
多 ALT 的 `TLOD`/`PSOMATIC`/`PGERMLINE`/`PARTIFACT` Number=A 向量现在按
concrete ALT 保持原始位置，站点级证据分别按 max/min 的最坏值归约；污染和
minimum-AF 使用全部 tumor ALT fractions，避免仅读取 ALT1。标量旧记录广播到
各 ALT，缺失向量元素不压缩位置，`AS_FilterStatus` 与站点过滤共用同一 allele
顺序，并由混合多 ALT regression 固定该边界。
Mutect2 的 candidate-emission boundary 与 HC 分离：`somatic_mode` 保留没有
PairHMM ALT owner 的 pileup candidate，并分别接受 GATK 的
`--initial-tumor-lod`（默认 2，active-region seed）与 `--tumor-lod-to-emit`
（默认 3，final emission）。当前 native activity engine 尚未独立消费前者，
因此 VCF/F1R2 写出前使用两者较严格值，并在 stats/manifest 记录原始与 effective
阈值；默认 call-set 与 GATK 默认发射边界保持不变。原始 `-L` selector 传入共享
caller，AssemblyRegion halo 候选不会泄漏到 core interval。恢复数量和候选摘要写入
 stats/manifest，作为可审计的 somatic pileup fallback。
Mutect2 的 `--force-active` 已按 pinned GATK 4.6.2.0 行为接入：它不改变
ActivityProfile 的坐标分段，只把原本 inactive 的 profile segment 标记为 active，
使 AssemblyRegion graph/PairHMM 消费完整 bounded window。该状态通过统一的
Kokkos activity-profile API 传递，stats/manifest 暴露 effective 值；
`fastgatk-mutect2-force-active-gatk-oracle` 比较 Java/native 的 22 行 IGV profile、
边界、active 状态和候选 site set，OpenMP/Serial 均通过。该切片只固定 activity/
writer 边界，完整 somatic posterior 仍为 release-specific fallback。
全局错配质量参数也保持 GATK 的边界：正数 Q 启用 per-read cap，负数以
`-Infinity` sentinel 关闭 cap，Q=0 在进入 assembly 前 fail-closed（GATK 的
`normalizeLikelihoods` 会拒绝该值），避免误用零宽 cap 改变 TLOD。disabled
sentinel 通过共享 Kokkos normalization API 保留 PairHMM 的合法 `-Infinity`
likelihood cells；`fastgatk-mutect2-mismapping-rate-boundary` 已在 OpenMP/Serial
与 pinned 4.6.2.0 Java 双侧通过，并包含 Java `[0.0,-Infinity,-3.0]` direct probe。
ClusteredEventsFilter 已按 GATK 语义消费 `ECNT`/`ECNTH`，支持
`--max-events-in-region`/`--max-events-in-haplotype`，并对缺失注释 fail-closed。

Mutect2 broad TLOD residual 已增加分层证据而不放宽 VCF 断言：新增
`fastgatk-pairhmm-results-oracle` C++/Kokkos replay target 与
`fastgatk-pairhmm-results-gatk-oracle` CTest。它读取 GATK 4.6.2.0
`--pair-hmm-results-file` 的真实可变长请求；706 行在 OpenMP/Serial 的
Kokkos bucketed PairHMM 均通过 `%e` 六位小数的固定 `5e-6` 输出边界（最大
`4.222168e-6`）。同一 fixture 的 native stats 还报告 724 request pairs、18
local haplotypes，GATK 为 706/18；request-set 差 18，且 native graph
haplotype considered/pruned 均为 0。该结果把已收敛的 recurrence 与尚未收敛
的 AssemblyRegion/ReadThreadingAssembler ownership 分开；完整 TLOD 最大差
`3.500717` 仍是 release-specific assembly 输入差异。下一步必须逐 region
对齐 haplotype/read request 集合，再复用 replay oracle；不以数值阈值掩盖差异。

Mutect2 writer 现复现 GATK 4.6.2.0 的调用阶段输出边界：site `QUAL` 保持为缺失值
`.`，不把 native somatic posterior 提前编码成最终质量；TLOD、PSOMATIC 等证据仍
写入 INFO，由后续 `FilterMutectCalls` 负责校准 QUAL。该行为已在 OpenMP/Serial
contract 与 pinned GATK oracle 中验证。

somatic depth aggregation 还会同步更新 `Result::candidates` 和已物化的
`Result::calls`。前者供 stats/候选 telemetry 使用，后者由 VCF writer 序列化；两者
若分叉会把同一个 locus 写成不一致的 FORMAT `DP/AD`。该同步已固定为回归门禁：
pinned GATK 4.6.2.0 chr17 fixture 的 site-set、GT、DP、AD 均为 6/6，TLOD
release-specific 数值差异（最大 3.500717）仍显式保留。

Mutect2 与 HaplotypeCaller 的 CLI Host 还显式启用 GATK 标准读过滤中的
`GoodCigarReadFilter` 和 `NonZeroReferenceLengthAlignmentReadFilter`；effective
mask 及 `read_filter_gatk_defaults` 会写入 telemetry/manifest。底层 shared Options
仍允许测试和自定义调用者使用宽松 mask，避免把 CLI policy 偷渡进 Kokkos API。
流式 region 的空 tile 不含任何 alignment 时，read-filter kernel 对空 batch 直接返回
空结果，不要求不存在的 CIGAR metadata，从而保持跨 contig tile 的 bounded queue
可用性。

Mutect2 的标准 mask 现在同时启用 `NonChimericOriginalAlignmentReadFilter`：HTSlib
将 OA/XM Z tags 以 flat offsets/payload 保存在 `ReadBatch`，跨 shard/tile merge 和
filter-select 均保持 record 对齐；Kokkos 在设备端提取 OA 首个逗号前的 contig 与 XM
比较，任一 tag 缺失时按 GATK 语义通过。API smoke 覆盖匹配、不匹配和缺失三态，
HC 仅在显式 `--read-filter` 时启用该类。

同时按 GATK `StandardMutect2ReadFilters` 启用 `ReadLengthReadFilter(30, ∞)`；长度
上下界沿用同一 Kokkos `ReadFilterOptions`，合成多等位/soft-clip 回归使用 ≥30bp
序列，避免测试数据绕过真实默认过滤。

该标准 mask 还显式启用 `MappingQualityAvailableReadFilter` 与
`MappingQualityNotZeroReadFilter`：设备端统一拒绝 MAPQ=255（HTSJDK unavailable
sentinel）和 MAPQ=0，且两项的 effective 状态会写入 telemetry；显式 `-RF/-DF`
仍可覆盖默认策略。

`FilterMutectCalls` 还接通了标准 `LearnReadOrientationModel` hand-off：
`--orientation-bias-artifact-priors`/`--ob-priors` 读取 `.orientation_priors`
tar.gz，FASTA Host 解析 3-mer 和 GATK FORMAT `F1R2/F2R1`，随后用 Kokkos
12-state beta-binomial kernel 计算 orientation artifact posterior，输出
GATK-shaped `ROQ` 与 `orientation` FILTER。未显式给阈值时使用 0.5 的确定性
adapter threshold；现在按 GATK ReadOrientationFilter 对全部非 normal tumor
genotype 以全 ALT AD depth 做 weighted-median 聚合，缺失 sample prior 仍以
posterior=0 参与权重；显式 --tumor-sample 可覆盖为单样本。GATK 的
multi-pass threshold learning 和完整 joint filtering 仍保留 fallback，并已接入
多样本 contract、dispatcher 和文件边界 benchmark。兼容的
CONSTANT/FALSE_DISCOVERY_RATE/OPTIMAL_F_SCORE threshold strategy 对 native
orientation posterior 采用两遍流式学习，完整跨过滤器 joint learning 仍保留
fallback；等长 MNP 现在按 GATK 逐碱基取最大 orientation posterior，indel 仍
fail-closed。orientation beta-binomial 现在将 AD allele-depth 与 F1R2/F2R1
orientation trials 分开消费，并把两者不一致的记录数写入 stats/OutputManifest
telemetry，避免由 AD 隐式重建 strand evidence。2026-08-31 又补齐了标准
Number=R F1R2/F2R1 的 REF+ALT 解析：orientation posterior 的 trial depth
严格使用 `f1r2[REF]+f2r1[REF]+f1r2[ALT]+f2r1[ALT]`，FORMAT/AD 只保留为
GATK weighted-median 的 sample weight；Number=A/scalar 旧 native 记录仍走明确
compatibility fallback。双 backend verifier 新增“相同 orientation vectors、不同
AD”回归，ROQ 必须相同，证明 AD 不再泄漏进 beta-binomial posterior。

**Kokkos kernel**

- 只处理复用的 assembly/likelihood/genotype primitives；somatic filter 的分支和
  证据解释先 Host。
- F1R2/likelihood histogram 采用 per-team local reduction，结果按 stable allele id
  恢复顺序。
- Tumor/normal PairHMM 的 per-read REF/ALT likelihood 进入独立
  `calculate_somatic_likelihood_kokkos` mixture-grid kernel，并由
  `calculate_somatic_posterior_kokkos` 结合 somatic/germline/artifact priors、
  F1R2/R1F2 orientation evidence 和 contamination correction；Host 负责
  GATK-shaped 参数校验与 FILTER/sidecar 编码，kernel 输出 TLOD、best/corrected
  allele fraction、三类 posterior 和 orientation posterior。
- FilterMutectCalls 的标量 evidence、slippage、normal-artifact、orientation 和
  ErrorProbabilities 边界也统一经过 `HostBatch -> KernelPlan.prepare -> Kokkos
  Views -> RangePolicy -> execute -> collect`；manifest 记录 execution space、
  kernel batches/observations 及 prepare/execute 时间，避免 scalar bridge 绕过
 统一生命周期。当前多 ALT record 已合并为一次 `RangePolicy` launch，单 ALT
保留低开销 scalar fast path；后续可再跨 record 做持久化 shard fusion，但不改变
当前 Host 的 GATK filter 顺序。
- FilterMutectCalls 的 GATK 通用 `--sites-only-vcf-output true|false` 也已接入：
  过滤阶段继续读取完整 FORMAT/sample 证据，最终 writer 通过 HTSlib sample mask
  输出 8 列 site-only VCF，并保留 INFO/FILTER、索引、stats 与 manifest；双 backend
  contract 及 pinned Java 4.6.2.0 shape oracle 已验证。
- FilterMutectCalls 现在还实现 GATK `GermlineFilter.germlineProbability` 的一个可审计
  后验边界：Host 提取 NLOD/POPAF/肿瘤 AD，Kokkos 计算 germline-het/hom-alt 与
  somatic likelihood 的归一化后验，并按 `QualityUtils.errorProbToQual` 的 [1,93]
  规则写出 Number=A `GERMQ`。POPAF 缺失时 fail-closed，不把可选 `PGERMLINE`
  当作 GERMQ；四条 pinned GATK 4.6.2.0 fixture 在 OpenMP/Serial 均通过。完整
  release-calibrated joint posterior/filter learning 仍是显式 fallback。
- FilterMutectCalls 的 `-L/--intervals`/`--region` 采用统一 HTSlib interval selector，
  支持 `--interval-set-rule UNION|INTERSECTION`，并在全部学习 pass 与最终过滤中保持
  相同子集；stats/manifest 记录规则、输入统计与跳过计数，双 backend contract 已覆盖。
- contamination filter 进一步复现 GATK 的 AD/POPAF 单污染者/多污染者 likelihood、
  SomaticClusteringModel posterior 与 ALT-depth weighted median；每个 tumour/ALT
  观测进入同一 Kokkos `RangePolicy` batch，`--max-contamination-probability` 决定
  AS/site 状态，`CONTQ` 写出 phred contamination quality。缺少 POPAF 或 AD 时保留
  legacy AF-floor，并在 manifest/stats 中区分两条路径。

**验收**

- Mutect2 integration、FilterMutectCalls、stats/F1R2/VCF/GVCF 配套文件；同一
  tumor-normal 输入严格/生物学结果比较。

Mutect2 的 `--f1r2-tar-gz` 现在按输出后缀提供两条明确的文件边界：`.tar.gz`
写出可直接交给 `LearnReadOrientationModel` 的标准 POSIX tar.gz，包含每个样本的
`.ref_histogram`、`.alt_histogram` 和 `.alt_table` 三成员；Host 从已过滤的
ReadBatch/参考序列重建 reference-locus histogram，并保留 Kokkos 产生的 ALT/
orientation 证据。显式 `.tsv` 路径保留旧版 fastgatk sidecar，manifest/stats
记录 `f1r2_standard_tar`、`f1r2_legacy_tsv` 和格式类型，避免把近似格式伪装成
GATK 输入。
`benchmark_mutect2.py` 对该标准 sidecar 路径报告重复运行的 wall/p50、VCF/sidecar
大小及 somatic likelihood/posterior execution-space telemetry。

Mutect2 的 tumor/normal HTSlib staging 现在也统一使用 `ReadBatch::bytes()`；
`AdaptiveController` 在每个安全 batch 边界调用 `HtsReader::set_batch_records()`，
stats/OutputManifest 记录 requested/effective batch 与 reduction 次数。当前
默认 assembly 仍保留稳定的 aggregate HostBatch；显式 `--stream-by-region` 才启用
indexed core/halo bounded path，因此默认路径仍不宣称 cloud spill。

本轮新增 Mutect2 `--stream-by-region N`：对有 `.bai/.crai` 的输入按 core/halo
区间重新打开 HTSlib reader，core 坐标作为 caller 输出域，tumor/normal somatic
likelihood、posterior、VCF 编码和 F1R2 sidecar 通过统一
`ThreeStagePipeline<Decoded,Computed,Encoded>` 进入唯一 sink。tile 超过
runtime byte budget 时按 midpoint 递归拆分；sink 只保留 F1R2 histogram/alt-row
归约和计数 telemetry，不保留完整 `Result`。chr17 fixture 的 aggregate 与
`--stream-by-region 500` 解压 VCF 逐字节一致，标准 F1R2 tar 三个成员内容也一致；
OpenMP/Serial contract 均通过。该路径要求 indexed BAM/CRAM，跨 tile 的 decoded-read
计数包含 halo 重叠，manifest 明确记录 `streamed_regions`/峰值 bytes；完整
Mutect2 release-specific somatic model、cloud spill 和全局去重统计仍不宣称完成。
本轮将该路径的 `calling::Options::intervals` 与 `locus_intervals` 明确分离：前者
只负责 core 输出过滤，后者把 halo 传入 activity/graph/PairHMM，避免边界证据丢失
或重复写出。500bp fixture 仍保持 aggregate 字节一致；100bp 压力 contract 只检查
有序、core-only、无重复 locus，因为 release-specific somatic posterior 在独立 tile
上并不具备可结合性，不能把该压力路径误报为 bit-identical。

Mutect2 的 `--tumor-sample`/`--normal-sample` 现在在同一 Host staging 边界生效：
HTSlib 按 `@RG ID→SM` 过滤记录后才进入 assembly、PairHMM、fragment grouping 和
F1R2 计数；选中的样本名以及 tumor/normal 过滤后 read 数写入 stats/manifest。
未知样本 fail-closed，避免多样本 BAM 被静默合并；未显式指定 tumor 时选择
header 中的第一个 sample，而没有任何 `@RG SM` 的 tumor 输入在聚合及
`--stream-by-region` metadata probe 阶段即 fail-closed（与 GATK 的
`samples cannot be empty` 一致），不再以合成 `TUMOR` 名称接收全部 reads。该
语义由双样本 synthetic regression、无样本 native/Java oracle 和双 backend
contract 覆盖。

标准 GATK 的单输入 tumor/normal 调用也已落到同一边界：只给一个多样本 `-I` 并
指定 `--tumor-sample` 与 `--normal-sample`/`-normal` 时，native 会复用同一
HTSlib reader 分别构造 tumor 与 normal ReadBatch；`--normal-input` 仍用于分离
BAM/CRAM 的兼容调用。双样本 synthetic contract 固定过滤后 read 数和
`TUMOR/NORMAL` VCF 列。aggregate 路径还支持重复 `-I` shard，先校验全局
`@SQ`/assembly metadata，再按选定 sample 过滤并做稳定 coordinate merge（argv 顺序
作为 tie-break）；indexed `--stream-by-region` 现在逐 shard 校验 index/dictionary，
并在每个 bounded core/halo tile 内复用同一 merge，避免静默漏片。
完整 joint-sample posterior 仍是后续实现项。

该 joint posterior 的一个实际数据边界已收口：tumor 与 normal 的
`AlleleLikelihoods` 保留独立 evidence width，normal coverage 不再被错误地要求
等于 tumor fragment 数。Kokkos posterior 从 normal candidate-major matrix 推导自身
行宽，并保持 normal penalty 每个 evidence unit 只消费一次；新增
`fastgatk-somatic-posterior-normal-count-gatk-oracle` 固定 GATK 4.6.2.0
`SomaticGenotypingEngine` 语义的 4-tumor/2-normal fixture，OpenMP/Serial 均通过。
其余 release-specific joint calibration 仍显式 fallback。

Mutect2 的重叠双端质量模式也与 Java 调用点保持分流：shared caller 在
`somatic_mode` 下仅封顶 concordant bases，保留冲突 mate 的原始质量（等价于
`setConflictingToZero=false`）；HC 仍按 GATK conflict-to-zero 语义执行。该边界
修正了重复 `-I` 合并后 broad oracle 的额外 indel，OpenMP/Serial 均恢复 6/6
site-set，bit-identical 仍未宣称。

Mutect2 的重复 `-L/--intervals` 现在同步接受 GATK `-isr/--interval-set-rule`
（`UNION` 默认、`INTERSECTION` 可选）。aggregate 与 indexed
`--stream-by-region` 共用同一 HTSlib Host interval-set 运算，避免 scatter tile
与一次性运行使用不同输入子集；实际规则写入 stats/OutputManifest，并由双后端
contract 覆盖。该边界只承诺 selector/IO 一致性，release-specific somatic posterior
仍保持显式 prototype/fallback。

Mutect2 还实现了 GATK 通用 `--sites-only-vcf-output true|false`：调用和
F1R2/统计计算不变，Host writer 在 site-only 模式只写 8 列（不写 FORMAT 与样本列），
同时保留 BGZF/Tabix、stats 和 manifest；双 backend synthetic contract 已覆盖开关、
INFO 保留和输出列数。

Mutect2 的 `-normal`/`--normal-sample` 现在支持重复 selector。Host 按 `@RG ID→SM`
逐样本构造 normal ReadBatch，合并后的 normal evidence 用于统一 somatic posterior，
但 writer 按 argv 顺序为每个 normal 保留独立的 GT/AD/AF/DP genotype 列；
`normal_samples` 数组与聚合 read 数同步写入 stats/OutputManifest，单 normal 旧调用及
`--normal-input` 分离文件形态保持兼容。双 normal synthetic contract 已覆盖
normal_reads=12、三列 VCF header 以及 OpenMP/Serial 两个 backend；完整 joint
posterior 与 release-specific bit-identical 仍显式 fallback。

Mutect2 输出 writer 现在支持 GATK `--create-output-variant-index[=true|false]`
optional boolean，默认开启。对压缩 VCF，`false` 会抑制 sibling `.tbi`，但不改变
caller/posterior、VCF 内容或 F1R2/stats，并将 `create_output_variant_index` 与实际
`vcf_index` 状态写入 stats/OutputManifest；OpenMP/Serial verifier、targeted CTest 和
dispatcher registry contract 已覆盖。未压缩 VCF 在开关开启时生成
HTSJDK/Tribble LinearIndex v3 sibling `.idx`，并由 GATK SelectVariants interval query
实读验证；关闭时不生成 `.idx`。

普通 VCF 索引已进一步抽为共享 `fastgatk-core` Tribble writer：
HaplotypeCaller、FilterMutectCalls、SelectVariants、VariantFiltration、
LeftAlignAndTrimVariants、ApplyVQSR 均统一输出 `.vcf.gz/.tbi` 或普通
`.vcf/.idx`，并由 OpenMP/Serial verifier 的 bundled GATK interval query 验证。
GenotypeGVCFs、CombineGVCFs 和 ReblockGVCF 现也复用该 API；三者的 OpenMP/Serial
普通输出已由 GATK `SelectVariants -L` 实读验证，关闭索引开关时两类 sidecar 均不生成。
GatherVcfs 同样复用该 API，普通 VCF 的 `.idx` 已通过 OpenMP/Serial GatherVcfs
oracle 和 GATK 区间查询验证。

2026-08-28 的 pinned GATK 4.6.2.0 broad oracle 进一步固定了当前边界：同一
chr17 fixture 上 VCF/Tabix、单样本列、F1R2 三成员 tar、GT/DP/AD 结构均可比较，
当前 bounded somatic assembly 已将 native call-set 提升到与 Java 相同的 6 个
site（site recall/precision=1.0），并通过 GATK 形状的 somatic `0|1` GT 与局部
读段相位恢复使 broad fixture 的 GT 达到 6/6 exact；Host 现在保存过滤后 read 的
坐标元数据，在 fragment grouping 后按 GATK informative-overlap margin 保留候选
证据，shared-site TLOD 最大绝对差由约 29.34 降至 5.908758（平均 1.145318），
后验仍未全部一致。该差异来自尚未完成的 GATK release-specific
assembly、read-orientation 和 somatic-prior calibration；Mutect2 继续标记为
`prototype/fallback`，不能把当前 Kokkos Dirichlet evidence kernel 误报为
bit-identical。低支持候选只在固定 haplotype/link budget 内采样，避免为提高召回
引入无界 `2^N` PairHMM 工作。

2026-08-29 又将 fragment 边界前移到 haplotype likelihood：Host 在通过同一
CIGAR/graph ownership gate 后，按 read name 聚合每个 fragment 内的 haplotype
log10 likelihood，再对 REF/ALT 做 max-marginalization，结果以可选的
candidate×fragment 矩阵交给 somatic Kokkos kernel。这样不再把“每条 read 先
边缘化、再求和”误当成 GATK `groupEvidence(...).marginalize(...)`；同一位点
的两个或更多 concrete ALT 还会把携带兄弟 ALT 的 haplotype 标为 unknown，避免
sibling ALT 被重复当作 REF。稀疏候选仍保留重叠 fragment 上的 pileup fallback，
并由 synthetic multi-ALT、pinned broad oracle 和 OpenMP/Serial contract 覆盖。

多等位 somatic evidence 现在由共享
`calculate_somatic_multiallelic_likelihood_kokkos` 统一计算：输入是
allele-major REF/ALT（可扩展末尾 `<NON_REF>`）log10 likelihood matrix，设备端按
GATK `SomaticLikelihoodsEngine` 的 flat-Dirichlet posterior 与 entropy 规则计算
all-alleles evidence，再逐 ALT 重算 without-ALT evidence。Host 只负责按 locus
合并 fragment/read-name rows 和把结果映射回候选；singleton 继续使用 biallelic
路径并有 parity smoke。该接口已在 OpenMP/Serial backend 通过 API smoke，Mutect2
synthetic multi-ALT contract 和全量回归；完整 release-specific assembly 与
prior calibration 仍是后续 oracle 门禁。

Strict 模式下，biallelic 与 multiallelic responsibility 归一化复现 GATK
`NaturalLogUtils.logSumExp` 的首个最大项/`1.0` 累加顺序及首最大值 tie-break，
并用 pinned Java `logEvidence(all)-logEvidence(without-ALT)` 数值固定 API
回归；这只保证该 Kokkos 数值 kernel 的运算顺序，不替代完整 Mutect2 assembly
和 joint posterior parity。

`LearnReadOrientationModel` 已加入 native hand-off：读取 Mutect2 生成的
F1R2/R1F2 TSV sidecar，或标准 `CollectF1R2Counts` tar.gz 的
`.ref_histogram/.alt_histogram/.alt_table` 三类成员；TSV 兼容路径接受 plain 或
`.tsv.gz`（zlib Host 解码）并使用 Kokkos
atomic reduction，标准 tar 路径实现 Java 的 12-state Beta-Binomial EM、
canonical/reverse-complement context 合并和 prior pseudocounts。标准 EM 的
每个 context observation batch 以及 legacy TSV aggregate 都通过统一
`HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect` 生命周期和
`RangePolicy` 执行，Host 仅按 Java 输入顺序归约，manifest 同时记录 execution
space、policy、observation/batch 数量及 prepare/execute 时间。输出为 64 个 context
的 GATK-shaped `.orientation_priors` 表，
以 POSIX tar.gz、样本 metadata 和 OutputManifest 交付。64 个 context 的输出行序
复现 Java `HashMap` 迭代顺序；无数据 context 行同时复现 Java `Double.toString`
格式，因此这部分可逐字节比较。四个真实标准输入/输出路径与 GATK 4.6.2.0
oracle 的 64 contexts 数值最大差为 4.5e-14（当前 OpenMP/Serial 后端），差异
来自 Commons-Math/Kokkos special-function 路径；legacy TSV 仍是显式兼容近似，
压缩包成员路径/gzip 字节、完整 bit-format identity 和下游完整
FilterMutectCalls 语义仍待完成。

`CollectF1R2Counts` 已补齐标准生产端：HTSlib Host 读取 BAM/CRAM，复用 Mutect2
默认 read-filter 表面并做 CIGAR-aware pileup 投影，Kokkos atomic kernel 聚合每个
locus 的 A/C/G/T 与 F1R2/F2R1 计数；该 kernel 统一走
`HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect`/`RangePolicy`
生命周期，Host 完成 median-MQ、indel、context 和
histogram/table 编码，输出标准三成员 `.ref_histogram/.alt_histogram/.alt_table`
POSIX tar.gz。chr17 fixture 与 GATK 4.6.2.0 的 64 个 ref histogram label、384 个
alt histogram label 和 3 行 alt table 逐字段计数一致；`verify_collect_f1r2_counts.py`
和 `benchmark_collect_f1r2_counts.py` 已接入。OA/chimeric edge cases、cloud 输入和
report byte identity 仍保持显式 fallback。指定 `-L` 且存在 BAM/CRAM index 时，
HTSlib 使用 indexed region iterator；无 index 时保留 sequential fallback，并把
`indexed_inputs`、`sequential_inputs` 与 `iterator_intervals` 写入 telemetry，保证
I/O 策略可观测且不改变计数 contract。多输入聚合 key 使用 `(sample, contig-name,
position)`，不依赖每个 BAM/CRAM header 的局部 TID；异序 contig header 的混合输入
已有双输入 regression。`-XL/--exclude-intervals` 现从 locus stream 剔除 site 而不丢弃
跨越该 span 的 read；标准 tar 三成员内容由 pinned GATK 4.6.2.0 oracle 对照。无 `-L`
的 whole-reference exclusion 暂显式 fail-closed。manifest/JSON telemetry 同时记录 execution space、kernel
batch/observation 数及 prepare/execute 时间。

### 4.9 GenomicsDBImport / GenotypeGVCFs adapter（P1/P2）

**C++ Host**

- 不从零重写 GenomicsDB/TileDB；复用 backend，负责 batch-size、reader 数、VCF
  buffer、fragment/consolidation、tmp/scratch/fd 预算，并写出确定性的
  `fastgatk-inputs.tsv` 输入索引供 native `gendb://` 读取。
- `--fastgatk-native-workspace`（或 `FASTGATK_GENOMICSDB_BACKEND=native`）提供无
  Java 的 sparse-index workspace（另写 `fastgatk-workspace.json`），仅承诺 native
  GenotypeGVCFs 可读；同时将本地输入物化到 `native-inputs/` 并写出
  `fastgatk-native-inputs.tsv`，因此 workspace 不依赖原始输入路径；manifest 明确
  标记非 TileDB、非 bit-identical，默认仍调用配置的 GATK GenomicsDB/TileDB backend。
- native workspace 保留命令行 `-V` 与 sample-map 的显式 callset/sample 顺序，只做
  稳定去重，不按文件名排序；因此重开后的 GenotypeGVCFs header/FORMAT sample 列
  仍与 GATK import 顺序一致。`verify_genomicsdb_import.py` 以逆序输入、重复输入和
  chr2-only shard 的 chr1 query 覆盖这一边界（S2 保留在 header、body 为 missing）。
  该保证只适用于当前 native sparse index；外部 TileDB/GenomicsDB workspace 仍由
  配置的 GATK backend/legacy bridge 负责，native 模式不宣称 TileDB 存储兼容。
- workspace 写入现在先取得原子 sibling-directory lock
  (`<workspace>.fastgatk.lock`)，覆盖外部 backend fork、native materialization 和
  manifest publish；并发的 Nextflow/SLURM retry 以 `RESOURCE_EXHAUSTED` fail-closed，
  正常/异常路径由 RAII 释放，lock owner 记录 pid/thread，manifest telemetry 记录
  lock path，避免把半写 workspace 当成可复用输入。
- native materialization 现在按输入文件原子更新
  `fastgatk-native-workspace.checkpoint`，记录 canonical path、size/mtime 与完成位；
  `--resume-native-workspace` 在重试时严格校验输入 fingerprint、复用已完成的
  `native-inputs/` copy，并只在 record index/metadata 成功发布后删除 checkpoint。
  这为 Nextflow `retry`/SLURM preemption 提供可审计的文件级 restart 边界，改变的
  输入或缺失的 checkpoint 均 fail-closed。
- native `gendb://` 读取若发现 sidecar，会 fail-closed 校验 schema/backend/input_count
  及本地输入文件；native 模式另写 `fastgatk-record-index.tsv`（每个输入的
  contig/start/end/record count），`GenotypeGVCFs -L` 会在保留所有 sample header
  的前提下跳过无关 shard body；无 sidecar 的旧 `fastgatk-inputs.tsv` workspace
  继续兼容。查询按 contig 名称匹配而非跨 header 复用 numeric RID，因此只含
  chr2 的 shard 在 chr1 查询中可安全跳过而 S2 仍保留在输出 header。首个 shard
  不含 `-L` contig 时仍保留 canonical selector name，CSI/TBI traversal 在每个
  shard 按名称映射 index tid（CSI 不复用 header-local RID）。该索引是 portable
  sparse index，不是 TileDB/GenomicsDB 存储替代；跨字典 shard 的 contig 会合并到
  全局 header，记录 RID 按名称重映射，FORMAT 数组延迟到全局 sample merge 后写出。
- 对真实 GATK TileDB/GenomicsDB workspace，若存在匹配的 GATK
  `libtiledbgenomicsdb.so`/`libjvm.so`，`fastgatk-genomicsdb-export` 通过旧 C++ ABI
  的独立进程桥接读取 `callset.json`、`vidmap.json`、`vcfheader.vcf` 和数组，按数组
  输出临时 VCF shards，再交给同一 native GenotypeGVCFs/Kokkos 路径；桥接器不进入
  主程序 ABI，缺失或版本不匹配时保持 fail-closed。`verify_genomicsdb_bridge.py`
  使用真实 GATK workspace 比较 GT/AD/PL，明确不宣称下游 QUAL/annotation
  bit-identical。
- `verify_genomicsdb_import_gatk_oracle.py` 进一步把 adapter 本身接到 pinned
  GATK 4.6.2.0：通过 wrapper 保留完整 `GenomicsDBImport` argv，验证真实
  `callset.json`/`vidmap.json`/`vcfheader.vcf`/array 发布、输入索引与 manifest，
  再用同版本 `GenotypeGVCFs -V gendb://...` 重开 workspace。该 CTest 与
  `fastgatk-genomicsdb-bridge-gatk-oracle` 在 OpenMP/Serial 各自通过；它只证明
  本地 GATK backend 的 workspace publication/reopen 边界，不提升 native sparse
  workspace 的 TileDB 兼容性声明。
- locus 迭代、sample/allele order、restart/resume 和 workspace 可读性由 Host。

**Kokkos kernel**

- 仅对已取出的 locus block 做 genotype likelihood/annotation 数值计算；不把
  sparse array 的随机访问伪装成 Kokkos bulk kernel。

**验收**

- GenomicsDB workspace 可被原 GATK 读取；数千 sample 的 fd/memory/tmp 压力；
  GenotypeGVCFs VCF 与 GATK 对照；同一 locus 的跨 shard concrete ALT union
  必须按 allele name 重映射 GT/AD/PL，`<NON_REF>` 不得被当作真实 ALT 或
  REF-only block 混入变体记录；reference-only `<NON_REF>` block 在存在
  concrete ALT 的 cohort locus 必须按 REF/<NON_REF> PL 投影为 hom-ref 样本；
  diploid GT/GQ 应由 remapped PL 决定，并输出 AC/AN/AF；`--heterozygosity` /
  `--indel-heterozygosity` 使用 GATK `GenotypePriorCalculator.assumingHW`
  的 SNP/INDEL/OTHER 公式；`--use-posteriors-to-calculate-qual` (`--gp-qual`)
  遵循 GATK `hasPosteriors` 条件：只有显式 posterior assignment 生成 GP 时才
  通过统一 Kokkos posterior kernel 计算可复现的独立样本 hom-ref posterior；
  PL/GQ-only gVCF 保持默认 QUAL，并写入 OutputManifest telemetry。该路径不改写 PL/GT/GQ，完整 GATK
  cohort AF likelihood 已按 GATK `AlleleFrequencyCalculator` 的 EM/Dirichlet
  pseudocount/convergence 规则移植到 Kokkos，并写出 MLEAC/MLEAF 与 cohort
  QUAL；共享 genotype-prior 的 cross-sample reference-confidence Kokkos
  posterior 已写出 `RCQ`/`RCP` 与 telemetry；release-specific calibration、完整
  joint posterior 和 GenomicsDB workspace 仍必须显式 fallback。

  `--genotype-assignment-method` 现在接受 GATK 的 `PREFER_PLS`（默认）、
  `USE_PLS_TO_ASSIGN`、`SET_TO_NO_CALL`、`SET_TO_NO_CALL_NO_ANNOTATIONS`、
  `BEST_MATCH_TO_ORIGINAL`、`DO_NOT_ASSIGN_GENOTYPES`，以及显式
  的 `USE_POSTERIOR_PROBABILITIES`。后者消费 GATK `GenotypePriorCalculator.assumingHW`
  的 per-site genotype prior，按同一 VCF genotype rank 在 Kokkos 中计算 posterior
  GT/GQ，并记录独立 telemetry；
  `SET_TO_NO_CALL_NO_ANNOTATIONS` 已按 Java `GenotypeBuilder` 语义清除所有
  genotype-level FORMAT 注释（保留 no-call GT；包括 GP/PG/PP、FT/MIN_DP 和
  caller-specific 字段），并记录清除字段计数；
  `BEST_MATCH_TO_ORIGINAL` 已保留原始 GT，并按 GATK GQ=0/PL[0]=0 规则转 no-call；
  `USE_POSTERIORS_ANNOTATION` 现在保留 FORMAT/PP、按 Number=G 重映射并使用
  Kokkos phred-posterior GT/GQ assignment。
  `USE_POSTERIOR_PROBABILITIES` 还写出 FORMAT/GP 与 FORMAT/PG 的 phred-scaled
  posterior/prior vectors，数值来自同一 Kokkos posterior kernel。

  reference-confidence block 物化现在在 Host 侧一次构建 locus→metrics 坐标索引，
  并将候选 REF span 按 contig 合并后执行二分 membership 查询；全参考 GVCF 不再对
  每个参考碱基线性扫描 observed loci/candidates。该优化不触碰 Kokkos
  likelihood/PL 数值，telemetry 记录 `reference_block_lookup_indexed=true`，并由
  `verify_indel.py` 的 interval/full-reference GVCF 回归门禁。

  对 `--include-non-variant-sites` 的纯 `REF/<NON_REF>` block，Host 在压缩
  writer-facing REF-only PL 之前保留源宽度 PL，并在跨 shard/sample 合并后按全局
  sample index 重建 source-width 矩阵，再调用共享 Kokkos cross-sample
  reference-confidence posterior，输出 `RCQ/RCP`。因此 dense block 不会把
  `PL=[0]` 当作 cohort certainty；零深度/缺失样本仍保持 no-call、无注释的
  fail-closed 语义。双样本 11-base regression 已在 OpenMP/Serial contract 与
  GATK GenotypeGVCFs verification 中覆盖，manifest 记录 22 个 sample-locus
  posterior 输入。

  cohort AF 还覆盖 GATK 对 PL-less reblocked genotype 的兼容边界：diploid
  hom-ref 且带 GQ 的样本由共享 Kokkos API 合成 0/GQ/10×GQ approximate PL，
  其它无 PL 或非 hom-ref 样本保持不可用于 AF 计算；approximation 数量和
  backend timing 写入 OutputManifest，并有 GQ-only regression。

  legacy AF-calculator 也有独立 oracle：`--use-new-qual-calculator false` 不
  启动新的 cohort EM，但 GATK 仍会将 HaplotypeCaller 携带的 Number=A
  `MLEAC/MLEAF` 从 concrete ALT+`<NON_REF>` 投影到最终 ALT。native 在 ALT
  union 和移除 symbolic allele 时复用 Kokkos allele-field remap；OpenMP/Serial
  各 3 条真实记录均与 Java 整行及 semantic header exact，且 manifest 明确
  `cohort_af_kernel_calls=0`。

GenotypeGVCFs 的跨 shard ALT union 与 `<NON_REF>` 投影已统一调用
`fastgatk-kernels` 的 Kokkos PL-remap API；该 API 支持缺失目标 allele、反向 ALT
顺序和固定任意 ploidy，相关多样本/triploid regression 与
`benchmark_genotype_gvcf.py` 已接入 `verify_all`；FORMAT/INFO AD Number=R 复用
companion allele-field Kokkos API 并记录独立 telemetry。输入带 CSI/TBI 且指定
`-L/--intervals` 时，Host 通过 `bcf_itr_queryi`/`tbx_itr_next` 做 indexed
traversal（TBI 先按 contig name 映射），重叠 query 先归一化；缺少索引时保留
sequential fallback，`indexed_inputs`、`indexed_interval_queries` 和
`interval_skipped` 写入 manifest。该边界只优化 decode/I/O，不改变下游 Kokkos
joint genotype、ALT union 或输出排序。
`benchmark_genotype_gvcf.py` 还用 HC 生成的真实 BGZF/TBI gVCF 测量 indexed interval
路径；OpenMP 当前 fixture 的 indexed p50 约 0.115 s，作为同机 I/O 回归基线，不作为
全基因组性能承诺。
GenotypeGVCFs 的 Host 初始化还会优先消费 `OMP_NUM_THREADS`，否则消费
`SLURM_CPUS_PER_TASK`，并通过 `ResourceSnapshot::effective_threads()` 做最终上限，
避免小 interval/低负载任务在未配置 OMP 时启动超出 SLURM allocation 的 OpenMP team；
manifest telemetry 记录 `kokkos_default_concurrency`。
重复 `-L/--intervals` 现在按 GATK `--interval-set-rule UNION|INTERSECTION`
解析：每个 selector 先在 Host 侧归一化为半开区间集合，再对 selector group 做
确定性并集或交集；indexed/sequential 和 `--stream-by-locus` 三条路径共享该结果，
并将规则写入 stdout、OutputManifest compatibility/telemetry。真实 HC 生成的 gVCF
交集区间已加入 `verify_genotype_gvcf.py`，OpenMP/Serial 均通过。
`--create-output-variant-index` 同样采用严格的 GATK optional-boolean 解析（裸参数
及 `true/false`、`=true/=false`），默认生成 `.tbi`；显式 `false/0` 不生成索引，
并在 `compatibility.vcf_index` 与摘要中记录。非法字面量直接 fail-closed。
本轮再把 materialized locus 的最终 annotation/serialization 阶段接入统一
`ThreeStagePipeline`：Host decode 只移动已排序的 `Record`，compute worker 顺序调用
cohort/reference-confidence/GT-GQ 等 Kokkos kernel，encode 负责 GATK-compatible 文本
格式化，sink 仍是唯一 HTSlib writer。队列按 record byte estimate 背压，manifest 暴露
三段 item/byte/peak；OpenMP/Serial 的 GenotypeGVCFs contract 与 GATK oracle 均通过。
输入 VCF 的全量 ALT union/跨样本合并仍在 pipeline 之前完成，因此这一步是输出阶段
的并发/确定性增量，不宣称已经把 joint staging 的峰值内存降为 tile 级。

本轮新增 GenotypeGVCFs 显式 `--stream-by-locus`：先用轻量的 concrete-variant
span probe 建立跨样本 reference-block 切分边界，再为每个输入建立一个 HTSlib
cursor；Host priority queue 只保留每个 shard 的当前 `Record`，同一 locus 的
concrete ALT union、`<NON_REF>` projection、sample FORMAT merge 完成后立即送入
统一 `ThreeStagePipeline` 并释放 `bcf1_t`。因此 joint staging 不再创建全基因组
`std::vector<Record>`，当前 locus 的队列只按 `O(input_shards + one_locus + queue_bytes)`
增长；span index 本身按 concrete variant 数量增长，但 block segment 通过 cursor
状态逐段生成，不会为长 block 创建整段 `Record` vector。`streamed_loci`、
`stream_max_inflight_records`、`streamed_peak_host_bytes` 和 `stream_span_probe_records`
写入 manifest。该路径与 aggregate 在 chr17 fixture、
多样本/反向 ALT 顺序、跨样本 reference block split 以及
`--include-non-variant-sites` 的 dense block 上逐点 lazy emission 均在解压后
与 aggregate exact；dense stream 要求 indexed `-R`，避免从 block 起始 REF 猜测
内部碱基。
`benchmark_genotype_gvcf.py --stream-by-locus` 的 128-locus 双 shard 回归在当前
OpenMP 构建 p50 约 0.139 s、`stream_max_inflight_records=4`、峰值 Host staged bytes
约 12.5 KiB；Serial p50 约 6.3 ms。该数字只用于生命周期/资源合同回归。

CombineGVCFs 现在增加显式 `--stream-merge`/`--stream-by-locus`：预读所有输入
header 后，每个输入只保留一个 decoded `Record`，用 Host priority queue 做稳定
coordinate k-way merge；同一 locus 的 ALT union、PL/AD remap、FORMAT sample merge
和相邻 REF/<NON_REF> block coalescing 仍沿用原有语义，唯一 writer 在 flush 前
完成 GT/GQ materialization。queue 按 Record byte estimate 背压并写出
`max_inflight_records`/`peak_queued_bytes` telemetry；无样本、多样本和
`--fastgatk-materialize-genotypes` fixture 均与 aggregate 输出 exact。现在
`--convert-to-base-pair-resolution`/`--break-bands-at-multiples-of` 也使用一个
原始 block + 当前 segment 的 lazy iterator，segment 数量不再在单输入 cursor
中预先展开，并由 `lazy_reference_blocks/segments` telemetry 审计。
重复 `-L/--intervals` 还支持 `--interval-set-rule UNION|INTERSECTION`：Host
先对每个 selector 归一化为半开区间集合，再按选择的集合运算过滤输入；aggregate
和 `--stream-merge` 使用同一结果，规则写入 summary/OutputManifest，并由真实
indexed fixture 覆盖交集只保留重叠位点。

SelectVariants 已补齐重复 site-expression 的 OR 语义、重复 `-L/--intervals`
区间筛选（按记录 span/`END` overlap）、任意样本 GQ/DP/AD/MIN_DP 比较和常用 genotype predicate，以及
contig/POS/REF/ALT-key concordance/discordance 子集；常用 genotype
`getAD()[i]`/`getPL()[i]` 与 `isHom` 以及 VariantContext
site methods（类型、transition/transversion、过滤状态、等位基因计数、INFO presence、坐标、type string）已显式
实现；常用字符串/正则（含 `=~`/`!~`）/INFO-null/VariantContext 方法已扩展，并支持 `--select` 中显式 `vc.getGenotype("S").*` receiver 与布尔方法的显式 `== true/false` 比较；完整 JEXL
方法、样本级 concordance 和 reference-backed interval 语义仍需单独门禁。
INFO 数组字段现在同时接受 GATK/HTSJDK 的
`vc.getAttribute("TAG").get(i)`、历史点索引 `vc.getAttribute("TAG").i`
以及兼容性简写 `TAG[i]`，按声明顺序读取 Number=A/Number=R 数值向量；
genotype AD/PL 也支持 `getAD().i`/`getPL().i`，并由 VariantFiltration 与
SelectVariants 的 native/GATK fixture 回归覆盖。
`--exclude-non-variants` 现按 `VariantContext.isVariant()` 语义丢弃仅含
hom-ref/no-call 的样本位点（无样本时保留带 ALT 的 annotation-only record），
`--min-indel-size`/`--max-indel-size` 按 REF/ALT 绝对长度差筛选 indel；纯选择
路径保留已有 AC/AN/AF INFO，只有 sample/allele subset 时才重算，三条边界均有
bundled GATK 4.6.2.0 oracle 和 manifest telemetry。

HTSlib BAM/CRAM Host reader 现在也统一处理 HC、Mutect2、BaseRecalibrator 的
重复 `-L/--intervals`/`--region`：literal、Picard interval-list/interval 和
BED（plain/`.gz`，忽略 UCSC `track`/`browser` header）均转换为零基半开区间，重叠范围先归一化，文件/记录计数写入 manifest；
HC 的 gVCF reference-confidence 输出消费同一组 disjoint intervals。

LeftAlignAndTrimVariants 已在 multiallelic split 后重算 GQ 与 AC/AN/AF，并支持
`--keep-original-ac` 写出 AC_Orig/AF_Orig/AN_Orig；复杂 symbolic、interval/cloud
和 bit-identical normalization 仍保留 explicit fallback。

GatherVcfs 的 `-L/--intervals` 已接受 literal、Picard `.interval_list`/`.intervals`/
`.list` 以及零基半开 `.bed`/`.bed.gz` 文件（忽略 UCSC `track`/`browser` header），按文件格式转换为统一的零基半开 Host interval，并按记录 span（含 gVCF `END`）做 overlap 选择，
并在 OutputManifest 记录输入文件与区间条数；dispatcher 已负责 HTTP(S)/presigned
URI 输入 staging，远程输出可在显式 upload mode 下按 sidecar→主文件顺序提交，默认
仍 fail-closed。`-CO/--COMMENT`、`--CREATE_INDEX` 和 `--COMPRESSION_LEVEL` 已
支持重复 selector 的 GATK `-isr/--interval-set-rule`（UNION 默认、INTERSECTION
可选），并在 summary/manifest 中记录实际集合规则；
对齐 Picard/GATK 的 VCF header/index/BGZF 参数边界，`benchmark_gather_vcfs.py`
提供 shard/record 吞吐基线。由于 HTSlib native writer 生成 BCF 2.2、而 pinned GATK/htsjdk 的
GatherVcfs 兼容边界是 BCF 2.1，native `.bcf` 输出现在明确
`BACKEND_UNAVAILABLE`，由 dispatcher 的显式 `--fallback` 路由到 Java；不会发布
下游 GATK 无法回读的文件。

VariantFiltration 已补齐重复 `-L/--intervals` 区间筛选（按记录 span/`END` overlap）、site/genotype
filter-expression inversion、按输入记录 span/`END` 计算的 VCF/BCF mask interval/extension、cluster-size/window
过滤，和现有 FORMAT/FT 规则一起写入 manifest；常用 VariantContext site
 methods（类型、transition/transversion、过滤状态、等位基因计数、INFO presence、坐标、type string、
contig/ID/reference/alternate 字符串、字符串方法、regex（含 `=~`/`!~`）、enum、null checks）已显式实现；布尔 site/genotype/allele 方法也支持显式 `== true/false` 比较；
`isPolymorphicInSamples`/`isMonomorphicInSamples`、`hasGenotypes()`、
`getCalledChrCount()`/`getNoCallCount()`、`getHomRefCount()`/`getHetCount()`/
`getHomVarCount()`、以及基因型 `isAvailable()`、`getAlleles().size()`、
`getFilters().size()`/`getFilters().contains()` 也已接通；
SelectVariants/VariantFiltration、LeftAlignAndTrimVariants、ReblockGVCF、
CombineGVCFs 和 GenotypeGVCFs 的 `-L/--intervals` 现在共享 Host
interval-list/BED 解析器，按记录 span/`END` overlap 选择，并在 manifest 记录解析计数；
ReblockGVCF、CombineGVCFs、GenotypeGVCFs、SelectVariants、VariantFiltration 和
LeftAlignAndTrimVariants 的重复 selector 还支持 GATK
`--interval-set-rule UNION|INTERSECTION`，在记录过滤/规范化前完成确定性半开区间集合运算；
CombineGVCFs 的 FORMAT PL/AD 投影和任意 ploidy GT/GQ 现在也复用共享 Kokkos
genotype APIs，manifest/benchmark 记录独立 kernel telemetry；
其标准 `--call-genotypes` 已加入（历史 `--fastgatk-materialize-genotypes` 保留为
native 实验模式），
并补齐 `--convert-to-base-pair-resolution`/`--break-bands-at-multiples-of` 的
reference-block 分段：按 1-based breakpoint 保留边界、禁止后续 merge 跨界合并，
提供 `-R` 时对分段内部 REF 碱基走 FAIDX；无参考索引时保持输入 REF 并显式记录参数。
常用数值 `AS_*` `Number=A` 过滤可通过 `--apply-allele-specific-filters` 写出每个 ALT
的 `AS_FilterStatus`；复杂 allele-specific 和完整 JEXL 仍待实现。
普通 INFO 数组也支持 GATK 的 `vc.getAttribute("TAG").get(i)` 与 `TAG[i]` 数值比较，
缺失或越界元素保持 missing-values contract；VariantFiltration 还支持带 Java
优先级的标量/INFO/genotype 数值算术（`+`、`-`、`*`、`/`、`%`），并由 GATK
4.6.2.0 oracle 覆盖；SelectVariants 对 HTSJDK 拒绝的向量算术保持 fail-closed。
完整 JEXL/复杂 allele-specific 仍保持 explicit fallback。

新增独立 `verify_variant_filtration_gatk_oracle.py`：固定多样本 VCF 由 GATK
4.6.2.0 Java 与 native 逐字段对照 QUAL/INFO/FORMAT/GT/FT/FILTER，覆盖
INFO/QUAL site predicate、compact genotype FORMAT predicate 及
`--invert-filter-expression`；比较器只归一化 INFO 数值短格式、FORMAT key 顺序、
尾部 missing FORMAT 与 FILTER label 顺序，避免 summary-only 假 smoke。native
同时修复 site filter append 时保留已有非 PASS FILTER 的链式语义；OpenMP/Serial
oracle 均 exact。完整 JEXL、annotation/plugin 与 cloud 语义仍明确 fallback。

HC/Mutect2 的本地 haplotype 组合上限已由
`--max-haplotype-combination-alleles` 显式控制（1..16），低于等于上限时保留
完整 `2^N` 组合，超过上限时按稳定坐标分块并在每个块内保留完整 `2^N`
组合；同时共享跨块 REF haplotype，并对相同 read×haplotype 的 PairHMM/SW
请求做稳定去重。候选×haplotype link budget 会在极端区域自动降低块宽，避免
资源上限隐藏在 Host 常量中；分块/去重计数进入 manifest telemetry。

PairHMM 的候选分组现在消费同一份 Host AssemblyRegion（含 halo）坐标：每个
region 独立构造 local reference、组合 haplotype、graph path 和 read×haplotype
请求，禁止跨 region 组合；路径坐标必须与 region 相交才进入该组。直接 API 或
graph projection 产生的 region 外候选仍保留一个显式 contig fallback，并由
`pairhmm_assembly_region_groups`、`pairhmm_unassigned_candidates` 和
`pairhmm_assembly_region_partitioned` 写入 manifest。分离 region 的跨边界回归已
加入 `fastgatk-calling-indel-smoke`；完整 GATK AssemblyRegionWalker 的 active
probability、downsampling、全量 assembly 参数仍按后续 oracle 继续补齐。

PairHMM request 现在在 Host 侧按当前 AssemblyRegion 的 CIGAR 投影执行
hard-clip：只保留 region 内首尾 reference-projected base 之间的连续 read slice，
并同步裁剪 base quality、BI/BD indel quality；裁剪和无法形成合法 slice 的数量由
`pairhmm_reads_clipped`/`pairhmm_reads_dropped_after_clipping` 写入 telemetry。
该路径已由 chr17 GATK oracle 断言，仍不等价于 GATK 完整 read-haplotype CIGAR、
uncertainty 和所有 soft-clip edge 语义。

PairHMM 输出后的 allele max-marginalization 与 read×haplotype BestAllele
reduction 已统一到 Kokkos `RangePolicy` API；`reduce_read_allele_best_kokkos`
返回 best/second likelihood 及 concrete allele owner，Host 仅保留
CIGAR/graph eligibility、稳定 row 构造和 BestAllele margin policy。regular/flow
调用共享同一数值归约，HC/Mutect2 多等位 depth/REF/ALT count 采用一次
`(locus, read)` owner，避免独立 biallelic 计数重复计深度，并在 API smoke、kernel
benchmark 和 HC oracle 中校验 determinism/telemetry；完整 haplotype posterior 及跨 AssemblyRegion 语义仍保持为
后续兼容性门禁，flow HMER 的双向 uncollapse、canonical sequence-remap 和
identical-group ownership 已在下条落地。
当前已先吸纳 GATK `LongHomopolymerHaplotypeCollapsingEngine.collapseBases()` 的
Host 等价物：`collapse_flow_homopolymers` 在 flow key 编码前执行首个 run 保留、后续
run 按阈值截断，并由 HC/Mutect2 选项与 telemetry 暴露；`-1` 复现 GATK 的
read-group `mc` 自动阈值选择。shared calling library 现在还提供按 GATK
双向 SW/CIGAR 规则的 `uncollapse_flow_hmers`，支持 partial/limit-to-threshold
模式；flow batch materialize post-likelihood restoration view 并写入
`pairhmm_flow_haplotypes_uncollapsed`。随后按恢复后的 bounded
`(tid,start,end,sequence)` 建立 canonical group，合并 candidate ownership/mask，
并让 Kokkos max-marginalization 消费 canonical entry，避免 identical collapsed
haplotype 重复计入；`pairhmm_flow_haplotype_remaps` 和
`pairhmm_flow_identical_haplotype_groups` 进入 HC/Mutect2 manifest。跨
AssemblyRegion posterior 仍是后续门禁。
启用 collapse 前另有 reference-window `needsCollapsing` gate，避免只因 alternate
haplotype 的长 run 就改变没有长 reference hmer 的普通区域。

候选进入 PairHMM 前还会对 anchored REF span 做同一份 region 边界判定：跨越
两个 disjoint AssemblyRegion 的 compound/graph allele 被显式抑制，region 外
候选保留为 contig fallback，并分别由 `assembly_cross_region_candidates` 与
`assembly_unassigned_candidates` 计数；这使跨 region 的局部图不会被误认为完整
haplotype，而不会把纯 CIGAR indel 因 activity 阈值而静默丢弃。

downsampling 也已从 kernel-only 参数提升到 HC/Mutect2 Host contract：
`--max-reads-per-locus`（兼容 `--max-reads-per-alignment-start`）和
`--downsampling-seed` 写入共享 `ReadFilterOptions`，在同一 mask 上消费 assembly、
activity、graph、SW/PairHMM 与 RCM；`downsampled_reads`、cap 和 seed 进入
manifest/stats/benchmark，并由重复运行的同 seed regression 验证稳定性。

### 4.10 MarkDuplicates / SortSam 外部内存模块（P2）

**C++ Host**

- 根据 ResourceBudget 生成 sorted runs；read-name grouping、duplicate key、
  representative selection、metrics、spill manifest 和多路 merge 全部由 Host
  控制。
- local SSD 优先；spill 采用 checksum、临时文件命名、失败清理和 checkpoint。
- MarkDuplicates 的 checkpoint 为版本化 Host 文本 sidecar：记录输入
  size/mtime/options fingerprint、effective group limit、首遍 counters/library
  metrics 和完整 spill-run index；`--resume-spill` 只在签名、run 完整性均通过后
  重用首遍状态，成功发布 output/metrics/index/manifest 后原子清理 sidecar 与 runs。

**Kokkos kernel**

- 批量抽取 coordinate/read-name/library/flags key；Kokkos sort/radix sort、key
  compare、压缩和 merge head buffer；不在 GPU 保存完整 read object。
- kernel 结果只给出 stable key/flag，duplicate 语义和输出顺序由 Host 决定。

**验收**

- duplicate flags、metrics、sort order、SAM/BAM index、内存上限、spill/restart；
  其中 restart 回归会注入首遍完成后的中断，验证 SLURM/Nextflow retry 不重新扫描
  输入且结果/metrics 与一次性运行一致；
  不能以“GPU sort 更快”替代语义验证。MarkDuplicates 的代表选择固定为
  Picard `SUM_OF_BASE_QUALITIES`（不混入 MAPQ），非配对 duplicate 不计 optical；
  默认 `ADD_PG_TAG_TO_READS=true`，有 pinned GATK/Picard jar 时，回归同时比较
  duplicate flag、`DT` tag 与标准 DuplicationMetrics 字段/duplicate-set
  histogram（估算 library size 排除 optical pairs）。SortSam 的 coordinate/queryname tie-break 也按
  HTSJDK comparator（strand、flag、MAPQ、mate、TLEN、HI）回归，且 benchmark 覆盖
  spill/k-way merge。

### 4.11 VQSR / CNV 数值模块（P2）

`DepthOfCoverage` 的 coverage-walker 边界同样遵循统一生命周期：Host/HTSlib 负责
read filter、CIGAR 投影和 sample partition，`HostBatch -> KernelPlan.prepare ->
Kokkos Views -> execute -> collect` 的 `RangePolicy` atomic kernel 负责逐 locus
depth/base/deletion 计数，Host 负责 GATK CSV/summary/histogram sidecar 写出。
manifest 与 benchmark 必须记录 execution space、batch/observation 数以及
prepare/execute 时间；native `COUNT_FRAGMENTS` 提供稳定 read-name+locus
去重扩展（固定 GATK 4.6.2.0 目前会拒绝 fragment counting，因此不宣称
bit-identical），gene、cloud 和无界 whole-reference 物化继续显式 fallback。

DepthOfCoverage 的区间边界现复用 HtsReader 的 include/exclude 语义：重复
`-L` 可用 `--interval-set-rule UNION|INTERSECTION` 组合，`-ip/-ixp` 分别扩展
include/exclude 区间，`-XL/--exclude-intervals` 在 Host 投影前切分目标 span。
排除后的 span 才进入 Kokkos dense locus Views，因而主 locus 表、interval
summary、quantile/cumulative 分母不会保留伪造的零深度行；双 backend 与 GATK
4.6.2.0 的排除/padding 输出已加入 contract oracle。

**C++ Host**

- VQSR：模型参数、收敛、随机种子、annotation schema、tranche/report writer。
- CNV：HDF5/TSV 读写、interval/sample metadata、segmentation/MCMC 状态和模型。

**Kokkos kernel**

- VQSR VBEM/GMM expectation/maximization、批量评分使用 Kokkos Kernels dense math。
- CNV CollectReadCounts 用 interval overlap/histogram；Denoise 用 Kokkos Kernels
  SVD/GEMM；ModelSegments/gCNV 的复杂 segmentation 首版 Host。
- 采用 block-wise matrix，控制 device memory；Strict 模式固定 reduction 和
  NaN/underflow 规则。

**验收**

- VQSR recal/tranche、CNV HDF5/TSV、SVD tolerance、segment boundaries、模型重启
  和跨 backend deterministic/fast 差异报告。

`CollectReadCounts` 已先落地 CNV 的 TSV/native 边界（native 默认 HDF5，与 GATK
默认一致；需文本交换时显式指定 `--format TSV`）：Host 解析 BAM/CRAM
interval 与 sequence dictionary，Kokkos 对 read-start 和 `-XL/--exclude-intervals`
切分后的 target 做批量 overlap 计数，输出带 SAM-style header、sample metadata 和
`CONTIG/START/END/COUNT` 的 GATK SimpleCountCollection TSV。默认过滤遵循
GATK 的 mapped/non-duplicate/MAPQ>=30 规则，`--include-duplicates` 可显式覆盖，
并接入 OutputManifest、CTest 和 dispatcher。GATK 4.6.2.0 的 CNV 参数校验要求
`OVERLAPPING_ONLY`，且拒绝 `ALL` 与非零 `-ip/-ixp`；native 对这些输入 fail-closed。
manifest/JSON telemetry 同时保留 requested、excluded interval 数与实际 merge 数，
确保重复、重叠、排除和相邻区间的边界可审计。
每个 read batch 的计数阶段统一走 `HostBatch -> KernelPlan.prepare -> Kokkos Views
-> execute -> collect`，并记录 `RangePolicy`、execution space、batch/record 数与
prepare/execute 时间；interval metadata view 在 batch 间复用，避免 accelerator
路径反复传输固定字典。
现在还可选写入 GATK
`HDF5SimpleCountCollection` 的 sample/dictionary/interval/counts schema；native
现在保留 HTSJDK `SAMTextHeaderCodec` 的 `@HD VN:1.6` 与全部 normalized `@SQ`
标签（包括 AS/M5/UR/SP），并以 Java GATK 4.6.2.0 生成的同输入 HDF5 做
dictionary/metadata 与 Java reader round-trip gate。该 slice 只证明
SimpleCountCollection 的 sample/dictionary/interval/counts 文件边界，不能把
TSV、局部 schema 或它误报为完整 CNV/HDF5 SVD-PoN 等价；cloud indexed
metadata 与 release-specific HDF5 byte identity 仍是 explicit fallback。

`DenoiseReadCounts` 已接入同一 TSV/HDF5 边界：Kokkos map kernel 执行 fractional
coverage、正 sample median、epsilon-floor `log2` 和 centered log median，并支持
interval-matched 四列 TSV `--panel-of-normals` 的 sample/PoN fractional-coverage
ratio；HDF5 输入和单样本 HDF5 PoN 使用同一 GATK SimpleCountCollection
schema；严格拒绝 all-zero sample 和未对齐 PoN。宽表 PoN
(`CONTIG/START/END/SAMPLE...`) 进一步通过 Kokkos GEMM-style `X^T X` 和确定性
power-iteration SVD，按 `--number-of-eigensamples`/`--svd-rank` 做低秩投影去噪，
并在 manifest 中记录 rank/样本数；所有 fractional/log2、Gram、投影和重构阶段
统一走 `HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect`，
以 `RangePolicy`/`MDRangePolicy` 执行，并记录 kernel batch/observation、execution
space 与 prepare/execute 时间。CTest、dispatcher contract 与
`benchmark_denoise_read_counts.py` 文件边界 benchmark 已接入。GATK-style
`GC_CONTENT` 的 101-bin correction 已接入；提供 PoN 时按 GATK 忽略
`--annotated-intervals`（除非 HDF5 PoN 自带 GC 模型），并在 manifest 记录实际是否
使用 GC correction；cloud metadata 仍显式 fallback。

`CreateReadCountPanelOfNormals` 现在补齐了 Denoise 的上游 PoN 边界：重复 TSV 或
`AnnotateIntervals` 已用 Host/FAIDX + Kokkos 实现 GC_CONTENT annotated-interval TSV，
并支持 htsjdk BEDCodec 坐标、score/NaN 规则和长度加权的 MAPPABILITY/
SEGMENTAL_DUPLICATION_CONTENT tracks；track 按 contig 建索引并在 RangePolicy 内二分查找，
所有 annotation 数值统一走 RangePolicy，
并由下游按 GATK 101-bin exponential GCBiasCorrector correction 消费。
`HDF5SimpleCountCollection` 输入经过 GATK 顺序的 fractional-coverage、GC correction、interval/sample
zero 过滤、interval fractional median、零值中位数填补、极端值截断和 centered
log2 标准化；Kokkos Gram/power-iteration SVD 生成
`HDF5SVDReadCountPanelOfNormals` v7 所需的 original/panel/interval/eigensample
chunk paths，并按最大 chunk cell 数分行写出；Gram 与 left-singular-vector projection
统一走 `HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect`，同时
记录 `MDRangePolicy`/`RangePolicy`、execution space 和 kernel timing。native HDF5 PoN 已由 Java GATK 4.6.2.0 `DenoiseReadCounts` 直接读回；
GC annotated intervals 已纳入 101-bin correction，cloud/Spark 和完整 release-specific
report 语义继续显式 fallback。契约见 `fastgatk-native/scripts/verify_create_read_count_panel_of_normals.py`。
文件边界 benchmark 见 `fastgatk-native/scripts/benchmark_create_read_count_panel_of_normals.py`。

`CountBasesInReference` 已补齐同一 reference/fasta walker 的最小边界：Host/FAIDX
负责 FASTA、literal、Picard interval-list 与 BED interval（含 `.gz` 压缩文本），
重叠区间先做集合合并，Kokkos `RangePolicy` atomic histogram 统计原始字节，stdout/`-O` 报告按 GATK
ascending-byte 格式输出；OutputManifest、Java 4.6.2.0 synthetic oracle、CTest、
dispatcher contract 和 benchmark 已接入。压缩 interval 文件、cloud reference 与
release-specific logging 仍显式 fallback。

`CompareReferences` 与 `CheckReferenceCompatibility` 已补齐 `tools/reference` 的
核心 CPU/I/O 边界：前者按相邻 `.dict` M5 或重算 MD5 复现 MD5-keyed 表、pair status
和 `FIND_SNPS_ONLY` 的 Kokkos mismatch mask；后者读取 BAM/CRAM/SAM 或 VCF/BCF
字典，复现 COMPATIBLE/COMPATIBLE_SUBSET/NOT_COMPATIBLE 表。两者均接入
GATK 4.6.2.0 fixture oracle、CTest、dispatcher、OutputManifest 与 benchmark；
FULL_ALIGNMENT/MUMmer、cloud 和完整 htsjdk/Tribble warning 仍显式 fallback。

`FastaReferenceMaker` 与 `FastaAlternateReferenceMaker` 已补齐 `tools/walkers/fasta`
的常用 reference 输出边界：Host/FAIDX 按 dictionary 顺序合并 literal、interval-list
和 BED 区间，复现数字 sequence ID、description、uppercase sequence、line-width 和输出 `.fai`/Picard `.dict`（M5）；Alternate
路径按 GATK 规则应用 simple SNP/anchored indel、SNP mask、mask priority 与 diploid
IUPAC。固定字节复制统一经过 Kokkos `RangePolicy` 的
`HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect` 生命周期，
并接入 GATK 4.6.2.0 oracle、CTest、dispatcher、OutputManifest 与文件边界 benchmark；
复杂 symbolic/overlapping indel、cloud 和 release-specific feature validation 仍显式
fallback。另已修正 `--use-iupac-sample` 的纯合投影：按 GATK
`getIUPACBase` 选择 hom-ref/hom-ALT 的 genotype allele（多 ALT 的 2/2 不再误用
ALT[0]），并由 pinned 4.6.2.0 Java/native oracle 与 OpenMP/Serial CTest 固定。

`ShiftFasta` 已补齐同一 FASTA walker 包的循环参考边界：Host/FAIDX 读取每个 contig，使用
Kokkos `RangePolicy` 按指定或默认半长 offset 旋转序列（该 byte permutation 为带宽受限
操作，故不强行套用 `Kokkos::Experimental::simd`），输出 GATK FASTA wrapping、`.fai`/`.dict`、
UCSC chain 以及 regular/shifted interval sidecar；OutputManifest 记录 execution space 和
prepare/execute telemetry。GATK 4.6.2.0 oracle、默认偏移双移 round-trip、CTest、dispatcher
和文件边界 benchmark 已接入 `verify_shift_fasta.py`/`benchmark_shift_fasta.py`；远程/cloud
存储及 release-specific validation 仍显式 fallback。

`IndexFeatureFile` 已加入 feature-file compatibility 边界：Host 用 HTSlib 检测格式，BGZF
VCF/GVCF/BED 调用 `tbx_index_build3` 生成 `.tbi`，BCF 调用 `bcf_index_build3` 生成 `.csi`，
并以 HTSlib decoder threads 控制 I/O；Kokkos `RangePolicy` 仅承担轻量元数据投影。未压缩
VCF 均由 Host 写出 HTSJDK LinearIndex v3；`.g.vcf` 使用 128000 bp 起始 bin、
occupied-contig 自适应合并、URI/size/mtime/byte offsets，并以 pinned GATK 4.6.2.0
corpus 做逐字节和实际区间查询 oracle。普通 VCF 按 HTSJDK DynamicIndexCreator 的
FOR_SEEK_TIME 评分选择稀疏 2000 bp linear candidate 或密集 75-record interval-tree，
空文件保留空 interval-tree；索引均可由 HTSJDK 查询。interval-tree 遍历顺序和稀疏
linear 的 leading-empty block 已在 pinned one-record VCF 上逐字节复现；dense
interval-tree 现复现 HTSJDK 红黑树插入/旋转后的前序遍历，并在 pinned dense corpus
上逐字节一致；未压缩 BED 也按同一 DynamicIndexCreator/LinearIndex v3 路径生成 `.idx`，
并以 GATK BED fixture 做逐字节 oracle。其它未压缩 Tribble
dynamic index、未知 codec 仍明确 fail-closed，dispatcher 通过 `--fallback` 保留 GATK 原始语义；contract、
OutputManifest 和文件边界 benchmark 已接入。

`CountReads` 与 `FlagStat` 已加入读文件 QC 边界：HTSlib Host 以 `--batch-records` 做有界
SAM/BAM/CRAM 流式读取，CountReads 的 Kokkos `RangePolicy` reduction 和 FlagStat 的
12-counter atomic kernel 统一经过 `HostBatch -> KernelPlan.prepare -> Kokkos Views ->
execute -> collect`。两者支持 `-R/--reference` 注入 CRAM reference、重复 `-L/--intervals`
和 `-XL`（有索引走 HTSlib iterator，无索引按 CIGAR span 顺序过滤；literal、interval-list/BED
选择器含 `.interval_list.gz/.bed.gz` 时由 HTSlib 解码），以及
`-RF/--read-filter`/`-XRF/--inverted-read-filter`/`-DF/--disable-read-filter` 的 26 类原生过滤器集合（包括配对、
方向、MAPQ 可用性、CIGAR、RG 与 fragment-length）；dispatcher 通过每个 registry 条目的
`native_read_filters` 显式声明该工具能力，未知类仍 fail-closed 并保留完整 argv 给
`--fallback`；FlagStat 同时检查
坐标缺失造成的 read/mate unmapped 状态，保持 HTSJDK 语义；文本输出、OutputManifest、
GATK 4.6.2.0 chr17 oracle、OpenMP/Serial contract、dispatcher 和 benchmark 已接入。
其中 `ReadLengthReadFilter` 的 `--min-read-length`/`--max-read-length` 已接入
HC、Mutect2、CountReads、FlagStat 的共享 typed Kokkos mask，并有 0-record
边界回归；`--disable-tool-default-read-filters` 仅清除隐式默认集，显式
`--read-filter` 保持生效并写入 HC/Mutect2 manifest；其余参数化 read-filter 仍按名称 fail-closed。
重复 `-I`/`--input` 现在按 GATK ReadWalker 语义逐 shard 顺序聚合，CountReads 的记录数与
FlagStat 的 12 项计数/百分比均跨输入累加；共享 Kokkos staging buffer 在输入之间复用，
OutputManifest 保留单输入 `input` 字段并新增 `inputs` 与 `input_count`。OpenMP/Serial
以及 GATK 4.6.2.0 重复输入 oracle 均通过，dispatcher registry 已声明该选项为 repeatable；
重复 `-L`/`--intervals`/`--region` 还支持 GATK `-isr/--interval-set-rule UNION|INTERSECTION`，
其中单个 interval-list/BED 文件内部作为一个 selector set，再与其它 selector 集合做交集；
该规则写入 OutputManifest，并由 OpenMP/Serial 与 GATK 4.6.2.0 11-record/12-field
intersection oracle 覆盖；
同时新增 `-ip/--interval-padding`、`-ixp/--interval-exclusion-padding`（按 contig 边界
clamp）和 `-imr/--interval-merging-rule ALL|OVERLAPPING_ONLY`。索引路径使用 HTSlib
multi-region iterator，在 OVERLAPPING_ONLY 的相邻/重叠 fragment 中对同一 alignment
去重，保持 GATK traversal 计数语义；include/exclude padding、manifest、非法 enum/负值
校验以及双 backend Java oracle 均已覆盖；
`ReadNameReadFilter` 与可重复 `--read-name` 也已接入，并强制 GATK 同等的 filter/参数
依赖检查；定点 read-name 的 CountReads/FlagStat 结果与 Java oracle 一致。
`ReadGroupReadFilter` 与可重复 `--keep-read-group` 也已接入，直接读取 BAM `RG` tag，
并覆盖相同的依赖检查与 Java oracle。
数值型 `ReadTagValueFilter` 也已接入，支持 NM/MQ/SM 等两字符 auxiliary tag 的六类
比较运算及 GATK 同等的参数依赖/非法值 fail-closed 校验。
`-XRF` 对同一原生 predicate 应用 GATK `InvertedReadFilter` 语义，并保持默认 filter
及普通 `-RF` 先行；pinned chr17 的参数化 MAPQ 反向边界在 OpenMP/Serial 均与
GATK 4.6.2.0 的 FlagStat 12-field 输出逐行一致。

`GetPileupSummaries` 的 LocusWalker filter mask 也支持 `-XRF/--inverted-read-filter`
（当前 14 个显式 native predicate）；默认和普通 `-RF` 先筛选、再运行 inverse，避免把
反向 filter 误作默认集的移除。Pinned Java 低-MAPQ pileup-table oracle 在 OpenMP/Serial
均逐字节一致，反向 filter count 与能力标记写入 manifest。
`ReadGroupBlackListReadFilter` 也已接入，支持重复 `--read-group-black-list ATTR:VALUE`
并通过 HTSlib 查询 `@RG` header 属性。
无 `-L` 的 `-XL` 也按 whole-reference minus excluded span 处理，跨越未排除区域的
alignment 不会被误删。
完整 Java/cloud dictionary 语义仍保持 prototype。
HC 输出也已收敛 provenance 开关：显式 `--add-output-vcf-command-line false` 时，native
不再把 `##source` 和 `##fastgatk_*` 诊断头写入 VCF/GVCF，普通 VCF 完整字节序列在
GATK 4.6.2.0 broad fixture 上 exact；GVCF reference-confidence block 的全量 RCM
数值/分块差异仍保持显式非 bit-identical 边界。
HC/Mutect2 还共享 GATK `--dont-use-soft-clipped-bases`：关闭 terminal/internal
soft-clip 派生的组装证据，但不屏蔽显式 CIGAR I/D；该策略在 shared calling
options 中统一，并进入 stats/manifest telemetry。
Mutect2 同时接通 GATK `--min-base-quality-score`/`-mbq`（兼容
`--min-base-quality`），默认 Q10，质量门在 assembly/activity/PairHMM 前统一生效。
`GetPileupSummaries` 同样接受 typed `ReadLengthReadFilter` 的
`--min-read-length`/`--max-read-length`，并在计数前复用该边界与 manifest telemetry。

`SplitIntervals` 已加入 scatter 编排边界：Host 读取 FASTA `.fai/.dict`，解析 literal、
interval-list/BED（含 `.interval_list.gz/.bed.gz`，压缩文本由 HTSlib 解码）与排除区间，按 base/count 或保持完整 interval 的模式做确定性分片；无显式
`-L` 时支持 `--min-contig-size` 过滤短 contig；同时支持 `--dont-mix-contigs`、稳定文件命名和 OutputManifest；总 bases 统计统一走
`HostBatch -> KernelPlan.prepare -> Kokkos View -> RangePolicy -> collect`，并记录
execution-space/policy/计时 telemetry；Picard dynamic scatter state machine 的五种
mode shard body 与 GATK 4.6.2.0 oracle 对照，并接入 OpenMP/Serial CTest、dispatcher、local
SLURM/Nextflow pipeline artifact 验证、`verify_split_intervals.py` 和
`benchmark_split_intervals.py`。

`FilterIntervals` 已承接 annotated-interval 和 TSV/HDF5 count 过滤：Host 解析
interval-list/BED/literal、GATK `OVERLAPPING_ONLY` merging/`UNION` interval-set rule、排除区间和 count collection，Kokkos mask kernel 对
`GC_CONTENT`、`MAPPABILITY`、`SEGMENTAL_DUPLICATION_CONTENT` 执行 GATK inclusive
bounds，并接入 low-count 与 Apache Commons Math 百分位 extreme-count 规则；Host
再执行 solitary-contig protection 并写 Picard interval-list/OutputManifest。现在接受
GATK 短别名 `-imr/-isr/-ip/-ixp`，并对默认/显式非 `OVERLAPPING_ONLY`、非零
`-ip/-ixp` 以及 `-isr INTERSECTION` fail-closed（GATK `FilterIntervals` 参数校验）；多 contig UNION 按
sequence-dictionary 顺序稳定归一化。annotation 和 count contract 已与 GATK 4.6.2.0 对照，并接入
`benchmark_filter_intervals.py`；native GATK SimpleCountCollection HDF5 已通过共享
HDF5 reader 接入；cloud/object-store 语义仍显式 fallback。

`verify_filter_intervals.py`、`benchmark_filter_intervals.py` 已纳入 `verify_all.sh`，覆盖
短别名、参数拒绝、字典顺序和 Java oracle。

`PreprocessIntervals` 已接入 CNV 上游：Host/FAIDX 读取可选 interval-list/BED/literal
（省略 `-L` 时按 whole-reference 处理），按 `OVERLAPPING_ONLY` 合并，并在 padding/
binning 前执行重复 `-XL/--exclude-intervals` subtraction；GATK 要求的
`--interval-exclusion-padding 0` 仍显式校验，非零值 fail-closed。之后采用 GATK
非重叠 padding/midpoint 规则生成 bins；Kokkos bounded batch kernel 删除 all-N bins，
输出可直接被 `CollectReadCounts` 使用的 Picard interval-list。contract 见
`fastgatk-native/scripts/verify_preprocess_intervals.py`（含 Java `-XL` oracle），
文件边界 benchmark 见 `fastgatk-native/scripts/benchmark_preprocess_intervals.py`；
cloud reference 和 release-specific interval 语义继续显式 fallback。

`CallCopyRatioSegments` 已加入该 CNV TSV 链路的下游：Kokkos
`Experimental::simd<double>` 分组计算 `2^log2` copy ratio（SIMD
width/groups 与 execute time 写入 manifest；非 SIMD execution space 由 Kokkos
自动退化为标量），该批处理统一走 `HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect`
生命周期和 `RangePolicy`，并记录 policy、batch/observation、prepare/execute telemetry；Host 按 GATK `SimpleCopyRatioCaller` 的 inclusive neutral bounds、两遍
length-weighted statistics 和 outlier/calling z-score 生成 called segment TSV，并对
singleton/空 neutral 集合进行显式数值契约处理，
并可输出 IGV legacy `.seg`，接入 CTest、dispatcher contract 和文件边界
benchmark。确定性的 TSV segmentation boundary 已由 `ModelSegments` prototype
接入；完整概率模型、HDF5/PoN/GC 状态仍显式 fallback。

`CollectAllelicCounts` 已加入 CNV 输入链路：HTSlib Host reader 解析 BAM/CRAM/SAM
和 FASTA，Kokkos atomic kernel 按 locus 聚合 ACGT，输出空 loci、跳过 reference-N
并按 GATK 规则计算 `REF_COUNT`/`ALT_COUNT`/`ALT_NUCLEOTIDE` 的 TSV，接入 CTest、
dispatcher contract 和 `benchmark_collect_allelic_counts.py`；完整 HDF5/cloud
container、sample discovery 和 warning/report parity 仍显式 fallback。

`ModelSegments` 已将该链路闭合到 `modelFinal.segments.tsv`：Host 严格解析
`DenoiseReadCounts` 的连续区间和 `CollectAllelicCounts` 的 `REF_COUNT/ALT_COUNT`，
将 minor-allele fraction 与 copy-ratio 一起放入选定 execution space 的 Kokkos
`RangePolicy` change-point kernel；等位基因失衡可以独立触发边界。异质位点过滤、边界、
embedding、bounded posterior 和多样本阶段统一走
`HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect`，并记录 kernel
batch/observation、execution space 与 prepare/execute 时间。Host 仅做坐标/字符串连续性、
running-mean 保守细化和确定性输出，结果保留 `MEAN_LOG2_COPY_RATIO` 兼容列，
并追加 GATK-shaped `NUM_POINTS_ALLELE_FRACTION` 及 copy-ratio/MAF posterior 10/50/90
近似列，manifest 记录 signal、阈值和 policy。显式 kernel 参数现在选择确定性的
bounded `KernelSegmenter` 路径：固定 seed 1216 anchor、linear/Gaussian CR/AF kernel、
repeatable window、changepoint cap/penalty，并以 O(N×approximation-dimension) 内存运行，
Serial/OpenMP 输出一致；默认无 kernel 参数时保持原阈值路径。无有效 allelic count 时保持
copy-ratio-only 行为。`--mode`、`--num-samples`、`--num-burn-in-iterations` 现在由
Kokkos `RangePolicy` 执行计数器驱动的 bounded deterministic random-walk Metropolis；
进入 bounded chain 前，native 按 GATK `CopyRatioSegmentedData` 的 midpoint 归属构造
原始点 Views，执行固定轮数的 Gibbs-conditional responsibility prepass：交替更新点级
outlier responsibility、每段条件均值、全局 variance 与 Beta(5,95) outlier probability。
固定顺序 Host reduction 只负责全局标量，避免 Serial/OpenMP reduction 顺序漂移；条件
均值与 variance/outlier probability 传入 bounded chain 的初始化和步长。Kokkos 只保留
固定 block，Host 维护每段固定容量 reservoir 后负责 quantile 汇总，manifest/contract
记录 sampler、conditional 模型、迭代次数、global 参数及 block/reservoir 上限。该路径保持
Serial/OpenMP 可重复，但不承诺 Java release-specific slice-sampler draws/report raw-bit
identity；完整多样本模型以及 HDF5/PoN/GC 状态仍显式 fallback。pinned
`fastgatk-model-segments-copy-ratio-conditionals-gatk-oracle` 比较 GATK 4.6.2.0 固定
partition/report schema，并以含离群点 fixture 验证条件模型没有被 outlier 拖离主模态。
相同的 count-backed allele-fraction 证据现在按 GATK `AlleleFractionSegmentedData`
形状进入 Kokkos Views，执行 binomial responsibility conditionals，更新每段 minor
fraction 以及 mean-bias/bias-variance/outlier 全局摘要，并用于 MAF chain 初始化/步长；
同一 pinned oracle 还验证 MAF posterior deciles、`.af.param` schema 和
allele-fraction conditional manifest 字段，CTest 以独立的
`fastgatk-model-segments-allele-fraction-conditionals-gatk-oracle` 名称暴露该门禁。契约见
`fastgatk-native/scripts/verify_model_segments.py`，文件边界 benchmark 见
`fastgatk-native/scripts/benchmark_model_segments.py`。

输出还会在提供 `--output-prefix` 时 materialize `.modelBegin.seg`、`.modelFinal.seg`、
`.cr.seg`、`.cr.igv.seg` 和 `.af.igv.seg` sidecar，并将路径写入 manifest，供
Picard/IGV/Nextflow 下游直接消费。
补充的 GATK advanced `--segments` Picard interval-list 现在真正作为固定 modeling
partition 使用：按 copy-ratio interval midpoint 归属、校验 sequence-dictionary 顺序
与非重叠，跳过 kernel segmentation，并在 manifest 标记
`segmentation_method=provided`。同时 materialize GATK-shaped
`.modelBegin/.modelFinal.cr.param/.af.param` 参数表（九个 posterior decile 列），
使 Nextflow/Picard 依赖的文件表面完整；参数值仍是 native bounded model 的保守摘要，
不宣称 Java MCMC raw-bit 等价。标准 copy-ratio/allele-fraction MCMC 与
smoothing 参数名会被解析并写入 manifest；显式提供 smoothing 控制时，native 已执行
一轮或多轮 adjacent credible-interval 合并和最终确定性 refit，并标记
`smoothing_applied=true`；未显式请求时保持历史 native contract，Java 多轮中间 refit
仍是后续 bounded fallback。
重复提供 `--denoised-copy-ratios` 与 `--allelic-counts` 时，native 现在按共同
interval key 建立 bounded joint evidence table，使用同一 Kokkos segmentation path
输出 Picard `.interval_list`；manifest 记录 sample/common-point 数量。完整 Java
多样本 posterior/modeling 仍显式 fallback。

`ModelSegments` 的 allelic 输入现在先经过 GATK
`NaiveHeterozygousPileupGenotypingUtils` 形状的过滤：支持
`--minimum-total-allele-count-case`、`--minimum-total-allele-count-normal`、
`--genotyping-homozygous-log-ratio-threshold` 和
`--genotyping-base-error-rate`，以 Beta 区间积分的 homozygous/heterozygous
log-ratio 判定 heterozygous loci。matched-normal 模式只保留 normal 判定为
heterozygous 的 case loci；case-only 多样本模式取所有 case 的 locus 交集；
allelic position 采用 interval overlap 而不是只匹配 bin 起点。single-sample
模式会输出 `.hets.tsv`，matched-normal 还输出 `.hets.normal.tsv`，并在 manifest
记录过滤前后 loci、阈值和 sample mode。该增量仍不等价于 Java 的完整 MCMC
modeling/parameter report，但不再把原始 homozygous pileup 静默送入分段。

`GatherTranches` 已补 VQSR scatter/gather 的 Host report 边界：输入支持 plain 或
`.tranches.gz`（zlib Host 解码），按
`minVQSLod` 合并 shard、重算加权 Ti/Tv 与 truth-site 计数，按 target truth
sensitivity 选择 tranche，写出 GATK-shaped version-5 CSV；VariantRecalibrator
的完整 Java raw-bit VBEM/covariance convergence 和 ApplyVQSR allele-level annotation 仍显式 fallback；
bounded diagonal-GMM scoring 已接入 VariantRecalibrator。输出使用 `.gz` 后缀时
也通过 zlib 流式写出，并在 manifest 记录压缩状态；解压后的 CSV 保持相同字节。
现在显式要求 GATK 同样要求的 `--mode`，兼容 `-tranche` 短别名，并将只有 version-6
表头的 scatter 空切片视作零行输入；pinned GATK 4.6.2.0 oracle 已逐字节覆盖该 empty-slice
merge。契约和
benchmark 分别为 `verify_gather_tranches.py` 与 `benchmark_gather_tranches.py`。

`VariantRecalibrator` 的 GATKReport v1.1 模型可在 Host 侧以 plain text 或
`.model.gz` 交换；压缩模型由 HTSlib/zlib 解码/编码后复用同一 Kokkos
VBEM/GMM scoring path，压缩只改变容器，不改变表格字节内容。

`ApplyVQSR` 已补上对应的 VCF native boundary：HTSlib Host 读取带 INFO/VQSLOD
的 recal VCF；scalar recal 支持 GATK dummy `N/<VQSR>` allele，按 `END`/坐标匹配，
AS recal 则按 contig/position/REF/ALT join。`--truth-sensitivity-filter-level` 选出的
tranches（plain 或 `.tranches.gz`，由 HTSlib Host 解码）按 GATK 反向 tranche walk 生成区间 FILTER、最低 tranche `+` 后缀和 PASS
边界；显式 cutoff 仍使用 `LOW_VQSLOD`。AS 模式按 ALT 长度识别 SNP/INDEL，支持
mixed multiallelic 的 `NA` 占位，per-ALT `AS_FilterStatus`/`AS_culprit` 按各自 score
分类，site FILTER 使用当前 mode 的最高 VQSLOD；scalar 输出按 GATK 写入
小写 `culprit`（同时接受旧的 `CULPRIT` 输入），AS 保持 `AS_culprit`；两条路径同步传播
`POSITIVE_TRAIN_SITE`、`NEGATIVE_TRAIN_SITE`。预先带 FILTER 的输入默认不重算，
`--ignore-filter`/`--ignore-all-filters` 显式覆盖；缺失参与 ALT 的 score 或 AS 字段
返回 `BACKEND_UNAVAILABLE`。BGZF/Tabix、OutputManifest、CTest、dispatcher contract
和 `benchmark_apply_vqsr.py` 已覆盖这些边界；OutputManifest 记录 recal/tranches
文件 provenance（size/mtime）、cutoff 来源以及 pre-filtered/ignored 计数；跨两次
ApplyVQSR 的联合 allele 状态、完整模型 provenance 和更复杂 recalibration-table 仍显式 fallback。
ApplyVQSR 还支持重复 `-L/--intervals`/`--region` 与
`--interval-set-rule UNION|INTERSECTION`；集合在 mode/scoring 前应用，空交集
fail-closed，summary/manifest 保存规则、输入统计和跳过记录；双 backend contract
与既有 GATK scalar/AS oracle 已通过。
`--create-output-variant-index` 也按 GATK optional boolean 解析（默认 true）；显式
`false`/`0` 时不创建 `.tbi`，并将该选择写入 summary/OutputManifest。
当前 scalar 3-record fixture 已与 GATK 4.6.2.0 逐记录比较（数值 INFO 格式归一化），
包括 GATK 小写 `culprit` 和 tranche `+` 边界；当前 AS oracle 再比较 2 条单 ALT
记录的 score/status；完整 header/浮点格式字节一致性和复杂 AS provenance 仍不宣称。
AS bounded slice 随后加入真实 multi-ALT 形状：GATK 4.6.2.0 的 one-ALT recal
记录以 scalar `culprit` 表示每个 ALT 的原因，native 按 contig/POS/REF/ALT join
并投影到输入的 Number=A `AS_VQSLOD`/`AS_FilterStatus`/`AS_culprit`；pinned
`verify_apply_vqsr_gatk_oracle.py` 与 `fastgatk-apply-vqsr-gatk-oracle` CTest
在 OpenMP/Serial 均通过。完整 VQSR model、recalibration-table 和跨运行 AS
provenance 仍显式 fallback。

`VariantRecalibrator` 已补上 VQSR scoring boundary：Host 读取 GATK-labelled
training/truth/known resource VCF，抽取数值 INFO annotations，在 Host 计算确定性的
training-only annotation normalization 并序列化真实 AnnotationMeans/Stdevs；resource
`prior` 按 GATK Phred prior 转换为 log10-odds，在最终 contrastive VQSLOD 边界加入，
并由 OutputManifest 记录范围；缺失 annotation 保留为 scoreable datum，training
归一化按每一维非缺失值统计，Kokkos likelihood 通过 missing mask 跳过缺失维度，
并由 OutputManifest 记录缺失记录数；对于缺失维度，Host 现在按 GATK 固定
随机种子执行缺失归一化，并在 good/bad 模型 pass 顺序下做 20 次高斯抽样边缘化；
完整记录继续由 Kokkos 执行
diagonal Gaussian log-likelihood ratio；`--max-gaussians > 1` 启用确定性 bounded GMM，
`--max-negative-gaussians` 控制负模型上限；K-means warm start 与
`--max-iterations`/`--k-means-iterations` 固定初始化，`--shrinkage`、`--dirichlet`、
`--prior-counts` 和 `--standard-deviation-threshold` 启用 Java-aligned 正则化 VBEM 更新与 outlier
gate；native 默认使用带正则化逆矩阵的 full-covariance GMM，`--full-covariance` 仅保留为兼容旧调用的拼写；VBEM 结束后按 GATK `evaluateFinalModelParameters` 重新计算责任加权最终模型；`--AS` 对 `Number=A` annotation 按 ALT 独立输出 `AS_VQSLOD`，scalar
annotation 直接 fail-closed。EM E-step 使用 Kokkos `RangePolicy`，M-step 使用
`TeamPolicy` team-local weighted tables，再由 Host 按 league 顺序确定性合并；训练
与评分 policy、`variational_normal_wishart` 都写入 OutputManifest；E-step、M-step
和 scoring 统一走 `HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect`
生命周期，并记录 policy、batch/observation、prepare/execute telemetry。`--output-model` 写出 GATKReport v1.1 风格的
AnnotationMeans/Stdevs、PMix、means/covariance 表，`--input-model` 按 annotation
name 校验并保留 serialized model 原有维度顺序后复用，可用于 scatter/reuse。工具输出带
`VQSLOD`/`CULPRIT`（AS 模式另写 `AS_VQSLOD`/`AS_culprit`）的 recal VCF、tranches、可选
R sidecar 和 OutputManifest；压缩 recal VCF 使用 .tbi，普通 recal VCF 复用 shared
Tribble LinearIndex v3 .idx，索引路径和 enabled 状态写入 manifest。完整
 Java raw-bit convergence/resource calibration、recalibration-table provenance
和完整 AS model semantics 仍显式 fallback；Host decode 已复现
`decodeAnnotation(..., jitter=true)` 的零点/端点抖动与 `--mq-cap` logit 变换，
使用固定 GATK seed 并在 manifest 记录 jitter controls。契约和 benchmark 为
`verify_variant_recalibrator.py` 与 `benchmark_variant_recalibrator.py`，并已纳入
`verify_all.sh`。

当前还补齐了一个可独立验收的 GATK model-artifact 边界：Host 可读取
GATK 4.6.2.0 的 GATKReport v1.1 `AnnotationMeans`/`AnnotationStdevs`、
`Positive/Negative Gaussian PMix`、means 与 covariance 表，按 annotation 名称
重排后复用 Kokkos scoring；pinned Java oracle 比较 VQSLOD、坐标/记录数、稳定的
正资源标签与模型表维度；GATK 的模型侧负位点选择不写入 GATKReport，因此作为
非 bit-identical metadata 边界记录。scalar recal VCF writer 对齐 GATK 的 `N/<VQSR>`、`END`、
小写 `culprit`、空 QUAL/FILTER 以及 `--sites-only-vcf-output` 的 8-column
形状。OpenMP/Serial 的 model oracle CTest 均通过；Java raw-bit VBEM 收敛、
完整 recalibration-table/provenance、完整 AS model parity 仍显式 fallback，
不将该 bounded slice 宣称为完整 VQSR。

### 4.12 Annotation / Funcotator / variant utility（P2/P3）

**C++ Host**

- variant normalization、left-align、header/schema、data-source interval lookup、
  VCF/MAF/SEG renderer、字符串和插件生命周期全部 Host；优先 memory-map/index
  cache。

**Kokkos kernel**

- 只有纯数值、固定字段的 annotation（depth、allele count、quality summary）
  可批量 `parallel_for/reduce`；字符串/区间查询不送 device。

**验收**

- VCF/MAF/SEG 字段顺序、缺失值、header、renderer 和 GATK plugin 结果；性能门槛
  必须证明 lookup/cache 才是瓶颈后才移植。

当前已落地 `VariantsToTable` 的 HTSlib streaming 子集：site/INFO、FORMAT、
allele-specific 和 allele-specific-genotype 字段，稳定 sample 排序、`NA`
缺失值、默认过滤记录跳过、`-SMA`、`--moltenize`、重复 `-L` 区间筛选（按记录
span/`END` overlap）以及 OutputManifest/CTest contract；标准 VariantContext getter（HET/HOM-REF/
HOM-VAR/NO-CALL/VAR/NSAMPLES/NCALLED/SAMPLE_NAME/TYPE/EVENTLENGTH/
TRANSITION）、Java allele-valued GT、QUAL/Filter 文本、INFO wildcard、
`-XL`/`-ip`/`-ixp` 区间边界和 ASGF A/R/G split
规则已由 pinned GATK 4.6.2.0 table oracle 逐字节回归，并有 file-boundary
benchmark；QUAL 缺失值（`-10.0`）、`*` spanning-deletion type、非标准字段的
`NA` 与 bare INFO/FORMAT key 语义也已加入 oracle 回归；文本 VCF INFO token
在 Host 边界保留长小数、指数拼写、`.` 和 flag，typed FORMAT/QUAL 使用共享
Java-style shortest serializer；typed BCF、复杂 annotation/plugin 及极端 Java
数值格式仍保持 explicit fallback。
新增 release-pinned `verify_variants_to_table_gatk_oracle.py` 固定了更深的
row-renderer 边界：`-SMA` 的 biallelic 记录也走 allele-specific 投影，`-ASF`
Number=R 先删除 REF 后按 GATK 的 List-string 格式保留非首项空格，`-EMD` 按字段
是否存在而不是渲染后的 `NA` token 判断缺失。`--moltenize` 还刻意复现 GATK 4.6.2.0
的顺序消费行为：只用 `-F/-GF` 标签写行，却在含 `-ASF/-ASGF` 时从完整提取列表
连续取值；这不是理想化接口，但为直接替换保留真实 Java 输出。
兼容参数面还接受 `--interval-merging-rule ALL|OVERLAPPING_ONLY` 和
`--variant-output-filtering STARTS_IN|ENDS_IN|OVERLAPS|CONTAINED|ANYWHERE`；前者
控制半开区间 selector 的相邻合并并写入 manifest，后者按 GATK
`VariantsToTable` 的实际行为做参数校验后仍采用普通 OVERLAPS table traversal
（该工具不安装 VCF writer 的输出过滤装饰器），相关模式已加入 native/GATK
命令级回归。

`VariantEval` 已补齐同一层的 Host/HTSlib evaluator 子集：`CountVariants`、
`TiTvVariantEvaluator`、`VariantAFEvaluator`、`ThetaVariantEvaluator`、`MendelianViolationEvaluator`、`VariantSummary`、`IndelSummary`、`MultiallelicSummary`、`IndelLengthHistogram`、显式
`GenotypeFilterSummary`、`PrintMissingComp`、`CompOverlap`、
`GenotypeConcordance`，并实际输出 Contig/Filter/VariantType/AlleleFrequency/Novelty/Sample
stratification；`IndelSummary` 按 concrete ALT 统计 insertion/deletion、large-indel
和 multi-allelic site，`MultiallelicSummary` 统计 SNP/INDEL site 比例及 comparison
partial/complete novelty，`IndelLengthHistogram` 输出 -10..-1/1..10 的 GATK molten
频率表。支持 repeatable `-eval`/`-comp`、`-EV`/`-ST`/`-S`，`-L` 对 gVCF 使用记录
span/`END` overlap，Sample stratification
会按真实 sample GT 输出 `CountVariants` rows（并以 GATK 4.6.2.0 oracle 校验），共享区间解析、
确定性 report/OutputManifest 和 dispatcher/CTest contract；重复 `-L/--intervals`
支持 GATK `UNION`（默认）或 `INTERSECTION`（`-isr/--interval-set-rule`），且
eval/comparison track 共享同一集合结果。显式 `--gatk-report` 已输出
GATKReport v1.1-compatible core tables，`benchmark_variant_eval.py` 已接入 `verify_all`；
完整 evaluator/stratifier plugin graph、全部 stratified GATKReport exact formatting
及远程 feature 输入仍保持 explicit fallback。

`Novelty` 的 bounded path 使用 comparison track 的坐标/等位基因匹配生成
`all`/`known`/`novel` CountVariants rows；`--no-st --stratification-module Novelty`
的 GATKReport 和 known/novel telemetry 已加入 `verify_variant_eval.py`，但不把
该增量扩大解释为完整 Java stratifier/plugin graph。

VariantEval 的下一步 bounded oracle 已收敛到可逐字节承诺的报告切片：
`CountVariants` aggregate 与 `Sample` rows 现在采用 GATKReport v1.1 的真实
`%d`/`%.8f`/`%.2e`/`%.2f` metadata、Java 数值格式和列宽对齐；GATK 的
`nProcessedLoci`、reciprocal rate、het/hom ratio、singleton 与
`nHomDerived` 字段也按 pinned GATK 4.6.2.0 修正。独立的
`verify_variant_eval_gatk_oracle.py` 同时比对 aggregate/Sample report bytes，
并覆盖 genotype-aware ignore-AC0 计数：含 FORMAT/GT 的 hom-ref/all-no-call ALT
记录计入 reference loci，sites-only 记录保持 site-level variant，同时固定
singleton 与 insertion/deletion ratio 边界。
并用截断 VCF 固定 HTSlib 读错误必须 fail-closed。该 oracle 只覆盖
`CountVariants` + `Sample` 的 bounded surface；其余 evaluator/stratifier 的完整
Java plugin graph、远程 feature 输入和所有报告表的 raw-byte parity 仍是
explicit fallback。

`VariantEval --keep-ac0`（及短别名 `-keep-ac0`）补齐了 GATK 的 optional-Boolean
入口，并在 native CountVariants aggregate 路径中恢复 genotyped AC=0 的 SNP/INDEL
site-level 计数。pinned GATK 4.6.2.0 oracle 对 aggregate report 做逐字节比较，
同时以 native-only checks 固定 bare/true/false 与短别名等价；该增量不扩大为其它
evaluator/stratifier 表的 Java parity 承诺。

`VariantEval -no-ev -EV ValidationReport -no-st` 已补齐 bounded evaluator：
Host 按 eval sample 集合对子代 comparison genotype 做子集，并按 GT 或 INFO/AC
区分 MONO/POLY；Kokkos 对 comparison-only locus 的 4x4 site-status matrix 做
确定性整数归约。pinned GATK 4.6.2.0 oracle 在 genotyped 与 sites-only/AC 两组
fixture 上逐字节比较完整 GATKReport，并固定默认无 Filter stratifier 时 filtered
record 在 evaluator 前被排除。多 eval/comp track、重复 eval locus 和
Filter-stratified ValidationReport 均 fail-closed，未扩张为完整 plugin graph。

`ValidateVariants` 已补齐 Host/HTSlib/FAIDX validation 子集：REF/reference
一致性、ALT 唯一性/非空、AN/GT chromosome-count、可选 dbSNP ID、重复 `-L`
筛选（按记录 span/`END` overlap）、`-XL` 后置排除和 `-ip/-ixp` padding、GVCF `<NON_REF>`/排序/overlap、excluded validation types、filtered skip、
warning mode、reference sequence-dictionary length/contig validation 和
INFO/FORMAT Number=A/R/G、固定长度与跨样本 GT ploidy validation、issue/telemetry OutputManifest，并接入 dispatcher/CTest；
重复 `-L/--intervals` 支持 GATK `UNION`（默认）与 `INTERSECTION`（`-isr`/`--interval-set-rule`），
并在验证前对半开区间集合做确定性运算，规则写入 summary/manifest；
`--disable-sequence-dictionary-validation` 现在显式生效。完整 Tribble/htsjdk
validation、dictionary-wide GVCF coverage 和原始异常文本仍保持 explicit fallback。
本轮补齐 GATK/Barclay 单短横线别名边界：`-D`（dbSNP）、`-isr`（interval set
rule）、`-do-not-validate-filtered-records`、`-warn-on-errors` 与
`-disable-sequence-dictionary-validation`；这些别名由 pinned GATK 4.6.2.0
command oracle 与 native 双 backend 对照，解析失败仍保持 fail-closed。
当 `--validate-gvcf` 同时提供 reference 且未指定 `-L` 时，现已额外校验每条
sequence 的 leading/trailing coverage 和 block `END` 不越界；更复杂的
GVCF annotation/posterior 语义仍保持 explicit fallback。`verify_validate_variants.py`
覆盖 `-XL` 后置排除与 `-ip/-ixp` padding，并与 GATK `IndexFeatureFile` +
`ValidateVariants` 做命令级 oracle；`benchmark_validate_variants.py` 已接入
`verify_all`，输出 records/s 与报告/manifest 大小。

`GetPileupSummaries` 已补齐 contamination workflow 的 native 子集：共享
HTSlib Host reader 读取 SAM/BAM/CRAM，按 AF 范围筛选 biallelic SNP，支持 literal/
Picard/BED 和 VCF-site `-L` selector，复用 GATK 默认 mapping/duplicate/primary/
vendor/mate/CIGAR/Wellformed read-filter mask；支持的内置 `--read-filter` 与
`--disable-read-filter` 名称会显式校验并改变过滤配置，`--disable-tool-default-read-filters`
从空过滤集开始。ref/ALT/other base 计数统一走 Kokkos kernel，输出 GATK-shaped
`SAMPLE` metadata + six-column table，并接入 OutputManifest、CTest 和 dispatcher
contract。重复 `-L/--intervals` 支持 GATK `UNION`（默认）与 `INTERSECTION`
（`-isr/--interval-set-rule`），先在 Host 侧归一化后再驱动 population-site
indexed traversal 与 Kokkos count kernel，规则写入 summary/manifest。完整 LocusWalker pileup/filter、cloud/dictionary 和 Java report formatting
仍保持 explicit fallback。
计数阶段现在也遵循公共 `HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect`
生命周期，manifest 记录 count kernel 的 execution space、batch/observation 数量以及
prepare/execute 时间。`benchmark_pileup_validation.py` 将 Get/Gather 与
ValidateVariants 的文件边界放在同一套 warmup/p50/p95 harness 中，并接入
`verify_all`；benchmark 只用于可重复性能基线，不替代 Java/GATK oracle。

`ReblockGVCF` 与 `SelectVariants` 的 Number=G allele-subset PL 重映射现在共享
`fastgatk-kernels` 的 Kokkos API，非参考 AD 清理也统一走
`HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect`/`RangePolicy`
生命周期并记录 prepare/execute telemetry；分别在 `benchmark_reblock_gvcf.py` 和
`benchmark_select_variants.py` 中报告 PL remap 与 GT/GQ prepare/execute telemetry；
固定非二倍体 ploidy 由 triploid regression 覆盖。ReblockGVCF native prototype
现在还按 sample-major 批次处理同一 header 的多样本 GVCF：AD/PL remap、GT/GQ、
逐样本 MIN_DP 以及 RAW_GT_COUNT 均保留 sample 行，reference block 合并要求每个
样本的 GQ band 一致；但 GATK 4.6.2.0 拒绝多样本 Reblock 输入，该扩展不计入直接
替换能力，跨样本 joint posterior/annotation 仍显式 fallback。

HaplotypeCaller 现在接受 `--sample-ploidy 1..8`，按 VCF Number=G 的非递减
genotype rank 由 per-read PairHMM likelihood 生成任意倍体 VCF/GVCF
PL/GT/GQ，并按实际倍体计算 AC/AN；共享 Kokkos RCM primitive 会保留每条
observation 的 REF/`<NON_REF>` likelihood，再生成 ploidy-specific reference block。
`verify_hc_ploidy.py` 和 `benchmark_hc_ploidy.py` 覆盖 haploid/triploid/octoploid；
无 reference 的 polyploid posterior 仍显式 fail-closed。
多等位 normal VCF 的 site QUAL 不再取独立 ALT QUAL 的最大值：writer 将完整
Number=G PL、SNP/INDEL Dirichlet pseudocount 和实际 ploidy 交给共享
`calculate_allele_frequency_kokkos`，按 GATK AFCalculator 的 flat-seed EM（阈值
0.1）计算 `P(non-ref)`；INFO/DP 同时只计入 PairHMM 保留的合格 evidence。复杂
tetra-ploid fixture 已锁定完整 candidate row（包括 QUAL、INFO、RankSum 与
FORMAT）的 GATK 文本一致性；更广泛的 annotation-engine provenance 仍是独立的
扩展边界；与 GATK 一致，
arbitrary-ploidy VCF/GVCF candidate 不伪造 diploid-only `ExcessHet` 字段。
benchmark 的 GVCF cases 读取 OutputManifest，额外报告 RCM execution space、
标准字段契约和 final-genotype block rebuild 状态。

HC `-ERC BP_RESOLUTION` 复用同一 Kokkos reference-confidence kernel，关闭
GQ-band block coalescing，逐碱基写出 REF/`<NON_REF>` 记录并保持
`GT:AD:DP:GQ:PL`/`INFO=.` contract；普通 `-ERC GVCF` 仍走 block 模式。HTSlib
从 `@RG SM` 收集 sample identity，默认使用首个唯一 sample；新增
`--sample-name` 在 Host decode 边界按 `@RG ID→SM` 过滤并将选中样本置于唯一
VCF/GVCF sample 列，未知样本 fail-closed，缺失 SM 才回退为 `FASTGATK`。
`verify_native.py`/`verify_indel.py` 已覆盖两种 ERC 模式、sample header 和多样本
过滤 contract。HC/HtsReader 还支持 `-XL/--exclude-intervals` 及 `-ip/-ixp`
padding：include/exclude 在 indexed 与 sequential 路径统一归一化，按 projected
reference locus 应用，跨边界 read 不会被整条丢弃；候选、PairHMM 和 GVCF block
writer 共用同一 exclusion mask，并由 `verify_hc_region_streaming.py` 的 GATK
fixture 回归覆盖。
重复 `-I/--input` 也已接入 Host 聚合和流式路径，按输入顺序合并多个 BAM/CRAM
shard，并校验序列字典及（未指定 `--sample-name` 时）样本集合一致；contig 流式
使用有界坐标 k-way merge，region 流式在每个 indexed halo tile 合并各 shard，避免
跨文件坐标重入并保持有界内存语义。
重复 `-L/--intervals` 还支持 GATK `-isr/--interval-set-rule` 的 UNION（默认）和
INTERSECTION，规则在共享 HtsReader 中统一应用并写入 telemetry/manifest。
提供 `-R` 且存在可读 `.fai` 时，HC Host 会在遍历前按 contig 顺序和长度校验
FASTA/FAI 与 BAM 字典；`--disable-sequence-dictionary-validation` 是显式可选布尔
逃生开关，跳过该检查并写入 summary/manifest。无索引的 plain FASTA 仅作为 smoke
模式保留，不伪称具有权威 sequence dictionary。
普通 GVCF 的 blocks 现在在最终 genotype selection 后重建：未成 call 的 SNP/MNP
候选回收到 `<NON_REF>` block，未成 call 的 indel 仍保留其完整 REF span，避免
candidate-site 与 reference block 重叠；零支持 assembly-only indel 按
`--min-alt-support` 回收到 reference block。普通 block 输出采用 GATK-shaped
`GT:DP:GQ:MIN_DP:PL`、`FILTER=.`、`INFO=END`，candidate 输出采用
`GT:AD:DP:GQ:PL:SB`；candidate 的 QUAL/FILTER、INFO 键序（含
`MLEAC/MLEAF/RAW_MQandDP`）及未封顶 Number=G PL 已与 GATK 对齐。manifest 明确记录
`gvcf_standard_fields`、`gvcf_candidate_standard_fields` 与
`gvcf_candidate_rebuild_after_genotyping`。bundled chr17 oracle 已验证 candidate
site 集合相等及 15 条 GVCF record 的结构契约；RCM 数值/最后一个 block 边界仍是
数值模型差异，bit-identical 继续保持显式未达成。

HaplotypeCaller 的 `GenotypePriorCalculator.assumingHW` prior 分类也已修正：
同长单碱基变异使用 SNP/log10(3) normalization，长度变化使用 INDEL
heterozygosity，MNP/complex/symbolic 等 OTHER 使用
`max(snpHet, indelHet)` 且不做 SNP normalization，避免把非 SNP allele 错算为
SNP prior。

`ReblockGVCF` 的 high-quality variant 现在补齐 GATK `RAW_GT_COUNT` INFO
（固定 hom-ref/het/hom-var 三元组），从最终 GT 计算并记录到 manifest telemetry。
`--do-qual-approx` 同时写出 variant-only AD 规则下的 `VarDP`/raw `AS_VarDP`，
并记录注释计数；完整 annotation-engine 聚合仍保持 fallback。
low-quality deletion/span-deletion/no-call ref-block conversion 也按 GATK
规则 trim REF 到首碱基并保留 END；high-quality deletion trimming 与提供 `-R`
时的单样本 deletion-gap reference block 已实现，复杂 multi-sample/overlap
writer 仍 fallback。
FILTER/FT 参数也已落地：reference block 清除 FILTER，variant 默认清除、
`--keep-site-filters` 保留，`--add-site-filters-to-genotype` 写 FORMAT/FT。
high-quality genotype 的 `<NON_REF>` AD 现在通过 Kokkos cleanup kernel 清零，
并按 GATK `removeNonRefADs` 从 FORMAT/DP 扣除；计数写入 OutputManifest telemetry。
stale fixed/dynamic `GVCFBlock*` INFO annotation 的 record/header 清理以及
`--annotations-to-keep`/`--annotations-to-remove` 参数也已接通，删除/保留计数写入
OutputManifest；完整 VariantAnnotatorEngine 聚合仍显式 fallback。
在 high-quality ALT compaction 后，native 现在按 GATK
`trimAlleles(..., false, true)` 规则对非 symbolic alleles 做共同后缀裁剪，
同步更新 variant REF span，并以 regression/manifest telemetry 覆盖；完整
annotation/posterior 聚合仍保持显式 fallback；提供 `-R` 时，单样本
`addRefBlockIfNecessary` deletion-gap reference block 也使用 FAIDX 参考碱基与
共享 Kokkos PL 投影生成，并记录 block/fallback telemetry，多样本与 overlap-aware
writer 语义仍保持 explicit fallback。

`CalculateContamination` 已补齐 Java `ContaminationModel` 的三轮 MAF/contamination
学习、matched-normal/tumor-only strategy cascade、四 genotype posterior、精确误差
二分搜索和 GATK-shaped `sample/contamination/error`/MAF 输出；每个 segment 的
likelihood 统一走 `HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect`/
`RangePolicy` 批量执行，manifest 记录 execution space、kernel batch/observation 数及
prepare/execute 时间；LikelihoodPlan 的输入/临时 Views 在多轮评估与可选 segmentation
重学习间持久复用，`likelihood_plan_allocations/reuses` 写入 manifest/benchmark。plain
或 `.gz` table 输入，以及 `.gz` 输出/segmentation，均由 Host zlib 流式边界处理。
四个真实 GATK 4.6.2.0 fixture
已接入 oracle（contamination delta <0.005、error delta <0.001）。当前 changepoint
已使用 seeded kernel approximation、deterministic Commons-Math/JAMA-style bidiagonal SVD、persistence
和 backward selection，窗口/惩罚参数与 Java `KernelSegmenter` 一致；对称核矩阵与
reduced-observation 投影也统一走 `HostBatch -> KernelPlan.prepare -> Kokkos Views
-> execute -> collect`/`MDRangePolicy`，并在 manifest/benchmark 中单独记录
segmenter kernel telemetry；尚未完成
Commons Math 3.5 SVD 的逐元素数值 parity，cloud 和 bit-format report 继续
explicit fallback。

`GatherPileupSummaries` 已补齐 scatter/gather 子集：读取多个 GATK-shaped
PileupSummary table，按 `.dict`/`.fai`/FASTA contig 顺序稳定排序，校验 sample
metadata，忽略空 shard，默认保留重复 locus（与 GATK 一致），
`--reject-overlaps` 时显式拒绝，
并在存在非空 shard 时输出原工具可读的 SAMPLE metadata/table 与 manifest；all-empty
输入则按 GATK 只写列名 header。现在也解析 `.dict`/`.fai`/
FASTA 的 contig length，拒绝越界位置和未知 contig；
`--disable-sequence-dictionary-validation` 时未知 contig 按稳定字典外顺序追加。
plain 与 `.gz` table 输入/输出均在同一 Host 文件边界完成，压缩输入数量、输出压缩
状态和完整性写入 manifest；AF 使用 binary64 shortest-round-trip 格式化，已由高精度
GATK 4.6.2.0 oracle 验证；cloud URI 与完整 Java dictionary edge semantics 仍显式
fallback。
样本名直接取 reads header 首个 `@RG SM`；当 header 没有任何样本时与 GATK 一样
fail-closed，不生成 `SAMPLE=UNKNOWN` 的不可消费表。

`GetPileupSummaries` 的 count kernel 也采用持久化 observation/ref/ALT/count Views，
跨 read batch 复用容量并把 allocation/reuse telemetry 写入 manifest，保持 Kokkos
execution-space 切换不改变 Host 过滤与输出顺序。
完整 cloud 和 Java sequence-dictionary edge semantics 仍保持 explicit fallback。

### 4.13 SV / BWA / PathSeq / Spark / Picard 长尾（P3，fallback-first）

**C++ Host**

- 第一阶段 registry 标为 `fallback(gatk|picard|spark)`，保持 Nextflow/SLURM 命令
  和输出，不在 native runtime 内复刻 Spark/Hadoop。
- 只对 profile 证明收益的孤立 primitive 建 native adapter，例如 BWA seed/extend、
  SV evidence count、PathSeq k-mer classify。

**Kokkos kernel**

- 每个孤立 primitive 先以 flat batch + Kokkos correctness kernel 验证，禁止直接
  把完整 Spark DAG 或字符串/图控制流搬进 device。

**验收**

- fallback 透明、退出码/日志可解释；native 仅在 tool-specific oracle 和资源
  benchmark 通过后进入 registry。

## 5. 工具级组合方案

### 5.1 HaplotypeCaller 主路径

```text
compat -> HTSlib reader -> TraversalPlan(core+halo)
       -> activity/pileup/downsample (Host + Kokkos primitive)
       -> local assembly (Host graph + Kokkos k-mer)
       -> Smith-Waterman
       -> PairHMM
       -> genotype/reference confidence
       -> deterministic VCF/GVCF writer + index
```

Host 保持 assembly region 生命周期和输出；Kokkos Plan 在 region/batch 内复用，
避免 PairHMM 已观察到的 prepare 串行瓶颈。HC 首个可替换版本应只接管固定参数子集，
未知参数和不支持模式 fallback。

### 5.2 Mutect2 主路径

```text
同一 traversal/assembly/PairHMM/SW
       -> tumor/normal likelihood
       -> somatic filters + contamination/orientation (Host)
       -> VCF + .stats + F1R2 + somatic GVCF
```

先做到 sidecar/output 完整，再优化 somatic 数值；不能只输出 VCF。

### 5.3 BQSR 主路径

```text
HTSlib read stream -> Host covariate descriptor
                  -> Kokkos local table batches
                  -> deterministic integer merge
                  -> Host recalibration report/table
                  -> Host ApplyBQSR stream
```

BQSR 的关键不是单 kernel 峰值，而是 table merge、report schema、重启和输出 read
顺序。因此 batch plan 必须能在 shard 结束时序列化。

### 5.4 GenotypeGVCFs/GenomicsDB

```text
GenomicsDB adapter (Host, budget-aware)
       -> locus blocks
       -> Kokkos likelihood/genotype math
       -> Host annotation + GVCF block/VCF writer
```

数据库访问和文件句柄限制由 Host 控制；设备只处理连续的小矩阵。

## 6. 集群与兼容汇总方案

### 6.1 Nextflow/SLURM 不变

现有 process 继续调用：

```bash
gatk HaplotypeCaller ...
```

launcher 在作业内读取 `SLURM_CPUS_PER_TASK`、`SLURM_MEM_PER_NODE`、`SLURM_TMPDIR`
和 GPU GRES；不自行 `sbatch`，不突破 cgroup，不把跨节点细粒度任务放到 Kokkos。
scatter 仍以 sample/interval/shard 为单位，native 输出可被原 gather/resume/cache
使用。

### 6.2 产物矩阵

```text
fast-gatk-cpu-x86   Kokkos OpenMP + scalar/AVX2/AVX-512 builds
fast-gatk-cpu-arm   Kokkos OpenMP/SVE-safe build
fast-gatk-cuda      Kokkos CUDA + host OpenMP
fast-gatk-hip       Kokkos HIP + host OpenMP
fast-gatk-sycl      Kokkos SYCL + host backend
```

每个产物共享 tool registry、schema、golden corpus 和版本。运行时选择设备；
设备初始化失败可回退 CPU native，再根据工具策略回退 Java。

### 6.3 输出和错误契约

所有 native tool 关闭 writer 前执行：文件存在、索引可读、header/contig 一致、
core interval 不重复、排序合法、sidecar 完整、CRC/manifest 正确。错误统一分类：

```text
BAD_INPUT
RESOURCE_EXHAUSTED
BACKEND_UNAVAILABLE
NUMERICAL_CONTRACT_FAILURE
OUTPUT_CONTRACT_FAILURE
FALLBACK_EXECUTED
```

只有主输出和全部要求的 sidecar 都通过 manifest，进程才返回成功。

## 7. 实施分期和模块依赖

### P0：可运行的兼容骨架和热点 kernel

依赖顺序：

```text
compat + runtime + HTSlib/model
              -> traversal/queues
              -> PairHMM + SW
              -> k-mer primitive + host assembly
              -> HC smoke tool
```

交付：native dispatcher、CPU OpenMP/AVX2/AVX-512、Strict oracle、BAM/VCF smoke、
HC 小数据集、fallback。

### P1：完整 HC/Mutect2/BQSR 基础路径

依赖 P0 的 output manifest、plan lifetime 和 batch schema；实现 activity/pileup、
GVCF/reference confidence、Mutect sidecars、BQSR table/ApplyBQSR。

### P2：GenotypeGVCFs、GenomicsDB、外部内存和 dense math

实现 GenomicsDB adapter、genotyping blocks、MarkDuplicates/Sort external-memory、
VQSR/CNV Kokkos Kernels；扩大真实 cohort 和 scratch/fd 压力测试。

### P3：GPU 优化和长尾

先将 PairHMM/SW/固定 k-mer primitive 迁移 Team/Vector/tile GPU，再处理 VQSR/CNV；
SV/BWA/PathSeq/Funcotator/Spark/Picard 只按 profile 进入 native registry，其余保持
fallback。

## 8. 统一验收门槛

| 层级 | 必须证明 |
|---|---|
| API/CLI | 参数、`@args`、错误码、help/version、fallback 透明 |
| Kernel | GATK Java/GKL/reference 结果；Strict raw-bit 或逐 CIGAR；Fast tolerance |
| Format | BAM/CRAM/VCF/BCF/index、header、compression、sidecar 可回读 |
| Tool | 同输入/参考/interval/seed 的 record、排序、GVCF block、metrics 一致 |
| Resource | cgroup/SLURM memory、device OOM、scratch/inode、fd、spill/retry |
| Pipeline | Nextflow scatter/gather、resume、cache、失败重试、SLURM allocation |
| Performance | decode、assembly、PairHMM、SW、genotype、compress/write/spill 分项；同时报告 kernel 与端到端 |

任何一层不通过，工具状态只能是 `fallback`，不能标记为 native-compatible。

## 9. 第一批实际编码任务

1. ✅ 把现有 `pairhmm-demo` 提炼成 `fastgatk-kernels` 的生产 PairHMM 源码，拆分 prepare/execute/collect，加入 workspace pool 和 prefix cache；demo 仅保留为客户端、oracle 和 benchmark，不再作为生产 kernel 的源码依赖。
2. 把 `kokkos-modules-demo` 的 SW、k-mer、BQSR、region primitive 改为库目标，接入统一 `HostBatch/DeviceBatch/KernelPlan`。
3. 建立 HTSlib BAM/CRAM/VCF reader/writer 和 `OutputManifest`，替换 SAM-like fixture。
4. 实现 `gatk` dispatcher、tool registry、`HaplotypeCaller` native smoke adapter，未知参数自动 fallback。
5. 用真实 GATK SW、HC、BQSR 测试资源建立 golden corpus；先 CPU Strict，再 Fast，再 GPU。
6. 在 Nextflow+SLURM 上运行 scatter/gather/resume/retry，并记录八段式性能 telemetry。

这六项完成后，才有足够证据把“架构和模式可行”升级为“第一个 GATK 工具可直接替换”。

## 10. 当前实现状态

已落地的第一步：

- `fastgatk-core` 已成为 PairHMM 和四类模块原型共同链接的 CMake library target；
  `HostBatch/DeviceBatch/KernelPlan` 已进入 prepare/execute 计时和 View accounting。
- `fastgatk-native` 已接入 HTSlib 1.22.1，真实 BAM/CRAM 记录会按 batch 解码到 flat
  Host arrays，再进入 Kokkos reduction kernel。
- `fastgatk-hc-smoke` 已验证 BAM 全文件、BAM region filter 和 CRAM + reference；
  已增加 GATK `-I/-R/-L/-O` 参数别名、VCF/VCF.GZ 容器 writer、OutputManifest/telemetry
  sidecar，以及 `gatk_compat.sh` launcher；VCF 当前明确是空的 smoke 输出，不宣称生物学调用。
- `fastgatk-hc-call` 已接通第一条真实组合链路：HTSlib flat reads → Kokkos
  high-quality pileup/CIGAR indel assembly → bounded local haplotype combinations →
  SW score/Host traceback (GATK overhang/tie policy) → persistent bucketed PairHMM likelihood marginalization →
  Kokkos genotype/GQ。它支持 GATK 参数别名并输出候选 calls、VCF/VCF.GZ、
  OutputManifest 和跨 batch/thread 稳定的整数 signature；仍不宣称完整
  HaplotypeCaller/VCF 等价，GATK pruning/traversal/prior/reference-confidence 继续
 由 oracle 和 fallback 门禁约束。固定 chr17 broad oracle 的三条普通 VCF record
  及 `INFO`/`FORMAT`/`FILTER` schema 定义已逐行匹配 GATK；command-line/date/
  contig-assembly provenance 和更广泛的 assembly corpus 仍是独立差距。
- `verify_pipeline_local.sh` 已在无集群环境验证 Nextflow process 等价路径、双区间
  scatter/gather、VCF+`.tbi` 输出、OutputManifest 读回和 SLURM 本地资源模拟；第二次
  runner 调用验证完整 shard 的 resume。`workflow/nextflow_smoke.nf` 与
  `workflow/nextflow_scatter_gather.nf` 在安装 Nextflow 的环境中会额外执行真实 DSL2
  process，后者使用 `-resume` 复用 Nextflow cache。
- `verify_gatk_oracle.py` 已将同一 chr17 区域交给 GATK 4.6.2.0 与 native prototype，
  检查共享 sentinel call、call-set overlap，并强制 manifest 明确
  `bit_identical_to_gatk=false`；当前该 fixture 的 GATK/native call-set 交集为 1/1，
  说明数据流已经可比对，但尚未达到工具等价。

尚未完成的下一层：除已落地的 `--min-pruning` read-kmer support floor 外的
GATK-equivalent haplotype pruning/traversal、官方 genotype priors
和 reference-confidence semantics、GATK dispatcher 的 native promotion、完整
HaplotypeCaller oracle 以及 Nextflow/SLURM 真实集群 end-to-end。

HC prior correction：`calling_pipeline.cpp` 已按 GATK
`GenotypePriorCalculator.assumingHW` 对 SNP heterozygosity/hom-var prior
应用 `log10(3)` normalization，并将 per-candidate values 写入 telemetry；
`verify_hc_genotype_priors.py` 直接编译/运行 bundled GATK Java class 做逐项
oracle，覆盖 diploid/polyploid 组合。完整 assembly/haplotype posterior 仍未宣称
bit-identical。

VariantEval 补充：`MendelianViolationEvaluator` 已支持 `--pedigree/-ped` 与默认 mvq=50，
输出完整 31 字段 GATKReport，并用 trio fixture 与 pinned GATK 4.6.2.0 逐字段对照；
缺少 pedigree 关系的样本不会被猜测，按 GATK 语义跳过并写入 telemetry。
`Family` stratification 也已接入：按同一 PED 关系输出 `all` 与各 family 的 CountVariants
rows，并与 GATK Family-stratification oracle 对照 rate 字段。

## 11. 统一模式的跨后端验收（2026-08-27）

后续算法统一采用 C++ Host + Kokkos kernel 的单源模式：HTSlib 解码、排序/分块、
字符串与文件格式、spill/checkpoint、错误/fallback 保持 Host；固定形状数值阶段
只通过 `KernelPlan.prepare -> Kokkos Views -> execute -> collect`，CPU SIMD 也
使用 Kokkos SIMD API。`FASTGATK_KOKKOS_BACKEND` 只选择 execution space，不产生
后端特有算法分支。

这一模式已由同一份 CTest 门禁验证：当前 OpenMP 81/81、Serial 81/81（均包含 pinned
GATK 4.6.2.0 Oracle；FilterMutectCalls contamination oracle 的 Serial 增量为 2/2）。CTest 会注入当前 `FASTGATK_NATIVE_BUILD` 和期望的
execution-space 标签，防止 portability build 误调用旁侧 OpenMP artifact。
随后 CombineGVCFs 标准参数/分段改动的受影响集合在 OpenMP 15/15、Serial 15/15
通过（含 `--call-genotypes` sample payload oracle）。因此可以把该模式复制到其它算法；但每个工具仍须完成自己的 GATK 语义 Oracle、
真实数据端到端、资源/云边界和 fallback 后，才可从 prototype 升级为
native-compatible。

2026-08-30 增量：FilterMutectCalls 省略阈值参数时现在采用 GATK
`M2FiltersArgumentCollection` 的默认 `MBQ=20`、`MMQ=30`、`MPOS=1`、
`max-median-fragment-length-difference=10000`，且 `max-n-ratio=+Infinity`；
显式负值仍可关闭对应过滤器，零 ALT-depth 时按 Java `NRatioFilter` 跳过。
`--mitochondria-mode`/`--microbial-mode` 也已接通：两者关闭 clustered、
multiallelic、fragment、haplotype 等 genomic-only filters；microbial 保留
polymerase slippage 且默认 MMQ=20，mitochondria 关闭 slippage，模式和禁用集合
写入 telemetry/manifest。
新增的默认 hard-filter verifier 在 OpenMP/Serial 均通过，随后两套完整 CTest
分别 **79/79**（OpenMP 1240.42 s，Serial 1228.45 s）通过。随后
`--max-n-ratio` 默认 `+Infinity`、零 ALT-depth 跳过语义和合法 JSON telemetry
增量在双 backend verifier 中通过，最新完整 CTest 仍为 **79/79**（OpenMP
1292.85 s，Serial 1225.15 s）。这些增量提高 CLI 兼容性和替换边界覆盖，但不
改变完整 Mutect2/FilterMutectCalls joint posterior 仍为 prototype/fallback 的判断。
随后 `--mitochondria-mode`/`--microbial-mode` 的 genomic-only filter gating、
microbial MMQ=20 默认及模式 telemetry 回归完成；最新完整 CTest 仍为
**79/79**（OpenMP 1257.58 s，Serial 1244.37 s）。

随后模式 CLI 又补齐 GATK optional-boolean 表面（bare、`=true|false`、
next-token），并在 mitochondrial mode 的未显式覆盖路径采用 GATK
`log-SNV=-2.5*ln(10)`、`log-indel=-3.75*ln(10)` 默认；verifier 断言有效先验、
microbial `MMQ=20` 及 false/true 解析。最新双 backend 完整 CTest 为
**79/79**（OpenMP 1260.76 s，Serial 1202.73 s）。这些增量仍不改变完整
Mutect2/FilterMutectCalls joint posterior 为 prototype/fallback 的判断。

2026-08-28 增量：MarkDuplicates 已把上述 checkpoint/restart 路径接入 native
实现。回归通过首遍完成后的注入中断验证 sidecar 与 spill run 可在同一输入/选项
下恢复，恢复输出的 duplicate flag、`DT` tag 和 DuplicationMetrics 与一次性运行
一致；OpenMP/Serial 两套 `fastgatk-mark-duplicates-contract` 均通过。由于完整
Picard release-specific metrics、复杂 read-group/optical 规则和大规模真实 BAM
尚未全部覆盖，registry 仍保持 `prototype`。

2026-08-30 Mutect2 线粒体模式增量：native Mutect2 接受 GATK
`--mitochondria-mode` optional-boolean 形式（bare、`=true|false`、next-token），
并在未显式覆盖时设置 `initial-tumor-lod=0`、`tumor-lod-to-emit=0`、
`population-af=4e-3`、
`graph-pruning-log-odds-threshold=-4*ln(10)` 和
`recover-all-dangling-branches`。stats/manifest 输出有效模式、默认值和 pruning
阈值；OpenMP/Serial Mutect2 verifier 及全量 CTest 均通过（79/79，1258.86 s /
1305.92 s）。`OriginalAlignment` 注释、完整 mitochondrial assembly 及
release-specific somatic posterior 仍需 oracle，故 registry 继续保持
`prototype`。

同日 HC CLI 边界增量：显式 `-ERC/--emit-ref-confidence NONE` 现在覆盖 `.g.vcf`
输出后缀的便利推断，保持普通 VCF；`GVCF`/`BP_RESOLUTION` 仍走各自的
reference-confidence 路径，未知枚举值 fail-closed。真实 chr17 fixture 在
OpenMP/Serial 均验证 `NONE` 不产生 `<NON_REF>`，manifest 的
`primary_output_kind=vcf` 与 `compatibility.gvcf=false`，而未提供 ERC 时的
历史后缀推断保持兼容。

同日 dispatcher 兼容增量：HaplotypeCaller registry 已补齐 native parser 中原先
遗漏的 assembly/streaming 参数（`--stream-by-contig`、`--stream-by-region`、
activity/region bounds、downsampling、adaptive graph pruning、linked-graph
toggles、error-correction、genotype-assignment 等）。advanced dry-run 和
`verify_dispatcher.py` 已验证这些参数可以从原始 Nextflow/SLURM 命令到达 C++ Host，
未知参数仍 fail-closed 并保留显式 Java fallback。

同日 HC reference-confidence/assembly 增量：RCM 对 `nInformativeReads=0` 现在使用
GATK all-zero indel PL cache，并执行 least-confident SNP/indel 选择，修复读尾部
`max-indel-size` 窗口的过高 GQ。默认 `min-pruning=2` 的孤立单读 pileup SNP 保留在
PairHMM 输入但在图构建有效且 site QUAL<30 时折回 RCM；Mutect2、显式
`--min-pruning 1` 与短读长无图 fallback 不受影响。chr17 69000--70000 GATK
4.6.2.0 broad oracle 在 OpenMP/Serial 保持 call-set、VCF 核心字段与 gVCF candidate
rows exact；新增 `candidate_low_support_suppressed` telemetry。完整
read-to-haplotype realignment、soft-clip pileup 和全量 posterior numeric parity
仍未宣称 bit-identical。

同日 realignment 增量：新增显式 `--use-haplotype-realignment-for-rcm`。PairHMM
完成后，Host 从每条 source read 的最佳 normalized haplotype 选择稳定 owner，复用
read×haplotype SW traceback，再用 haplotype-vs-reference SW CIGAR 将 matched read
bases 投影回 reference coordinate；haplotype insertion/无法定位的 read base
保持 -1 并逐读回退到原始 CIGAR 投影。RCM 仍只接收 primitive observation arrays，
不把 CIGAR/HTS 对象带入 Kokkos；`rcm_haplotype_realignment_used`、
`rcm_realigned_observations`、`rcm_realignment_fallback_observations` 写入
OutputManifest/telemetry。该开关在 pinned chr17 fixture 上产生 3,265 个真实重映射
observation，开启/关闭的解压 GVCF 文本 hash 相同；OpenMP/Serial 回归脚本
`verify_hc_rcm_realignment.py` 已接入 CTest/verify_all。完整 GATK assembly
haplotype ownership、soft-clip pileup 规则和大规模 oracle 仍需完成，因此默认保持
关闭以维持当前 strict compatibility。

2026-08-31 gVCF block-boundary 修正：ReferenceConfidenceModel 的 GVCFBlockCombiner
现在在相同 GQ band 内仍区分 covered 与 depth-0 run，避免 read 结束后的 uncovered
tail 与前一段 covered locus 合并后把 median `DP` 错写为 0。纯 insertion CIGAR
回归现保留 `17:11-20` covered block，并输出 `17:21-END=32,DP=0`；`verify_indel.py`
及相关 OpenMP/Serial HC contract 已通过；修复后的全量 CTest 在 OpenMP 与 Serial
两个 backend 均为 81/81 通过。

2026-08-30 BQSR 数值路径增量：BaseRecalibrator 仍由 Host 负责 RG/PU、Context、Cycle
和 known-sites 的字符串/坐标语义，但每个 decode batch 的 `CovariateKey` 现在先压缩为
稠密整数 ID，再调用共享 Kokkos `count_bqsr_covariates_kokkos` 的 `RangePolicy` 整数
归约；Host 只在 collect/sink 阶段按稳定 key 顺序合并计数。新增
`BqsrCovariateObservationBatch`/`BqsrCovariateCountResult` primitive API、OpenMP/Serial
API smoke、manifest 的 covariate kernel policy/prepare/execute/observation telemetry，
并在 `verify_bqsr.py` 中断言这些字段；真实 chr17 五张 BQSR 表、ApplyBQSR SAM 核心字段
仍与 Java GATK 逐行/逐记录一致。该改动扩大统一 Kokkos 数值边界，但不改变当前 BQSR
large-corpus、object-store/cloud 语义仍待完成的状态。

同日外部内存兼容增量：SortSam 与 MarkDuplicates 的 `--add-pg-tag`、
`--create-output-bam-index`/`--create-output-variant-index` 现在统一按 GATK optional
boolean 解析（bare/true/1 与 false/0），非法字面量在打开输入或创建 spill 前
fail-closed；对应 registry、SortSam/MarkDuplicates verifier 已增加拒绝非法值的
断言。OpenMP/Serial SortSam 与 MarkDuplicates verifier 及 pinned Picard 语义 oracle
均通过；这改善 SLURM/Nextflow 参数契约，但不改变 Picard 全量 release-specific
metrics/cloud 语义仍为 prototype 的判断。

同日 CLI 边界兼容增量：LeftAlignAndTrimVariants、CombineGVCFs、VariantFiltration、
SelectVariants、FilterMutectCalls 的 `--create-output-variant-index` 统一迁移到
共享 `optional_boolean.hpp`，接受 bare/`true|1` 与 `false|0`，非法字面量在读取
输入前 fail-closed；SelectVariants 的 `--keep-original-ac/dp` 与 FilterMutectCalls
的 mitochondrial/microbial mode 继续使用同一 optional-boolean registry 语义。
五个 native binary 的 `=maybe` 定向检查均返回退出码 2，registry JSON 通过重复键与
语法校验；该增量只收敛 CLI 契约，五个工具仍按其 oracle 覆盖范围保持 prototype。

2026-08-30 BQSR 并发内存增量：`count_bqsr_covariates_kokkos` 改为受内存上限约束的
Kokkos `TeamPolicy` 局部 dense histogram，并按 league 顺序在 Host 确定性合并；
dense key 空间超过 4M 个计数槽（约 32 MiB）时回退到 `RangePolicy` 全局 atomic。
结果对象、benchmark 与 BaseRecalibrator manifest 记录实际 policy、workspace bytes
和 team-local 标志，API smoke 覆盖正常与宽 key-space 回退。OpenMP/Serial smoke、
benchmark verifier、BQSR GATK oracle 均通过，五张报告表和 ApplyBQSR 记录保持
bit-identical；该优化不改变 large-corpus 与 object-store/cloud 语义尚未完成的范围。

2026-08-30 CLI parser consolidation：HC、Base/ApplyBQSR、Mutect2、GenotypeGVCFs、
ReblockGVCF、VariantRecalibrator、ApplyVQSR 及 VCF 工具的 optional-boolean 分支
统一调用 `optional_boolean.hpp`。bare、`=true|false|1|0` 和 next-token 形式拥有
相同语义，`--flag=` 空值在打开输入前拒绝；registry 补齐 ReblockGVCF，15 个
native binary 的 inline/empty/next-token 非法值 contract 在 OpenMP/Serial 均通过；
SortSam/MarkDuplicates/GatherVcfs/DepthOfCoverage 也迁移到该 helper，contract 当前
覆盖 17 个 binary 名称、31 个 option spelling（registry 为 16 个 tool、47 个声明），
受影响的 12 项 OpenMP/Serial CTest 也均为 12/12。该增量只收敛 CLI/dispatcher
表面，不改变完整算法等价与大规模生产语义仍待完成的判断。

同日 GenotypeGVCFs optional-boolean 补齐：`--use-new-qual-calculator`/`--new-qual`
及 `--use-posteriors-to-calculate-qual`/`--gp-qual` 现统一调用
`optional_boolean.hpp`，空 inline、非法 inline 和非法 next-token 均在输入前
fail-closed；registry 补齐四个别名。GenotypeGVCFs contract 与普通、多样本、
多等位、inbreeding GATK oracle 在 OpenMP/Serial 各 5/5 通过，整体仍保持
prototype，未把 CLI 解析修复误报为 joint-genotyping 算法等价。

同日 DepthOfCoverage optional-boolean 补齐：六个 output-omission 开关（含
`--omit-intervals`/`--omit-sample-summary` 旧别名）现统一调用
`optional_boolean.hpp`，非法/空值在打开输入前 fail-closed；DepthOfCoverage 与
扩展 parser contract 在 OpenMP/Serial 均通过。`--do-impute-zeros` 保持必填值
boolean，不属于 optional parser。

同日 registry/dispatcher 收尾：含 optional-boolean 的 16 个 registry 工具均满足
`optional_boolean_options ⊆ flag_options`，JSON 语法/重复声明与 dispatcher dry-run
共 138 项检查通过；GenotypeGVCFs、GatherVcfs、SortSam、MarkDuplicates、
DepthOfCoverage 及共享 parser 的新增双 backend CTest 共 6/6。该项修复关闭
dispatcher 预检误拒绝，整体算法加权完成度仍约 53%。

同日 CalculateContamination 数值收敛：native `ContaminationModel` 现在采用与 GATK
`OptimizationUtils.max` 相同参数的 Apache Commons Math Brent 最大化，并逐字复刻
`MathUtils.binarySearchFindZero` 的 1e-6 误差搜索；四个真实 GATK 4.6.2.0 fixture
在 OpenMP/Serial 的 contamination/error/MAF 严格 oracle 中分别达到
`1.4e-17`、`1.2e-17`、`7.1e-14` 以内。表格输出使用 binary64
shortest-round-trip；segmentation sidecar 复现 Java `HashMap` contig bucket/head-insert
文件顺序，manifest 记录 optimizer 与 error-search；四个 pinned fixture 的
contamination/error 已 raw-bit exact、MAF 差小于 `5.1e-14`，但完整 segmentation
likelihood/report 仍保留 `bit_identical_to_gatk=false` 的显式边界。

2026-08-31 MarkDuplicates 增量：原有输入状态不应影响新一轮 Picard 判定。native 现在在每条输出 record 上先清除旧 `BAM_FDUP`/`DT`，再写入新 duplicate 标记，并为每条输出 record 写入 `PG` 标签；这对齐 GATK/Picard 4.6.2.0 的 `CLEAR_DT=true` 与 `ADD_PG_TAG_TO_READS=true`。独立 pinned oracle `fastgatk-native/scripts/verify_mark_duplicates_gatk_oracle.py` 与现有 spill/restart contract 分离，对比 SAM core fields、duplicate flag、DT/PG、DuplicationMetrics 行和 duplicate-set histogram，并测试 `REMOVE_DUPLICATES` 及仅移除 optical duplicate 的 `REMOVE_SEQUENCING_DUPLICATES`。该路径已在 OpenMP/Serial 构建上通过，但不升级 registry status；UMI/flow/cloud 和 Picard 全量 release-specific 语义仍必须 fallback。


### Scheduler failure/retry/resume contract（2026-08-31）

scatter_gather_smoke.sh 现在对每个 shard 和 gather 使用独立 attempt 目录：
失败 attempt 的 partial output 不会进入最终 outdir；完整 primary/index 与
manifest 按 manifest-last 的 rename 顺序提交。--retries N 提供有界本地/fake-SLURM
模拟重试，重试耗尽时 fail-closed 且不执行 gather。成功 bundle 保存
input/reference/binary/gather binary/threads 及 shard manifest 的 SHA-256 workflow
signature；下一次运行只有在 signature 和 complete manifest 同时匹配时才 resume，
不重算 shard 或 gather。fastgatk-scheduler-retry-contract 用 fake executables
覆盖 shard/gather partial failure、retry、无执行 resume 和 retry exhaustion；它明确
报告 real_slurm_multi_node=false，不能替代真实 SLURM allocation/preemption 测试。

2026-09-01 并行收口：GenomicsDBImport native sparse workspace 的 materialized
input index 现在保留显式 `-V`/sample-map 顺序并稳定去重；reopen/query contract
覆盖 header 顺序和缺失 genotype，但该 workspace 仍不是 TileDB/GenomicsDB 格式。
HC 的 `--dont-use-soft-clipped-bases` 过滤边界与 pinned GATK 四位点 oracle 对齐，
只抑制 soft-clip 派生 assembly evidence。FilterMutectCalls 新增 Kokkos
GermlineFilter posterior 与 Number=A `GERMQ`（NLOD/POPAF/tumor AD，[1,93]，缺
POPAF fail-closed）。CombineGVCFs 新增 `--sites-only-vcf-output`：在完整
allele/FORMAT merge 后用独立 writer header 和 `bcf_subset_format` 发布 8 列 site-only
VCF，reference block 不合成 DP；GATK 4.6.2.0 shape/payload oracle 双 backend exact。

同轮最后三项 bounded slice 已收口：VariantsToTable 修正 `--moltenize` 的
`-ASF/-ASGF` 字段消费顺序、biallelic `Number=R` 投影及 `-EMD` presence 语义；
CalculateContamination 对全部位点被 `MIN_COVERAGE` 过滤的合法输入输出
`0.0/1.0` 与 header-only segmentation，并支持 `-matched/-segments`；GatherTranches
强制 `--mode`、支持 `-tranche`，version-6 header-only scatter shard 作为零行
contributor。三者均有 pinned GATK 4.6.2.0 oracle；VariantsToTable 专用 oracle
使当前 CTest 清单达到 129 项，未改变 registry 的 prototype/adapter 状态。
2026-09-01 CollectAllelicCounts read-filter slice：native 现在按 GATK 4.6.2.0
按 GATK 规则解析过滤器集合：LocusWalker 的 Wellformed+Mapped，加工具附加的
Mapped/NonZeroReferenceLength/NotDuplicate/MappingQuality>=30；secondary、
supplementary 与 QC-fail flag 默认保留。`--disable-read-filter/-DF`、
`--read-filter/-RF` 和 `--disable-tool-default-read-filters` 支持重复/显式
Boolean，并按 GATK 的“先移除默认、再加入显式 filter”顺序解析；Wellformed
边界拒绝 CIGAR `N`。新增 pinned GATK 4.6.2.0 oracle，覆盖 default、全部默认
过滤器禁用、显式重启 NotDuplicate/MAPQ 三种结果；OpenMP/Serial contract 与
oracle 均通过。该项只收敛 pileup/filter boundary，multi-sample/cloud、完整
malformed-input 与 CNV 全链路仍显式 fallback；共享树当前 CTest 为 129 项。

2026-09-01 BQSR/SelectVariants/CNV bounded slices：BQSR 新增 BI/BD per-base
indel-quality 读取、quality/event 分组的 RecalTable1、GATK Bayesian report 对齐与
截短 tag fail-closed；SelectVariants 实现 `--set-filtered-gt-to-nocall`，在
allele compaction 前按 FT 过滤并刷新 AC/AN/AF；CollectAllelicCounts 对齐默认
read-filter 集合及 `-DF/-RF/--disable-tool-default-read-filters` 解析顺序。三个
pinned GATK 4.6.2.0 oracle 均在 OpenMP/Serial 通过，完整 CTest 125/125 通过
（OpenMP 766.81 s、Serial 1024.16 s）。对应 registry 仍保持 prototype；完整
模型、真实 GPU/多节点和云端语义继续显式 fallback。

2026-09-01 CallCopyRatioSegments degenerate-statistics slice：native 现在保留 GATK
4.6.2.0 SimpleCopyRatioCaller 对空/单例 copy-neutral 集合的 IEEE-754 语义：统计
结果为 NaN，非 neutral segments 的比较全部落入 neutral；不再将单例标准差误设
为 0，也不把 outlier 过滤为空回退到原集合。`-O` 默认自动生成 `<output stem>.igv.seg`，
`--legacy-output` 仍可覆盖路径；manifest 对非有限统计写 JSON null。pinned GATK
called/legacy bytes oracle 在 OpenMP/Serial 均通过；full ModelSegments/Cloud semantics
remain explicit fallback。

2026-09-01 SelectVariants FILTER boundary：native 增加
`--exclude-filtered[=BOOL]` 及 `--exclude-filtered-variants[=BOOL]`，按 GATK 语义
处理 bare/true/1、false/0、`PASS`/`.` 与多标签 FILTER 排序；manifest 同步记录
选择值及 `filtered_records`。新的 pinned GATK 4.6.2.0 oracle 与原有 SelectVariants
contract 在 OpenMP/Serial 均通过；ID/pedigree/random fraction/cloud/完整 JEXL 仍为
prototype/fallback。

2026-09-01 GenotypeGVCFs dense reference writer boundary：在
`--gatk-compatible-annotations` 下，`--include-non-variant-sites` 的单态参考位点
遵循 GATK 4.6.2.0 cleanup 语义：MIN_DP 恢复为 DP，GQ 移为 RGQ，并清除 PL/MIN_DP。
新增 pinned Java oracle 对 aggregate 和 `--stream-by-locus` 两条路径的 101 条数据行
逐行 exact；OpenMP/Serial 均通过。默认 native diagnostic profile 继续保留 RCQ/RCP、
GQ/PL，不把兼容 profile 的 writer 变化误报成完整 GenotypeGVCFs 算法等价。
### GetPileupSummaries interval-exclusion boundary (2026-09-01)

Native `GetPileupSummaries` now applies `-XL/--exclude-intervals` after `-L`, supports `-ip`/`-ixp` padding and separate-token Boolean options for default-filter/dictionary switches. The pinned GATK 4.6.2.0 oracle validates normal and padded exclusions plus disabled-default-filter byte identity on both backends.

### ModelSegments supplied partition boundary (2026-09-01)

`ModelSegments --segments` is now a real native input boundary. A validated
Picard interval list (dictionary order, non-overlap) becomes the fixed
modeling partition and bypasses kernel segmentation; copy-ratio observations
use GATK midpoint ownership. The pinned Java 4.6.2.0 oracle checks exact
boundaries, 3+3 midpoint counts, and `.cr.seg` data rows on OpenMP/Serial.
Native bounded posterior/model reports remain explicit fallback scope.

### DenoiseReadCounts HDF5 PoN interval-identity boundary (2026-09-01)

`DenoiseReadCounts --count-panel-of-normals` now performs the GATK 4.6.2.0
pre-denoising SimpleInterval equality check: case and PoN original interval
lists must have identical contig/start/end values in identical order, before
the native path subsets filtered panel intervals. A pinned Java oracle creates
a native HDF5-SVD PoN v7 and verifies valid acceptance plus subset, superset,
and coordinate-mismatch rejection on OpenMP and Serial. This prevents a
successful-looking case from using a different fractional-coverage
denominator; full CNV model/Spark/cloud behavior remains explicit fallback.

### 2026-09-01 parallel bounded slices

BaseRecalibrator now resolves its default read-filter plugin set with GATK's
default-removal/explicit-addition ordering (`--disable-tool-default-read-filters`, `-DF`,
`-RF`); HaplotypeCaller accepts repeated `--kmer-size` and attempts the normalized
`[10,25]` list through the shared Kokkos graph path; FilterMutectCalls uses GATK-bounded
`ROQ` phred serialization `[1,93]`. The independent pinned GATK 4.6.2.0 oracles pass on
OpenMP and Serial. CTest is now 134 entries; these are bounded compatibility gates, not
claims of complete BQSR/HC/Mutect2 algorithm or bit identity.

VariantFiltration additionally implements GATK 4.6.2.0's optional-Boolean
`--invalidate-previous-filters`: only site FILTER labels are reset before new
predicates, while existing genotype FT labels remain. The pinned Java oracle
passes on OpenMP and Serial and records the behavior in compatibility and
telemetry metadata; unsupported full JEXL/annotation semantics stay
fallback-first.

GenotypeGVCFs assignment-method oracle note: pinned GATK 4.6.2.0
`GenotypeGVCFsEngine.createMinimalArgs` forces effective `PREFER_PLS` for
the inherited `--genotype-assignment-method`/`--gam` option. The native
GATK-compatible writer now mirrors this (no GP/PG or posterior-GQ side
effects), while the diagnostic profile continues to exercise richer Kokkos
assignment modes. Java/native rows and aliases pass on both host backends;
this remains an output boundary rather than full cohort parity.

GenotypeGVCFs output-allele subset slice (2026-09-01): after the shared
Kokkos cohort AF posterior, native now applies GATK GenotypingEngine's
per-ALT standard-confidence filter (default 30; `--stand-call-conf` alias)
and remaps Number=A/Number=G PL, AD, and GT fields before publication. The
pinned GATK 4.6.2.0 two-sample unsupported-sibling oracle is registered as
`fastgatk-genotype-gvcf-assignment-gatk-oracle` and passes on OpenMP and
Serial, including `PREFER_PLS` and `USE_PLS_TO_ASSIGN` aliases. Full cohort
calibration and raw bit identity remain explicit follow-up scope.

2026-09-01 full regression closeout: bundled CTest completed the active 134-test
matrix on both CPU execution spaces with OpenMP 134/134 (2184.97 s) and Serial
134/134 (2146.46 s), zero failures. This closes the current host-backend regression
gate only; it does not promote prototype tools to complete GATK algorithm/bit identity,
nor does it substitute for CUDA/HIP/SYCL and multi-node SLURM/Nextflow validation.

Post-closeout bounded additions: the active CTest inventory is now 140. The
ModelSegments supplied-partition, CopyRatio conditional-model, and allele-fraction
conditional-model oracles plus the FilterMutectCalls `-OVI` variant-index alias oracle
each pass on OpenMP and Serial; the prior 134-test full-matrix timings remain the
latest complete-suite baseline. ModelSegments also emits the four
GATK-shaped `.modelBegin/.modelFinal.cr.param/.af.param` sidecars while reporting
`smoothing_applied=false` for the not-yet-implemented Java multi-round smoothing.

GenotypeGVCFs multi-sample reference-confidence slice (2026-09-01): for
`--include-non-variant-sites`, finalized INFO/DP now sums the merged
sample-major FORMAT/DP vector, including disjoint source blocks, instead of
copying the first shard's depth. The pinned GATK 4.6.2.0 oracle compares 101
rows for a two-sample cohort in aggregate and `--stream-by-locus` modes, with
OpenMP/Serial both exact. This closes a production multi-sample annotation
boundary; complete joint-model calibration remains follow-up scope.

GenotypeGVCFs GP-input cleanup slice (2026-09-01): an upstream gVCF may carry
FORMAT/GP and GATK's phred-scaled FORMAT/PG vectors, but the Java
`GenotypeGVCFsEngine` forces `PREFER_PLS`; its `AlleleSubsettingUtils` therefore
removes those pre-`<NON_REF>` Number=G fields before final output. Native now
clears GP/PG at the final GATK-compatible writer boundary, after all cohort
consumers have run, so `--gp-qual` cannot leave stale posterior vectors or
change the PL-based result. The pinned `fastgatk-genotype-gvcf-gp-input-gatk-oracle`
uses a two-sample GP/PG-populated fixture and compares rows exactly in aggregate
and `--stream-by-locus` modes; OpenMP/Serial CTest #71 both pass. The diagnostic
profile intentionally retains posterior fields for explicit posterior-assignment
experiments; full release-specific joint posterior calibration remains fallback.

### CreateReadCountPanelOfNormals degenerate-SVD boundary (2026-09-01)

The native PoN creator now rejects a multi-sample standardized matrix with no
singular value above the GATK `EPSILON` threshold whenever eigensamples were
requested. This matches the pinned GATK 4.6.2.0
`HDF5SVDReadCountPanelOfNormals.create` failure path and prevents an identical
normal panel from being persisted as a silently non-denoising HDF5 container.
`fastgatk-create-read-count-panel-of-normals-degenerate-svd-gatk-oracle` covers
native/Java rejection and partial-output cleanup on both Kokkos CPU backends.

### CombineGVCFs PL-less gVCF payload (2026-09-01)

`CombineGVCFs` now detects whether any merged sample actually carries
FORMAT/PL. For a PL-less reference-confidence cohort, native leaves PL out
of the output FORMAT payload instead of emitting a fabricated missing PL
vector. Pure reference blocks also keep only `END` in INFO; FORMAT/DP is not
copied to INFO/DP. This is implemented at the post-merge writer boundary,
after ALT/sample projection, so PL-bearing remap and genotype-materialization
paths are unaffected. The pinned GATK 4.6.2.0 two-sample oracle
`fastgatk-combine-gvcfs-plless-gatk-oracle` passes on OpenMP and Serial;
complete joint-genotyping posterior/GenomicsDB behavior remains fallback.

### 2026-09-01 新增边界回归

当前 CTest 配置为 152 项：在既有 140/140 CPU 分段基线之外，新增 PL-less
CombineGVCFs（含 `--stream-merge`）、FilterMutectCalls orientation-joint、
CreateReadCountPanelOfNormals 退化 SVD、GenotypeGVCFs maximum-ALT、ModelSegments
smoothing、Mutect2 normal-evidence、FilterMutectCalls normal-artifact、DenoiseReadCounts integer-input、HC AssemblyRegion boundary、HC RCM PL-range、Mutect2 global-mismapping-rate 与 PoN v7 sample-metadata 共十二项 GATK 4.6.2.0 oracle；OpenMP/Serial
均通过。PL-less stream writer 的 `INFO/DP` 差异已在 CTest 首次失败后修复并复核。
这些是可独立复现的直接替换边界，完整 joint-genotyping、release-specific
somatic calibration 与 CNV SVD 生产规模仍不宣称完成。

### ModelSegments credible-interval smoothing slice (2026-09-01)

With smoothing controls explicitly supplied, native applies the GATK
`SimilarSegmentUtils` adjacent-segment rule in a Kokkos `RangePolicy`: posterior
median differences are compared with both 10--90% credible widths, then merged
posteriors use the Java inverse-variance normal approximation. A deterministic
final refit conditions the bounded chain on the merged partition. The pinned
`fastgatk-model-segments-smoothing-gatk-oracle` checks merge and strict
non-merge boundaries against GATK 4.6.2.0 on OpenMP and Serial. Intermediate
refits between multiple smoothing rounds remain a bounded follow-up; the
default native path is unchanged unless smoothing is explicitly requested.

### GenotypeGVCFs maximum-ALT likelihood subset (2026-09-01)

The GenotypeGVCFs joint-locus path now applies GATK's
`AlleleSubsettingUtils.calculateMostLikelyAlleles` semantics before the AF
calculator. `calculate_allele_likelihood_scores_kokkos` scans each sample's
Number=G PL row, chooses the first lowest-PL genotype, and accumulates its
hom-ref PL distance for the contained proper ALTs; Host performs a stable,
deterministic score ordering and the existing Kokkos remappers project all
genotype/allele fields. A dropped-ALT original genotype is retained as a
no-call with missing GQ and a numerically normalized PL row. The pinned
`verify_genotype_gvcf_max_alternate_alleles_gatk_oracle.py` is exact for seven
concrete ALTs reduced to six and covers aggregate plus `--stream-by-locus` on
both CPU backends.

GenotypeGVCFs discovered-ALT annotation slice (2026-09-01): native now
implements the optional GATK `--annotate-with-num-discovered-alleles` Boolean
at the joint-locus writer boundary. INFO/NDA is captured before max-ALT
reduction and counts concrete discovered ALTs (excluding `<NON_REF>`); pure
reference-confidence blocks remain unannotated. The pinned GATK 4.6.2.0
seven-ALT oracle asserts `NDA=7`, exact row/INFO ordering and manifest
telemetry for aggregate and `--stream-by-locus` on OpenMP and Serial.

### FilterMutectCalls NormalArtifactFilter threshold boundary (2026-09-01)

The native NormalArtifactFilter subset now compares the Kokkos-computed NALOD
posterior with the effective FilterMutectCalls error threshold. A posterior
that is positive but below a CONSTANT/learned threshold must remain PASS; the
old `> 0` adapter incorrectly emitted `normal_artifact`. The pinned
`fastgatk-filter-mutect-calls-normal-artifact-oracle` uses GATK 4.6.2.0 with
NALOD=-3, normal AF 0.02, tumour AF 0.10 and CONSTANT=0.1 to fix this boundary.
OpenMP and Serial both pass; full release-specific empirical clustering and
cross-filter posterior calibration remain explicit fallback.

### DenoiseReadCounts integer COUNT input boundary (2026-09-01)

`DenoiseReadCounts` now parses TSV `COUNT` with the same lexical
non-negative signed-32-bit rule as GATK 4.6.2.0 `DataLine.getInt()`. The
native reader no longer accepts `10.0`, `1e1`, negative counts, or values above
`INT_MAX`, which previously could alter the fractional-coverage denominator
or be silently interpreted as a double. Valid integer fixtures retain the
same native/Java standardized values. The pinned
`verify_denoise_read_counts_integer_input_gatk_oracle.py` is registered as
`fastgatk-denoise-read-counts-integer-input-gatk-oracle` and passes on OpenMP
and Serial. The unrelated `--normalization-target-coverage` compatibility
switch remains explicit fallback because it is absent from the pinned
DenoiseReadCounts interface.

### Per-tool progress-score audit (2026-09-01)

The existing progress-score contract now produces a machine-readable
`tool_audits` array for all 48 dispatcher registry entries. Each audit row
contains registry `status`, five gate scores (`api_cli`, `oracle`,
`format_sidecar`, `resource_io`, `e2e_perf`), explicit verifier evidence,
tool-specific SHA-256 fixture/evidence provenance, and a non-empty fallback
reason. The mapping is explicit because registry names and native aliases
differ; missing/stale mappings, absent evidence files, malformed registry
metadata, or a missing pinned digest manifest fail closed. The pre-existing
workflow calculation remains unchanged (`global_score=0.5385`, rounded 54%).

### CreateReadCountPanelOfNormals sample metadata boundary (2026-09-01)

`CreateReadCountPanelOfNormals` now passes the absolute CLI input paths into
the native preprocess stage. This fixes both HDF5 v7 sample arrays:
`original_data/sample_filenames` and `panel/sample_filenames` remain aligned
with their source rows after sample filtering, matching GATK 4.6.2.0 rather
than using `@RG SM` names as filesystem paths. The independent sample-metadata
oracle runs the native and pinned Java writers and checks both path values on
OpenMP/Serial. The mathematical SVD approximation and release-specific
Spark/cloud HDF5 metadata are still explicit fallback boundaries.

### ModelSegments allele-fraction initializer slice (2026-09-01)

For count-backed allelic observations, native now initializes each segment's
minor fraction with the GATK beta-integrated alt-minor/ref-minor responsibility
at `f=0.5`, followed by the Java flat-prior pseudocount. The calculation is a
Kokkos `RangePolicy` over segment offsets and is shared by Serial/OpenMP and
other enabled execution spaces; the result is exposed as
`allele_fraction_initial_segment_means` in the manifest. The pinned
4.6.2.0 oracle validates `(ALT,REF)=(2,1),(1,2)` → `29/64`, plus Java/native
partition and output schema. This materially improves the initial state while
the release-specific slice sampler and latent Gamma-bias MCMC remain fallback.

### ApplyVQSR default cutoff boundary (2026-09-01)

The pinned GATK 4.6.2.0 threshold oracle fixes a direct-replacement edge:
`--tranches-file` alone does not activate tranche walking. Without an explicit
truth-sensitivity target, Java uses effective cutoff `0.0` and `LOW_VQSLOD`.
Native now follows this path; `fastgatk-apply-vqsr-default-cutoff-gatk-oracle`
passes on both CPU backends. Explicit target values retain tranche semantics.

### CallCopyRatioSegments interval validation boundary (2026-09-01)

Native now rejects overlapping or unsorted segment partitions during parsing,
matching GATK `CopyRatioSegmentCollection`. The pinned
`fastgatk-call-copy-ratio-segments-interval-validation-gatk-oracle` passes on
both CPU backends. Full CNV/ModelSegments parity remains explicit fallback.

### VariantRecalibrator Java VBEM update boundary (2026-09-01)

Native `VariantRecalibrator` now follows the pinned GATK 4.6.2.0
`GaussianMixtureModel` control flow for the fitted model: Java-Random seeded
mu/Sigma initialization, K-means empty-cluster random restart, empirical
`I/200` covariance prior, zero empirical mean after normalization, Wishart
shrinkage with one shrinkage factor on the mean outer product, raw pMix
normalization, and post-E-step pMix convergence. Normalized rows beyond
`standard-deviation-threshold` are retained for output but excluded from
training and worst-variant selection, matching `VariantDataManager`. The new
pinned `fastgatk-variant-recalibrator-vbem-gatk-oracle` compares a two-Gaussian
positive/negative model to Java 4.6.2.0: PMix/mean/covariance maximum delta is
`5.23e-9` and six VQSLOD values differ by at most `5e-5`; the old squared
shrinkage-factor path differed by `1.22744`. It passes on OpenMP and Serial.
Java raw-bit VBEM and full resource/recalibration-table provenance
remain separate fallback boundaries.

### VariantRecalibrator zero-variance normalization boundary (2026-09-01)

GATK 4.6.2.0 rejects a training annotation when its standard deviation is
below `1e-5` during `VariantDataManager.normalizeData`. Native now enforces
the same fail-closed boundary in `compute_annotation_normalization` instead of
substituting unit scale and publishing a model with different likelihood
geometry. The pinned constant QD/MQ Java/native oracle and CTest #106 pass on
OpenMP and Serial; this does not claim full VQSR raw-bit VBEM parity.

### ApplyVQSR sites-only writer boundary (2026-09-01)

GATK 4.6.2.0's `--sites-only-vcf-output` is a writer-only optional Boolean.
Native ApplyVQSR now keeps FORMAT/sample state through score/filter evaluation,
then strips it only at publication using a duplicated writer header. The
`fastgatk-apply-vqsr-sites-only-gatk-oracle` compares eight-column output,
VQSLOD/FILTER rows, and manifest metadata on both host backends.

### VariantFiltration mask reverse boundary (2026-09-01)

Native `VariantFiltration` now covers GATK's indexed VCF/BCF mask direction
contract. `--filter-not-in-mask` reverses the normal overlap test, so records
outside the mask receive the mask FILTER label while overlapping records pass;
existing site FILTER labels remain intact. `--mask-description` emits the
custom FILTER header description, with GATK's normal/reverse default text when
omitted. The pinned Java 4.6.2.0 oracle indexes the mask via
`IndexFeatureFile` and compares both directions, chained FILTER history,
custom header description, and manifest metadata on OpenMP/Serial. Complete
JEXL/annotation-engine, cloud mask, and complex allele-specific semantics stay
explicit fallback.

### CombineGVCFs interval/reference-block boundary (2026-09-01)

Pinned GATK 4.6.2.0 confirms that a `REF/<NON_REF>` block whose `END` lies
after an interval end is dropped (with Java's "cuts in the middle" warning),
while a block whose `END` is contained is emitted in full. Native now applies
this rule before both aggregate and `--stream-merge` traversal. The independent
`fastgatk-combine-gvcfs-interval-refblock-gatk-oracle` covers the drop and
contained-END cases on both backends.
### Mutect2 `--normal-lod` (2026-09-01)

Implemented the matched-normal normal-evidence emission gate with GATK's
default 2.2. The native log10 normal likelihood boundary is applied only to
ALT calls with a matching normal candidate and is recorded in stats/manifest.
`fastgatk-mutect2-normal-lod-gatk-oracle` uses a deterministic synthetic BAM
fixture and pinned Java 4.6.2.0 to check the low-threshold emit and zero
threshold suppress cases on both Kokkos OpenMP and Serial.

### ModelSegments default kernel dispatch and linear feature (2026-09-01)

GATK selects `KernelSegmenter` by default when no advanced `--segments`
partition is supplied. Native now follows that dispatch, while retaining the
native `--change-point-threshold` extension as an explicit legacy threshold
mode. The copy-ratio-only linear kernel uses its exact rank-one feature
(`x -> x`), matching the SVD of GATK's linear outer-product kernel and making
the changepoint cost independent of random anchors. The pinned
`fastgatk-model-segments-default-kernel-gatk-oracle` covers multi-segment step,
gradient, and blip fixtures on OpenMP and Serial; Gaussian/multidimensional
projection and full Java SVD/MCMC bit identity remain explicit fallback.

### GenotypeGVCFs spanning-deletion ownership boundary (2026-09-01)

Pinned GATK 4.6.2.0 confirms that an orphan `*` ALT is removed after the
joint AF/QUAL calculation, while a `*` covered by a concrete deletion span is
preserved. Native now tracks deletion spans across shards and performs the
bounded final union/remap in both aggregate and `--stream-by-locus` paths;
affected genotypes become no-call with normalized PL/GQ. The dedicated
`fastgatk-genotype-gvcf-spanning-deletion-gatk-oracle` is exact on OpenMP and
Serial. Full release-specific joint posterior calibration remains fallback.

### DepthOfCoverage read-filter boundary (2026-09-01)

Implement the GATK 4.6.2.0 default Wellformed + Mapped + NotDuplicate +
NotSecondary set in the Host projection boundary. Preserve QC-fail and
supplementary records, and reject whole records containing CIGAR `N`, as
required by CigarContainsNoNOperator. The pinned Java/native oracle covers
both cases on OpenMP and Serial; broader gene/cloud and release-specific
coverage reports remain fallback.

### DepthOfCoverage deletion-site boundary (2026-09-01)

GATK 4.6.2.0 treats `--ignore-deletion-sites` as an effective override of
`--include-deletions`: deletion-only loci are counted when include-deletions
is active, but contribute zero depth when the ignore switch is true. Native
now applies the same rule in its Host-to-Kokkos projection and records both
requested and effective settings. The pinned
`fastgatk-depth-of-coverage-ignore-deletion-sites-gatk-oracle` is exact on
OpenMP and Serial.

### Nextflow output/resource boundary correction (2026-09-01)

The generic DSL2 compatibility process now declares output files from the
caller-selected `params.output` stem.  The former fixed `result*` declaration
could report a successful native command as a missing Nextflow output when an
existing process used another filename.  The local verifier now exercises
`custom-counts.txt`.  The scheduler config also exposes optional
`container`/`cpus`/`memory` settings and a `slurm_gpu` profile; GPU requests are
translated to SLURM `--gpus` but do not claim a CUDA/Kokkos device runtime.
`slurm_smoke.sh` accepts `FASTGATK_GPU_COUNT/TYPE` and an explicit
`FASTGATK_CUDA_VISIBLE_DEVICES` override, with fake-sbatch coverage for mapping
and invalid values.  Real container, GPU runtime, and multi-node/preemption
validation remain deployment gates.

### ModelSegments first-site oriented AAF (2026-09-01)

组合 CR+AF 分段现在为每个 copy-ratio interval 保留坐标序第一个 het site 的
oriented ALT fraction，并与 posterior/modeling 使用的 folded MAF 分离；single/
multisample Kokkos kernel 均消费该独立信号。固定 GATK 4.6.2.0 fixture 在区间
平均 AAF 恒 0.5、MAF 恒 0.1 时仍精确复现 `1-1600/1601-3200` 边界。完整
per-locus AF posterior 与 Java MCMC/report raw-bit parity 仍为后续项。

### ApplyVQSR filter-control optional Boolean（2026-09-01）

`--ignore-all-filters` 与 `--exclude-filtered` 已按 pinned GATK 4.6.2.0
改为 optional Boolean：bare/`true`/`false` 都有独立 Java/native oracle。
前者决定已有 FILTER 的输入是否重新评分，后者决定新产生的 VQSR FILTER 是否
写出；`false` 不再被误当作未知参数或 bare `true`。完整 VQSR model/provenance
仍保持 fallback。
