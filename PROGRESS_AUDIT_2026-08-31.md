# Fast-GATK 当前进度审计（2026-08-31）

> **A-line 已知限制（2026-09-03，Round 18-36 专项诊断结论）**：HC 在**连续
> 多 AssemblyRegion 单 -L 窗**下的 call 完整性与等位正确性存在架构级缺陷：
> 首 region 噪声吃满全局候选预算（64）、长图窗单种子 DFS 深度不够、以及
> 连续链式多窗的等位合并错误（实测 20228 真实 ALT=C 被误判为 G>A 等）。
> 分离 interval / 隔离单窗行为与 read 证据及 GATK 4.6.2.0 一致；full 窗
> 修复曾出现错误 call 回归，故所有实验性修补已回退，保持 canonical
> 210/210 × 双后端基线。缺陷已完整归档（NEXT_PHASE_TASKS.md R18-36 证据链），
> 修复需架构级重构（逐 AssemblyRegion 独立图/装配/基因型 + 窗间等位合并），
> 属 A-line 专项；在修复前，多 region 连续窗输出以 GATK/read 证据复核为准。

> **Round 5–11 收口（2026-09-03）**：registry 全部 48 工具定稿——**47
> contract-compatible / 1 adapter / 0 prototype**（R5 批量晋级 BQSR 四件套、
> Mutect2 链路五件、CNV 五件、VQSR 两件等最后 16 个工具，每工具补
> `compatibility_scope` + `fallback_boundaries` + 事实 note；`verify_dispatcher.py`
> 142 项检查通过；双后端 progress-score/fixture-digest CTest 通过）。同期定位并
> 解决 Mutect2 SIGSEGV 事故：主 `build/` 为跨 08-21~09-02 增量构建的陈旧对象
> 混链（Makefile 不追踪 flag/头文件变更），全量重建两后端后 OpenMP 210/210、
> Serial 210/210（修复 api_smoke 用例漏填 indel_counts 字段与 jshell prefs 只读
> HOME 两处测试代码问题）。**B1（VariantRecalibrator VBEM raw-bit）完成**：
> commons-math3 3.5 FastMath.log/digamma 精确移植（`vqsr_fastmath.hpp`，探针
> 17/17 位级一致）+ GATK 捆绑 JAMA 1.0.3 LU decompose/det/solve 移植
> （`vqsr_jama_lu.hpp`，随机 SPD 探针 300/300 位级一致）+ M-step 结合序/
> E-step 与 score kernel log10 域组合 + **VQSLOD 按 GATK `%.4f` 文本写出**
> （GATK 以字符串而非 float 写 VQSLOD）→ pinned fixture **VQSLOD delta 0.0**，
> oracle 收紧为 token 精确相等，双后端 7/7 VR oracle 通过；模型 report 参数残余
> ~5e-9（≤1e-6 gate，正模型迭代路径边际 raw-bit 化收益低，记录为已知边界）。
> 剩余主线为 Track A（A1 装配图 raw-bit、A2 PairHMM flow、A3/A4/A5 后验与
> joint、A6 WGS benchmark）、B2–B5、Track D（GPU/SLURM 需相应环境）。

> **第三波收口（2026-09-01）**：PairHMM fragment-first 聚合新增 Kokkos 有序
> fragment×haplotype 求和，保留 Java `-Infinity`（零概率 mate）并在后续
> max-marginalization 中继续消费；六 cell pinned GATK 4.6.2.0 oracle 在 OpenMP/
> Serial raw-bit exact。ModelSegments 仅修正 manifest 的 total/burn-in/retained
> provenance，并明确 conditional-point estimate 与完整 Gibbs/MCMC/RNG parity 的
> 边界，未改变数值 baseline。GenomicsDBImport native workspace 明确声明
> `fastgatk-portable-sparse-v1`、非 TileDB、非 GenomicsDB query-compatible；真实
> TileDB 仍由 external adapter/legacy bridge 提供，native-storage boundary oracle
> 双后端通过。现场 CTest 清单已更新为每个 backend 209 项（新增 somatic likelihood
> oracle 条目后复核）。

> 2026-09-01 ApplyBQSR ContextCovariate 报告维度：BaseRecalibrator 虽已把
> `mismatches_context_size` 写入 GATKReport，consumer 过去仍硬编码 2-mer，读取
> 非默认报告时会错失所有 Context rows。现从 `Arguments` 严格解析 1..13 并用于
> ApplyBQSR Host key construction，随后复用 Kokkos RangePolicy quality transform。
> 新增 `fastgatk-bqsr-context-size-gatk-oracle`：size-3 synthetic 先比较五张报告表，
> 再比较 200/200 Java-decoded records；fixture 的 context posterior 确实改变 QUAL，
> 因而不能靠空表或取整偶合通过。完整自定义 covariate/plugin 仍 explicit fallback。

> 2026-09-01 HaplotypeCaller/Mutect2 CIGAR indel activity ownership：发现 GATK 会把
> 仅由 CIGAR `I/D` 表达的 indel 作为 active locus；native 原先只在 mismatch pileup
> 中播种 activity，候选因此落入 fallback 并可能被 PairHMM 全部 disqualify。现由统一
> Host 投影 indel anchor 到 Kokkos ActivityProfile；交替重复 deletion 左规范化也按等价
> edit 左移。新增 `fastgatk-hc-cigar-indel-activity-gatk-oracle`，OpenMP/Serial 均
> `pass`：8x `60M2D58M` + 8x `120M` 输出与 GATK 4.6.2.0 相同的 `chr1:250 GCA>G`、
> GT/AD/DP/GQ/PL/QUAL，且 assembly region owned、PairHMM unassigned=0。该项只收敛
> activity/ownership，不宣称完整 HC/Mutect graph parity。

> 2026-09-01 FastaAlternateReferenceMaker IUPAC 边界：GATK 对
> `--use-iupac-sample` 的纯合 genotype 选择 genotype 对应 base；native 已修正
> 多 ALT 2/2 从 ALT[0] 错误回退为 ALT[1] 的差异。pinned Java/native oracle
> `fastgatk-fasta-alternate-iupac-hom-gatk-oracle` 在 OpenMP/Serial 均通过，覆盖
> hom-ref、hom-ALT 与多 ALT 选择；复杂 indel/no-call 仍显式不扩大承诺。

> 2026-09-01 AnnotateIntervals feature-query boundary：native 新增并记录
> GATK `--feature-query-lookahead`（含 0/负值），现有 pinned interval Java/native
> oracle 双 backend 复核通过；该参数只控制 FeatureManager 查询缓存，不改变 annotation
> 数值，云端 track/query 仍按 registry 显式 fallback。

> 2026-09-01 VariantRecalibrator annotation-order gap：native 现按 GATK
> `VariantDataManager` 在归一化后的 training/non-training 均值偏移对 annotation 维度做稳定
> 重排，再进入 Kokkos VBEM；避免 CLI 顺序改变多高斯随机协方差初始化。新增 pinned
> `fastgatk-variant-recalibrator-annotation-order-gatk-oracle`，OpenMP/Serial 均 pass，
> VQSLOD 最大绝对差约 4e-5。完整多高斯 VBEM/raw-bit parity 仍显式 fallback。

> 2026-09-01 BaseRecalibrator 长读边界：GATK CycleCovariate 默认上限为 500；native
> 新增 `--maximum-cycle-value/--max-cycle`，超限时 Java-equivalent fail-closed，并将
> 该维度写入 checkpoint/report/manifest。BQSR projection 修正 CIGAR insertion bases：
> 按 GATK `isSNP=0` 纳入 substitution covariates。新增 `fastgatk-bqsr-long-read-gatk-oracle`
> 对 1,500M 与 `50M2I48M` synthetic 做默认失败、显式 1,500 的 Quantized/RecalTable0/1/2
> 全表比较；OpenMP/Serial 均 pass。该项不扩大完整云端分布式 BQSR parity 承诺。

> 2026-09-01 PairHMM indel-quality boundary：GATK 无 BI/BD 标签时原始质量为 Q45；
> Conservative PCR 以原始 read 坐标仅调整非终位，clipping 后末位若原本是内部重复位点
> 可为 Q20/Q40，`NONE` 模式保持 Q45。native HC/Mutect2 已共用该规则，706-row
> pinned 4.6.2.0 输入 oracle 双后端通过，显式 BI/BD 路径未变。

> 2026-09-01 ModelSegments 多样本 kernel：修正先取样本均值导致反相关事件抵消的错误，
> 改为完整 per-sample feature vector（线性拼接、Gaussian anchor 求和）并复用 Kokkos
> RangePolicy。pinned GATK 4.6.2.0 anti-correlated/concordant synthetic 在 OpenMP/
> Serial 分段边界 exact；完整 posterior/SVD raw-bit parity 仍显式 fallback。

> 2026-09-01 VariantRecalibrator retry/iteration boundary：pinned GATK 4.6.2.0
> 将 `--max-attempts`（模型构建重试，默认 1）与 `--max-iterations`（VBEM 迭代
> 上限，默认 150）明确分开；native 已分别传递并在 OutputManifest 记录。新增
> `fastgatk-variant-recalibrator-attempt-iteration-gatk-oracle`，OpenMP/Serial 均通过。
> 该项仅收敛参数/重试契约，不宣称完整 VBEM/raw-bit 模型 parity。

> 2026-09-01 VariantRecalibrator zero-variance guard：训练 annotation 的标准差低于
> GATK 阈值 `1e-5` 时，native 现在在模型拟合前 fail-closed，不再以单位尺度继续训练。
> pinned GATK 4.6.2.0 常量 QD/MQ fixture 在 OpenMP/Serial 的 Java/native 退出边界均
> exact；新增 CTest #106。

> 2026-09-01 LearnReadOrientationModel orientation CLI boundary：native 已支持
> GATK `--QUIET/--quiet`、`--tmp-dir`、`--verbosity/--VERBOSITY` 与 JDK
> deflater/inflater Boolean aliases；两样本 tar writer pinned oracle 在 OpenMP/Serial
> 均通过，manifest 记录 utility controls，未扩大 gzip/member byte parity 承诺。
> 同一 pinned 双样本 fixture 还覆盖 `--max-depth 2`：源 histogram 的高深度 bin
> 保留用于 examples 计数，仅跳过 EM 观测，native 与 Java 最大概率差约 1.9e-12。

> 2026-09-01 CollectReadCounts CNV CLI boundary：native 已支持 GATK utility 的
> `--QUIET`、`--read-validation-stringency/-VS`、`--disable-bam-index-caching/-DBIC`、
> `--tmp-dir` 与 JDK codec Boolean/value 参数；TSV interval/exclusion Java oracle 在
> OpenMP/Serial 的 contract+oracle 均通过，manifest 记录 validation/utility choices。

> 2026-09-01 GatherBQSRReports CLI boundary：native 支持 GATK utility 的
> `--QUIET/--quiet`、`--tmp-dir`、`--verbosity`、JDK deflater/inflater Boolean/value
> 参数；pinned GATK 4.6.2.0 五表 oracle 在追加这些 controls 后仍通过，输出表保持
> bit-identical，空报告与不兼容维度继续 fail-closed。

> 2026-09-01 SortSam CLI boundary：native 已支持 Picard/GATK 4.6.2.0 的大写
> input/output/reference、sort-order、bounded spill、CREATE_INDEX、QUIET、JDK codec、
> compression、validation 与 verbosity 参数。pinned Java/native oracle 比较坐标
> header/records，并验证 compression/validation telemetry 及 `CREATE_INDEX=false`；
> OpenMP/Serial 均通过。完整 Picard sort-order/cloud/MD5/byte-identical 语义仍为
> explicit fallback。

> 2026-09-01 GatherVcfs CLI boundary：native 已支持 GATK/Picard 大写
> `--INPUT/--OUTPUT/--REFERENCE_SEQUENCE/--QUIET/--VERBOSITY`，并保留
> `--COMMENT/--CREATE_INDEX/--COMPRESSION_LEVEL` writer controls。pinned GATK
> 4.6.2.0 oracle 在 OpenMP/Serial 对两 shard 数据行、comment、compression=0 和
> 禁用 index 均通过；完整 cloud/TMP/MD5 语义仍保持 explicit fallback。

> 2026-09-01 LeftAlignAndTrimVariants CLI/interval boundary：native 已支持重复
> `-XL/--exclude-intervals`、`-ip`、`-ixp`，使用记录 END span 做 include/exclude
> padding 与过滤，并补齐 `-isr/-imr/-OVI` aliases 及 trim/split/keep-original 的
> separated Boolean。pinned GATK 4.6.2.0 oracle 在 OpenMP/Serial 均通过，manifest
> 记录 padding、跳过计数和 index policy；复杂 symbolic/cloud/bit-identical
> normalization 仍是 explicit fallback。

> 2026-09-01 AnnotateIntervals validation boundary：native 现在与 GATK
> 4.6.2.0 一致要求 `--interval-merging-rule OVERLAPPING_ONLY`、
> `--interval-set-rule UNION`，并拒绝非零 `--interval-padding` /
> `--interval-exclusion-padding`；include/exclude selector 与 Java 错误诊断已在
> OpenMP/Serial contract/oracle 复核通过。该项只收敛 CNV interval 参数校验，不改变
> cloud/remote-output fallback 或完整 CNV 模型边界。

> 2026-09-01 ValidateVariants/LeftAlignAndTrimVariants bounded CLI 增量：
> ValidateVariants 对 symbolic ALT（如 `<DEL>`/`<CNV>`）不再执行 concrete ALT 的
> called-genotype usage 误报，真实 concrete ALT 仍严格校验；LeftAlignAndTrimVariants
> 新增 interval include/exclude、padding、optional Boolean 与 index 组合边界，均由
> pinned GATK 4.6.2.0 oracle 在 OpenMP/Serial 通过。两项只收敛已证明的 CLI/validation
> 边界，复杂 symbolic normalization 与完整 Tribble 语义仍是 fallback。

> 2026-09-01 ValidateVariants symbolic-ALT boundary：native 已修正
> `VariantContext.validateAlternateAlleles()` 语义，仅对 concrete ALT 要求其出现在
> called genotype 中；`<DEL>`、`<CNV>` 等 symbolic ALT 即使为 0/0 或 no-call 也不再被
> 错误拒绝。pinned GATK 4.6.2.0 symbolic/mixed/concrete-failure oracle 在 OpenMP/
> Serial 均通过，具体 ALT 的严格失败检查保持不变。

> 2026-09-01 VariantFiltration genotype no-call boundary：native 新增 GATK
> `--set-filtered-genotype-to-no-call`，支持 bare/分离 Boolean（embedded `=` 与
> GATK 一致 fail-closed）；带非-PASS FORMAT/FT 的 called genotype 改为 unphased `./.`，
> 并仅在 genotype 实际变化时重算 AC/AN/AF。pinned GATK 4.6.2.0 oracle 覆盖新旧 FT、
> bare/true/false、计数和 writer 语义，OpenMP/Serial 均通过。

> 2026-09-01 SelectVariants sites-only writer boundary：GATK 4.6.2.0 的
> `--sites-only-vcf-output` 在最终 writer 移除 FORMAT/sample；native 新增独立
> writer header 和 `bcf_subset_format`，选择、样本/等位基因重映射仍在完整 record 上
> 完成。pinned Java/native oracle 比较 8 列 site rows 与 false 模式完整样本输出，覆盖
> bare、分离 true/false 及 embedded `=` fail-closed；OpenMP/Serial 各 `1/1` 通过，manifest
> compatibility/telemetry 均记录该选择。

> 2026-09-01 ApplyVQSR exclusion interval：native 新增 GATK `-XL/--exclude-intervals`
> 与 `-ixp/--interval-exclusion-padding`。先应用 `-L` include，再按 exclusion union
> 移除记录 span；仅给 `-XL` 时仍遍历隐式全输入。pinned GATK 4.6.2.0 Java/native
> oracle 覆盖 include+exclude、exclude-only、padding，site/filter/VQSLOD 行与 manifest
> telemetry 一致；OpenMP/Serial CTest #94 均 Passed。该项只收敛 ApplyVQSR interval
> traversal，不改变完整 VQSR 模型的 fallback 边界。

> 2026-09-01 HaplotypeCaller sites-only writer：native 新增 GATK
> `--sites-only-vcf-output` optional Boolean。VCF/GVCF 先完整消费 FORMAT/sample
> 证据完成 genotyping/RCM，再在最终 writer 输出 8 列 site-only 记录；false/bare/
> separated forms、manifest telemetry 与 pinned GATK 4.6.2.0 Java/native oracle
> 已在 OpenMP/Serial 通过。GVCF reference-block 分段差异仍不被该 writer gate
> 掩盖，concrete candidate rows 保持 exact。

> 2026-09-01 HaplotypeCaller GVCF `--floor-blocks`：native 现在对普通 reference
> block 执行 GATK 4.6.2.0 的实际 writer 语义：GQ 向下取到 `-GQB` 分箱下界并输出
> `GT:DP:GQ`，移除 `MIN_DP/PL`；candidate 和 BP_RESOLUTION 行不变。扩展的 pinned
> GQ-band oracle 对 header、行、FORMAT/sample shape 和 manifest 分别在 OpenMP/Serial
> 通过；它不扩大完整 RCM/bit-identical 声明。

> 2026-09-01 Mutect2 `--normal-lod` emission boundary：native 新增 matched-normal
> ALT log10-odds 阈值解析（默认 2.2），在 somatic likelihood 物化后控制候选输出，
> 不改变 PairHMM/后验计算。pinned Java/native fixture 验证低阈值保留位点、零阈值
> 抑制位点及 stats/manifest telemetry；OpenMP/Serial oracle 均通过。完整
> SomaticGenotypingEngine posterior calibration 仍保持 explicit fallback。

> 2026-09-01 GenotypeGVCFs exclusion interval：native 新增 GATK `-XL`/`--exclude-intervals`
> 解析与 Host interval predicate；先求 `-L` include set，再移除与 exclusion union
> span 相交的记录，aggregate/`--stream-by-locus` 一致。pinned GATK 4.6.2.0 fixture 将
> concrete variant 放入排除区间，Java/native 两种 traversal 的数据行和 manifest 均
> exact，OpenMP/Serial 通过。该项仅收敛 CLI/interval traversal 边界，不宣称完整
> GenotypingEngine joint posterior 等价。

> 2026-09-01 CollectF1R2Counts/GatherPileupSummaries CLI alias：核对 GATK
> 4.6.2.0 Barclay help 与实际命令后，确认两工具使用 `--I`/`--O` 大写长别名；
> native parser、dispatcher registry 已补齐，现有 Collect archive contract、Gather
> table contract/GATK oracle 均增加 alias 调用并在 OpenMP/Serial 复核。该项只收敛
> 直接替换的参数名边界，不扩大 F1R2 cloud/read-filter 与 Gather cloud/dictionary
> fallback 承诺。

> 2026-09-01 CollectF1R2Counts exclusion interval：native 现实现 `-XL/--exclude-intervals`
> 的 locus-level subtraction；跨过排除 span 的 read 仍可为相邻 site 提供观察。pinned
> GATK 4.6.2.0 oracle 比较 `.ref_histogram`、`.alt_histogram` 和 `.alt_table` 成员内容，
> OpenMP/Serial 均通过。无 `-L` 的 whole-reference exclusion 仍明确 fail-closed。

> 2026-09-01 LeftAlignAndTrimVariants sites-only writer：native 已补齐
> `--sites-only-vcf-output` optional Boolean。规范化、拆分和 genotype/annotation
> 投影仍消费完整样本字段，最终 writer 才裁剪 FORMAT/sample，Java/native 输出均为
> 8 列 site-only shape；pinned GATK 4.6.2.0 oracle 覆盖 bare、true、false，OpenMP/
> Serial 均通过，registry/verify_all/manifest 已同步。复杂 symbolic/bit-identical
> normalization 仍保持 explicit fallback。

> 2026-09-01 BQSR CLI alias：BaseRecalibrator 新增 GATK 4.6.2.0 `-ics`/`-mcs`
> covariate 短别名；BI/BD indel oracle 使用 `-ics 4` 对 Java/native 四张报告表
> 逐行比较，OpenMP/Serial 均通过。PU→RG 归一与可选 I/D 事件模型已有独立 oracle，
> 本次未重复修改。

> 2026-09-01 CallCopyRatioSegments z-score 参数边界：GATK 4.6.2.0 的完整参数名
> `--outlier-neutral-segment-copy-ratio-z-score-threshold` 和
> `--calling-copy-ratio-z-score-threshold` 已加入 native parser/registry；pinned
> Java/native oracle 与现有 contract、interval-validation CTest 在 OpenMP/Serial
> 均通过。短别名仍兼容，完整 ModelSegments segmentation 未因此宣称完成。

> 2026-09-01 VariantFiltration missing-value Boolean boundary：GATK 4.6.2.0
> 将 `--missing-values-evaluate-as-failing` 声明为 Boolean，支持 bare/分离 true/false，
> 且拒绝 embedded `=`；native 已统一其 separated-Boolean parser。pinned Java/native
> oracle 同时验证缺失 INFO 数值与 `== null` 控制路径，OpenMP/Serial 均通过。

> 2026-09-01 DenoiseReadCounts 双输出 writer：GATK/Picard 4.6.2.0 的
> `--denoised-copy-ratios` 与 `--standardized-copy-ratios` 都是 required。native
> 原先允许只给 `-O`，现已在读取输入前 fail-closed；integer-input oracle 新增
> Java/native missing-output guard，并同步更新现有 contract、dispatcher 和 benchmark。

> 2026-09-01 DenoiseReadCounts HDF5 metadata：SimpleCountCollection 的 HDF5 输入将
> 完整 sequence dictionary 放在 `/locatable_metadata/sequence_dictionary`，GATK 会将
> 其中全部 `@SQ`（含 AS/M5/UR/SP）写入 standardized 与 denoised TSV。native 已在
> HDF5 解码后保留并传播该 header；新增 pinned GATK 4.6.2.0 oracle 比较两份输出的
> header 与数据行，OpenMP/Serial 均通过。其余 HDF5/Spark/cloud metadata 仍是显式
> fallback。

> 2026-09-01 GatherVcfs Picard alias boundary：GATK 4.6.2.0 的 GatherVcfs 暴露
> `--REFERENCE_SEQUENCE`（`-R`）大写长别名；native 已补齐该合法拼写，并以同一
> pinned Java/native shard gather oracle 验证记录一致。

> 2026-09-01 LearnReadOrientationModel 多样本 writer：已修复标准
> CollectF1R2Counts tar 中多个样本 histogram/table 被按 context 错误碰撞的问题，
> native 现在逐样本输出 `<sample>.orientation_priors`。pinned GATK 4.6.2.0 双样本
> oracle 在 OpenMP/Serial 均通过；输出成员内容在 1e-10 EM 容差内一致，tar 前缀/顺序
> 与 gzip bytes 不作 bit-identical 承诺。

> 2026-09-01 MarkDuplicates `CREATE_INDEX` 默认边界：GATK/Picard 4.6.2.0 默认不
> 创建 BAM index；native 原先默认创建 `.bam.bai`，已修正为 false。pinned
> Java/native verifier 新增默认无 index 与显式 true 有 index 的副作用检查，
> OpenMP/Serial `fastgatk-mark-duplicates-contract` 和 `...-gatk-oracle` 均通过。

> 2026-09-01 MarkDuplicates `TAGGING_POLICY=OpticalOnly`：已修复一个实际语义缺口，
> native 原先只接受 `DontTag|All`。现在 optical duplicate 写 `DT:Z:SQ`，普通
> library duplicate 不写 `DT`，且 duplicate flags、metrics、histogram 仍复用同一
> classification 路径。pinned GATK 4.6.2.0 Java/native 14-record oracle 与双
> backend 独立 CTest 均通过；registry/dispatcher/verify_all 已同步。UMI/flow/cloud
> 等完整 Picard 语义仍需 explicit fallback。

> 2026-09-01 MarkDuplicates paired ReadEnds key：已修复无 `MC` 标签及 terminal
> S/H clipping 时 native 仅依赖 mate fields、将一对 mate 拆成两个 duplicate
> fragments 的问题。native 现在首遍按 read name 暂存 mapped mate metadata，成对后
> 使用两条实际 unclipped 5′ 坐标构造统一 key，并将必要的 key metadata 带入 spill
> checkpoint/resume。pinned GATK/Picard 4.6.2.0 4-record oracle（无 MC、2S4M/4M2S）
> 与双 backend CTest 均通过；该项不改变 UMI/flow/cloud 等 explicit fallback。

> 2026-09-01 PreprocessIntervals interval-rule validation：GATK 4.6.2.0 虽在通用
> 参数中将 interval merging rule 默认显示为 `ALL`，但该工具的 validateArguments
> 对省略参数和显式 `ALL` 均失败，并要求 `OVERLAPPING_ONLY`。native 已补齐默认/显式
> `ALL` 的 fail-closed 行为；现有 contract 增加 pinned Java/native 双边界断言，
> OpenMP/Serial 均通过。

> 2026-09-01 GatherBQSRReports empty-report boundary：GATK 4.6.2.0 的
> `RecalibrationReport.gatherReports()` 会过滤无 recalibration row 的报告；当所有
> 输入仅含 Arguments/Quantized、RecalTable0/1/2 均为空时，以
> `there is no usable data in any input file` 失败且不发布输出。native 原先会发布
> 误导性的零行 merged report，现已在读表后按 Java `isEmpty` 语义 fail-closed，且不
> 写 primary/sidecar。pinned Java/native oracle、现有 Gather CTest 在 OpenMP/Serial
> 均通过；AnalyzeCovariates 当前 pinned CSV 已 bit-identical；本轮另补齐 GATK
> 4.6.2.0 的小写 `-bqsr` short alias（native 原仅接受 `-BQSR`）。专用 alias
> oracle 固定 Java lower/long CSV 与 native lower/legacy-upper CSV 等价，并验证
> `-bqsr=...` 被拒绝；OpenMP/Serial 双 backend CTest 均通过。该项只收敛 CLI parser
> 边界，不扩展为完整 AnalyzeCovariates R 绘图或全量 report 数值等价。

> 2026-09-01 CombineGVCFs interval/reference-block boundary：`-L` 在 REF/<NON_REF>
> block 中间结束时按 GATK 丢弃该 block，END 被区间包含时保留完整 block；aggregate
> 与 `--stream-merge` 共享 predicate。pinned GATK 4.6.2.0 oracle、OpenMP/Serial
> CTest 均通过。

> 2026-09-01 ReblockGVCF overlap block boundary：已确认并修复 native 原先会让
> concrete variant 穿过 hom-ref block 的覆盖错误。native 现在在 block merge 前按
> Java `ReblockingGVCFBlockCombiner` 的 trim/split 逻辑裁剪或拆分重叠 block，拆分后
> 从 indexed reference 更新 REF/POS；无 reference 时需要移动 block 会 fail-closed。
> 新增 pinned GATK 4.6.2.0 synthetic oracle/CTest，固定 Java/native 输出为
> `69000-69004 / 69005 / 69006-69010`，OpenMP/Serial 均通过，并写入
> `overlapping_ref_block_trimmed/split/dropped` telemetry。该项只收敛 overlap
> coverage，不宣称完整 ReblockGVCF annotation/posterior/deletion parity。

> 2026-09-01 VariantRecalibrator full-covariance 默认模型：GATK 4.6.2.0 的
> `MultivariateGaussian` 始终保存并评估完整协方差矩阵；native 已将默认模型从
> diagonal 修正为 full-covariance，`--full-covariance` 仅保留为旧调用兼容拼写。
> pinned VBEM score（最大 VQSLOD 绝对差约 4e-5）、Java model-artifact load 与
> broad contract 在 OpenMP/Serial 双 backend 均通过。该项不宣称完整多高斯 VBEM
> 的 raw-bit 收敛、完整 AS provenance 或大样本模型等价；这些仍是显式 fallback。

> 2026-09-01 VariantRecalibrator zero-variance normalization boundary：pinned
> GATK 4.6.2.0 在训练 annotation 的标准差小于 `1e-5` 时于
> `VariantDataManager.normalizeData` fail-closed；native 原先替换为 unit scale
> 并继续建模。现已在真实归一化路径加入同一阈值守门，constant QD/MQ fixture 的
> Java/native 均以 zero-variance 错误退出；新增 `fastgatk-variant-recalibrator-zero-variance-gatk-oracle`
> CTest #106，OpenMP/Serial 均通过。该项只收敛输入模型有效性边界，不宣称完整
> VQSR raw-bit VBEM 或资源校准等价。

> 2026-09-01 BaseRecalibrator `--mismatches-context-size` 边界：GATK
> 4.6.2.0 的 ContextCovariate 允许 1..13 的 mismatch context；native 之前把
> substitution context 固定为 2。现在参数会进入 Host covariate construction、
> GATKReport Arguments、checkpoint signature 与 manifest telemetry。pinned
> Java/native read-filter oracle 增加 size=3 场景，报告表逐行一致，OpenMP/Serial
> 均通过；完整 BQSR 模型/大样本/cloud 语义仍保持 explicit fallback。

> 2026-09-01 GatherPileupSummaries all-empty writer boundary：GATK 4.6.2.0
> 会在移除所有空 shard 后允许冲突的空表 SAMPLE metadata，并输出仅含列名的
> table header；native 已同步该 writer 与“只对非空输入校验 sample”的顺序。
> `verify_gather_pileup_gatk_oracle.py` 新增 pinned Java/native 逐字节检查，
> OpenMP/Serial contract 与 oracle 均通过。该项只收敛空 shard/metadata 边界，
> cloud 与完整 Java dictionary edge semantics 仍保持 explicit fallback。

> 2026-09-01 VariantFiltration mask reverse boundary：native 新增 GATK
> `--filter-not-in-mask` 与 `--mask-description`。正常模式过滤 mask overlap，
> reverse 模式过滤 mask 外记录，均保留既有 FILTER；默认/自定义 mask FILTER
> header description 与 manifest metadata 已对齐。pinned GATK 4.6.2.0
> `IndexFeatureFile` mask oracle 比较 normal/reverse/chained history/custom
> description，OpenMP/Serial 双 backend 均通过；完整 JEXL、cloud mask、复杂
> allele-specific semantics 仍 explicit fallback。

> 2026-09-01 VariantRecalibrator `--sample-every-Nth-variant` 增量：native 现与
> GATK 4.6.2.0 一样按零基遍历计数保留第 1、1+N、1+2N 条记录，并在 writer
> 边界真正省略未采样记录（此前仅不计算分数却仍输出全量输入）。新增
> `verify_variant_recalibrator_sample_every_gatk_oracle.py`，固定 6 条输入在
> `N=2` 时输出坐标 `[1,3,5]`、8 列 sites-only shape、VQSLOD/训练标签和
> manifest sample/count；OpenMP/Serial 均由 pinned Java oracle 验证。该项收敛
> scatter/recalibration-table sampling provenance，不把完整 full-covariance VBEM
> 或全量 VQSR 模型宣称为 bit-identical。

> 2026-09-01 AnalyzeCovariates 参数比较边界：native 现按 GATK 4.6.2.0
> `RecalibrationArgumentCollection.compareReportArguments()` 的固定字段比较
> before/after report；`indels_context_size`、`covariate` 和 report 路径等 Java
> 明确忽略的字段不再误触发失败，旧 report 缺失字段使用 Java 默认值。新增
> `indels_context_size` 差异 oracle，OpenMP/Serial 均 CSV bit-identical；PDF
> 仍只保证确定性语义，不宣称 R 绘图字节一致。

> 2026-09-01 ValidateVariants 参数兼容增量：补齐 GATK/Barclay 的 `-D`、`-isr`、
> `-do-not-validate-filtered-records`、`-warn-on-errors` 与
> `-disable-sequence-dictionary-validation` 单短横线别名；pinned GATK 4.6.2.0
> command oracle 与 native OpenMP/Serial 均验证退出码一致，invalid Boolean 与互斥
> 组合仍 fail-closed。该项只收敛 CLI 兼容性，不改变 ValidateVariants 未覆盖的
> Tribble/htsjdk、dbSNP overlap 与 dictionary-wide GVCF 语义。

> 最新复核（2026-09-01）：双后端当前 CTest 清单各 209 项。最新增量加入
> VariantRecalibrator sampling provenance、GenotypeGVCFs GP/PG 输入兼容与 ModelSegments allele-fraction initializer 三项
> pinned oracle；随后加入 ApplyVQSR default-cutoff、GenomicsDBImport sample-map
> 与 CallCopyRatioSegments interval-validation oracle；前一批 PoN v7
> sample-metadata 与 HC RCM PL-range 也已在 OpenMP/Serial 双 backend 通过。当前
> 清单实际为 209 项，新增边界共 63 项，仍以
> 140/140 作为连续 CPU 基线。
> 随后新增 VariantEval `--keep-ac0` aggregate、GenomicsDBImport
> `--genomicsdb-update-workspace-path` 与 VariantRecalibrator retry/iteration
> 参数边界，三项 pinned oracle 均已在 OpenMP/Serial 双 backend 通过。
> dense-reference INFO/DP writer 边界后，按同一版本分段执行的当前完整清单合计
> OpenMP/Serial 均 140/140：#1–#59 为 59/59，#60–#67 为 8/8，#68–#140 为
> 73/73（总耗时分别 2213.42 s、2266.39 s）。随后新增 CombineGVCFs PL-less、
> FilterMutectCalls orientation-joint 与 CreateReadCountPanelOfNormals 退化 SVD oracle，
> OpenMP/Serial 均通过。此前 134/134 是修复前历史矩阵。

> FilterMutectCalls normal-artifact、DenoiseReadCounts integer-input、HC AssemblyRegion
> boundary 与 Mutect2 global-mismapping-rate oracle 也已在双 backend 定向通过。
> BQSR preserve、BQSR read-filter、HC k-mer-list、FilterMutectCalls orientation、SelectVariants FILTER、CallCopyRatioSegments 与 GenotypeGVCFs dense
> writer 切片均已双 backend 定向通过；此前新增 BQSR BI/BD indel、
> SelectVariants filtered-GT、CollectAllelicCounts read-filter，以及此前的 BQSR CRAM、ApplyVQSR、
> SortSam、HC min-pruning、HC soft-clip、FilterMutectCalls GermlineFilter 与
> VariantRecalibrator model bounded oracle。

> 2026-09-01 全量复核：VariantRecalibrator 输入模型顺序 contract 修正后，OpenMP
> 与 Serial 各自独立连续运行此前 194 项 CTest，均为 194/194 通过（总耗时
> 1198.16 s / 1184.40 s）；随后新增 DepthOfCoverage deletion-site 与 MarkDuplicates pair-key oracles，双后端专测通过。
> 该结果证明当前配置没有已知测试回归，但不改变各工具
> 完整算法、bit-identical、GPU 或多节点生产语义仍部分为 fallback 的结论。

> 2026-09-01 GenomicsDBImport sample-name-map bounded oracle：native sparse workspace
> 严格校验两列/唯一 sample 名，并在私有 materialized VCF/GVCF 副本中完成样本名
> 重写；源输入不变，native GenotypeGVCFs header 与 pinned GATK 4.6.2.0
> `callset.json` sample name 一致。OpenMP/Serial 均通过；该项只覆盖 sample-map
> 与重试物化语义，opaque TileDB/完整 GenomicsDB storage 仍 fallback。

> 2026-09-01 GenomicsDBImport incremental workspace bounded oracle：外部 adapter
> 透传 `--genomicsdb-update-workspace-path`，要求已存在 workspace，保留完整
> `fastgatk-inputs.tsv` 索引并复现 GATK 4.6.2.0 update 时沿用初始 interval set；
> native sparse workspace 对本地输入执行 deterministic atomic input-index rebuild，
> 保留旧/新输入顺序与 sample-map，manifest 标记 non-TileDB。OpenMP/Serial CTest
> 均通过；完整 TileDB append/storage 仍由 GATK backend 负责。

> 2026-09-01 GenomicsDBImport native interval-span bounded oracle：修复 native sparse
> workspace 忽略 `-L/--intervals` 的真实缺口。现在按 contig 名称和完整记录 span
> 执行 GATK UNION/INTERSECTION，跨 selector 的 REF confidence block 保留完整，
> 过滤/样本名重写后不复制 stale VCF/CSI index，checkpoint 固定初始 selector。
> pinned GATK 4.6.2.0 bridge 导出的 `20:69491-69500` block 与 `20:69511-69511`
> variant 两例均 rows exact，OpenMP/Serial CTest 均通过；完整 TileDB storage/query
> 仍为 external backend/bridge 边界。

2026-09-01 CallCopyRatioSegments interval validation 增量：按 sequence dictionary
顺序拒绝 contig 倒退，并拒绝同 contig 重叠 segment；pinned GATK 4.6.2.0
`verify_call_copy_ratio_segments_interval_validation_gatk_oracle.py` 在 OpenMP/Serial
均通过。该项只收敛输入区间校验，不宣称完整 CNV caller 等价。

2026-09-01 HC bounded region contract 更正：多 tile/低内存二分会独立重建
activity/graph context，因此不能把局部 DP、候选集合或 VCF 证据与 aggregate 宣称
byte-identical；contract 现验证坐标排序、唯一性、区间归属和 backpressure telemetry。

本文件是一次只读、可复核的仓库审计，目的不是重新定义模块计划。它把“代码已
存在”“契约测试通过”“可以直接替换 GATK”三个容易混淆的概念分开。当前主计划
仍是 [FAST_GATK_EXECUTION_PLAN.md](</home/turing-agents/Documents/fast-gatk/FAST_GATK_EXECUTION_PLAN.md>)；
[MODULE_IMPLEMENTATION_PLANS.md](</home/turing-agents/Documents/fast-gatk/MODULE_IMPLEMENTATION_PLANS.md>)
和 [EXECUTION_PLAN_MODULES_DRAFT.md](</home/turing-agents/Documents/fast-gatk/EXECUTION_PLAN_MODULES_DRAFT.md>)
保留了历史增量/任务分解，不能把其中较早的状态段落当成最新快照。

2026-09-01 ReblockGVCF 增量：修正 `--rgq-threshold` 从整数解析为 GATK 兼容的
有限非负 double；`10.5` fractional threshold 的 Java/native 完整输出 oracle
在 OpenMP/Serial 均通过。此边界不宣称完整 reblocking 算法等价。

- 2026-09-01 HC RCM PL range bounded oracle：共享 Kokkos
  `ReferenceConfidenceModel` Number=G materializer 不再套用 PairHMM candidate
  envelope 的 PL=999 上限，Host indel RCM 分支也仅受 signed-int 表示范围约束。
  pinned `fastgatk-hc-rcm-pl-range-gatk-oracle` 使用
  `--indel-size-to-eliminate-in-ref-model 1` 的高深度 BP_RESOLUTION
  `17:69646-69647`，Java 4.6.2.0 与 native 两行逐字段 exact，最大 PL=1080；
  OpenMP/Serial CTest #50 均通过。该项只收敛 RCM PL 表示边界，不宣称完整 HC
  assembly、跨 region 或 GVCF 全文件 bit identity。

BQSR Gather 的当前独立边界证据为
`verify_gather_bqsr_gatk_oracle.py`：两个 disjoint Java BaseRecalibrator 报告以
`quantizing_levels=4` 合并时，native 与 GATK 4.6.2.0 的 Arguments、Quantized、
RecalTable0/1/2 五张表逐行一致；indel-context 维度不一致时两者均 fail-closed。
该证据已在 OpenMP/Serial 双 backend 通过，不扩展为 BQSR 大样本、cloud 或完整模型
等价声明。

本次 BQSR I/D 边界新增
`fastgatk-bqsr-indel-gatk-oracle`：在 pinned GATK 4.6.2.0 同一批 reads 上注入
确定性的 per-base `BI`/`BD` Z tag，native 与 Java 的 `Quantized`、
`RecalTable0/1/2` 均逐表一致，并固定 I={10,20,30,40}、D={11,21,31,41} 的
quality-keyed rows。缺失 tag 使用 GATK Q45/default fallback；非空截短 tag 被
native Host reader fail-closed。OpenMP/Serial 各 `1/1` 通过；该证据只覆盖 I/D
event-table/report serialization slice，不扩展为完整 BQSR 模型等价。

本次补充 ApplyBQSR preserve-threshold 边界：`fastgatk-bqsr-preserve-gatk-oracle`
使用 pinned GATK 4.6.2.0 report 和 Q0..Q40 混合质量的同一输入，验证
`--preserve-qscores-less-than 10` 下 Java/native 493 条记录逐行一致，并单独
断言全部低于 Q10 的碱基质量未改变。OpenMP/Serial 各 `1/1` 通过；该 bounded
证据不扩展为完整 BQSR 模型、大样本或 cloud 等价。

本次新增的 HaplotypeCaller k-mer-list 边界为
`fastgatk-hc-kmer-list-gatk-oracle`：native 按 GATK 4.6.2.0 支持重复
`--kmer-size`，默认 `[10,25]`，并在 Host 规范化为排序去重列表；共享 Kokkos graph
orchestration 按该列表尝试后再执行 cyclic/non-unique retry。默认与乱序 `[25,10]`
的 VCF 在 OpenMP/Serial 均与 Java 4.6.2.0 字节一致，manifest 同时固定 normalized
list 与 legacy scalar 语义。该证据只覆盖 ReadThreadingAssembler 参数/attempt-order
boundary，不证明完整 assembly graph/path 或全基因组算法等价。

本轮新增 HC AssemblyRegion boundary bounded oracle：GATK 4.6.2.0 的默认
`assembly-region-padding=100` 与 `max-prob-propagation-distance=50` 已同步到 native；
Host 读段裁剪后按原始 `AlignmentUtils.unclippedReadLength`（去除 soft-clip）执行
10 bp gate，保留长读段在小 AssemblyRegion 边界留下的短片段。新增
`verify_hc_assembly_region_boundary_gatk_oracle.py`，在 `--max-assembly-region-size 50`
的 chr17 fixture 上与 Java VCF 完整 byte-identical，17:69368 的 `AD=20,22 DP=42`
也固定一致；OpenMP/Serial 的 `fastgatk-hc-assembly-region-boundary-gatk-oracle`
均通过。该证据只覆盖 active ownership/trim/PairHMM 输入边界，不扩大为完整 HC
assembly/graph 等价声明。

本次新增的 [progress_score.json](</home/turing-agents/Documents/fast-gatk/progress_score.json>)
冻结了 `weighted-gates-v1` 的权重、五项门禁、工作流分数、证据路径和 registry
快照。它把当前专家估计重放为 `global_score=0.5385`、四舍五入 **54%**；分数仍
不是“已完成工具数量”的比例，也不会把未验证的 GPU/多节点能力记为完成。

FilterMutectCalls 的 GATK CLI alias 边界也已固定：pinned 4.6.2.0 help 声明
`--create-output-variant-index,-OVI <Boolean>`，native 与 dispatcher 均支持
`-OVI false`/`-OVI true` separated form；Java/native oracle 对压缩 VCF 的
`.tbi` 缺失/存在逐项一致，OpenMP/Serial 均通过。`-OVI=true` 是 GATK
4.6.2.0 不接受的 inline 变体，保持在契约外。

> **FilterMutectCalls contamination joint-learning 增量（2026-09-01）**：统一的
> Kokkos contamination AD/POPAF batch 现在同时服务经验阈值 pre-pass 和最终
> filtering；Java `ErrorProbabilities` 的 `NON_SOMATIC=max(Germline,Contamination)`
> 合并被固定。新增 `fastgatk-filter-mutect-calls-contamination-joint-oracle`，在
> pinned GATK 4.6.2.0 的 10 条固定 AD 梯度上验证 contamination FILTER 集合，
> Java 阈值为 0.193、native 阈值为 0.2381429495（差异来自其余 release-specific
> somatic calibration），并断言 native 不再退化为 1.0；OpenMP/Serial 均通过。
> 这是一个真实联合后验边界，不等同于完整 FilterMutectCalls bit-identical。

> **FilterMutectCalls orientation joint-learning 增量（2026-09-01）**：启用
> `--orientation-bias-artifact-priors` 时，Kokkos F1R2/F2R1 weighted-median
> `ReadOrientationFilter` 后验同时进入经验 `OPTIMAL_F_SCORE` pre-pass 的
> `ARTIFACT` reduction 和最终 ALT 联合后验。新增
> `fastgatk-filter-mutect-calls-orientation-joint-oracle`，pinned GATK 4.6.2.0
> 验证 orientation FILTER pattern、GATK 阈值 0 和无 priors baseline 的阈值差异；
> OpenMP/Serial 均通过。其余 release-specific calibration 仍为 fallback。

## 结论

主计划顶部的“总体约 54%”现在可以由版本化的 `progress_score.json` 和
`compute_progress_score.py` 重放；但它仍是门禁分数的专家赋值，不是代码行数比例、
工具数量比例或 CTest 通过比例。当前可被命令直接复核的事实如下：

仓库根目录没有统一 `README.md`；`fastgatk-core`、`fastgatk-runtime`、
`fastgatk-kernels` 和 `fastgatk-native` 的 README 主要描述接口/能力，不提供全局
完成百分比。全局百分比的现行机器可读出处是根目录的 `progress_score.json`；README
和主计划只负责解释方法和边界，不能替代该快照或测试输出。

| 指标 | 当前事实 | 能证明什么 | 不能证明什么 |
|---|---:|---|---|
| **计数说明** | **209** | 上方“第三波校正”及现场 `ctest -N` 的 209 项为权威当前值；下方历史行中的 207/208 项仅保留原始快照文字 | 不将测试数量等同于算法完成度 |
| native registry 条目 | 48 | 已为 48 个入口声明 native binary、fallback 和 backend 元数据 | 48 个都已实现 GATK 算法；registry 中 31 个仍明确是 prototype |
| registry 状态 | 31 prototype、1 adapter、16 contract-compatible | 只有 17/48（35.42%）已离开 prototype 状态 | 35.42% 不是总体完成度；`contract-compatible` 也不等于算法等价 |
| native binary 声明 | 48/48 | 每个 registry 条目有可解析的 native target 路径 | 路径存在不等于该工具的完整语义已覆盖 |
| fallback 命令 | 48/48 | 每个入口都有显式 GATK fallback 方向 | fallback 不是 native 完成 |
| required input/output 元数据 | 45/48、39/48 | registry 对主要工具有输入/输出契约描述 | 某些工具天然没有统一的 required input/output；不能直接当失败率 |
| 最新 CTest 配置快照 | OpenMP/Serial 当前各 209 项；此前 194 项独立连续运行 194/194 通过（1198.16 s/1184.40 s），新增 deletion-site、pair-key、first-oriented-ALT、count-backed AF likelihood、HC CIGAR indel-activity、GenomicsDB native interval-span、CallCopyRatio compensated-sum、GenotypeGVCFs STARTS_IN、VariantEval ValidationReport、BQSR context-size、SortSam duplicate comparator、Mutect2 NLOD、PairHMM fragment aggregation、GenomicsDB native-storage-boundary 与 somatic likelihood/NLOD oracles 已双后端专测通过；此前 140/140 为历史基线 | 当前测试清单、回归与新增边界证据可复核 | 完整算法、GPU、多节点仍未完成 |
| CTest 条目（OpenMP） | 当前配置 207 项；此前连续版本基线 140/140（分段 #1–#59、#60–#67、#68–#140），其后新增 HaplotypeCaller floor-blocks/CIGAR indel activity、GenotypeGVCFs orphan spanning-deletion/STARTS_IN、ModelSegments default-kernel/multisample/first-oriented-ALT/count-backed AF likelihood、VariantRecalibrator annotation-order/zero-variance/retry-iteration guards、CallCopyRatioSegments non-finite/compensated-sum、MarkDuplicates OpticalOnly/pair-key、DepthOfCoverage default read-filter/deletion-site、VariantEval keep-ac0/ValidationReport、GenomicsDBImport update-workspace/native interval-span、Mutect2 PairHMM indel-quality/NLOD、FastaAlternateReferenceMaker homozygous-IUPAC、BQSR context-size、SortSam duplicate comparator、workflow-local 与既有各项边界均已定向通过 | 该构建配置的契约、oracle、benchmark、smoke、progress-score、fixture-digest、resource-limit、SLURM/scheduler retry、ValidateVariants/GenomicsDBImport/bridge、LeftAlign/VariantEval、BQSR/HC/Mutect/CNV/Genotype bounded 覆盖面 | 通过条目不代表全工具或全数据集等价 |
| CTest 条目（Serial） | 当前配置 207 项；此前连续版本基线 140/140，其后新增边界均已在第二 execution space 定向通过 | 同一 Kokkos 源码的第二 execution space 有对应回归入口；上述双后端兼容性边界均可重放 | CUDA/HIP/SYCL 或真实 GPU 尚未被此数字证明 |
| 第三波校正 | 现场 `ctest -N` 已更新为 OpenMP/Serial 各 209 项；新增 PairHMM fragment aggregation 与 GenomicsDB native-storage-boundary，ModelSegments provenance 仅修正可观测声明，somatic likelihood/NLOD oracle 复核通过 | 209 项中的新增边界均已定向双后端通过；完整 Gibbs/MCMC、TileDB native storage 仍 no-go/fallback | 不把新增测试数当作算法完成度 |
| 回归套件（本次） | OpenMP/Serial 此前 194/194 连续通过（1198.16 s/1184.40 s），新增 deletion-site 与 pair-key oracles 双后端专测通过；另有 140/140 分段历史基线 | 覆盖既有全部门禁以及本轮 ModelSegments 条件模型、GenotypeGVCFs 多样本 reference-confidence、FilterMutectCalls contamination-joint 等新增边界；通过条目不代表全工具或全数据集等价 | 真实 GPU/多节点仍未作为本轮发布门禁 |
| dispatcher contract（本次） | 142 checks pass | registry JSON、参数声明和 dispatcher dry-run 当前一致 | native tool 的输出生物学等价 |
| Kokkos backend matrix（本次） | OpenMP/Serial pass，checksum 相等；backend-config gate 双 backend 各 1/1 | 统一 Kokkos API、checksum 和 benchmark schema 7 在两种 host backend 一致；显式 CUDA/HIP/SYCL 选择、Host Serial 保留和架构/编译器 provenance 可在 CPU-only CI fail-closed 检查 | 设备后端性能、真实 GPU context 和真实集群吞吐 |

本轮兼容性边界的聚焦回归还覆盖 BQSR Java report/ApplyBQSR oracle（四张报告表与
493 条 Java-`PrintReads` 稳定记录，含 `OQ`/`--use-original-qualities`）、HC
reverse-strand indel 的 reverse-complement/anchored CIGAR/left-normalization、
dispatcher launcher contract（`--java-options`、`--gatk-config-file`、递归 `@args`、
`--` 与未知选项 fail-closed），以及 runtime 的 cgroup/SLURM hard/free host budget、
device memory、scratch/local-SSD、remote-input、byte backpressure 和
`RESOURCE_EXHAUSTED` 语义。远端输入先做原子 staging；普通 VCF/Tribble `.idx` 由共享
writer 在 primary output 校验后原子发布。SelectVariants Java 4.6.2.0 semantic
oracle 覆盖 type/sample/JEXL/concordance/allele-subset，并校验稳定的
GT/AD/PL/GQ/AC/AN/AF payload；`fastgatk-runtime-output-smoke` 覆盖完整
primary/index/manifest 发布、缺索引、临时残留 fail-closed 与失败提交回滚。上述
五十四项聚焦门禁在 OpenMP/Serial 均为 `54/54`；它们不改变未完成的完整算法、真实
GPU 与多节点发布边界。

此前一次重建后的 OpenMP 运行曾得到 `82/87`；失败项均是 HC/Flow/indel synthetic
输入契约（缺少 `@RG/RG` 或端到端 flow 读长不足），并非源码断言被放宽。补齐输入
契约后最近一次完整双 backend `ctest --output-on-failure` 均已通过 **134/134**（OpenMP 2184.97 s、Serial 2146.46 s）；本次新增 ValidateVariants、GenomicsDBImport/bridge、LeftAlignAndTrimVariants、VariantEval、HDF5 metadata、VariantFiltration/MarkDuplicates Java oracles、BQSR CRAM/BI-BD indel、ApplyVQSR、SortSam、HC min-pruning/soft-clip、VariantRecalibrator model、FilterMutectCalls GermlineFilter、VariantsToTable molten/split oracle、SelectVariants filtered-GT、CollectAllelicCounts read-filter、GatherBQSRReports/FilterMutectCalls/CollectReadCounts GATK oracles、scheduler retry、legacy AF-calculator oracle、SLURM wrapper 与
Kokkos backend-config gates 后，GenotypeGVCFs 定向矩阵为双 backend **6/6**。
该结果仍只证明已登记的双 host backend 回归，不覆盖
真实 CUDA/HIP/SYCL GPU 或多节点发布门禁。

本次命令输出中的关键值：OpenMP 与 Serial 的 kernel benchmark 均为
`schema_version=7`、`sw_checksum=335360`、`genotype_checksum=34048`；OpenMP 的
SIMD width 为 8（PairHMM）/16（SW），Serial 为 1/1；float32 PairHMM 的
`max_abs_error=1.32669e-06`，这条路径仍然是显式近似模式，不能作为 GATK strict
raw-bit 证据。

GenotypeGVCFs 的文件边界 benchmark 现在可用
`--records-per-shard N --repetitions N` 扩展 workload；默认 128 records/shard
仍被标记为 `workload_class=micro`，并在 JSON 中固定
`speedup_claim_allowed=false`。该 benchmark 不运行 Java/GATK 对照，因此任何
128-record p50 只能证明 native 回归和流水线 telemetry，不能支持相对 Java 的
加速结论。需要做吞吐比较时，应使用至少 4096 records/shard、相同输入范围和
匹配的 Java 命令，并同时记录 wall time、RSS、I/O 与 kernel telemetry。
ReblockGVCF 文件边界 benchmark 也支持通过 `FASTGATK_REBLOCK_BINARY` 切换
OpenMP/Serial binary，并在 JSON 记录实际 execution space；当前 512-record
synthetic smoke 的 OpenMP p50 约 0.76 s、Serial p50 约 0.0058 s，只能说明
小 workload 的 launch/线程池开销，不能外推真实集群吞吐。

本次还实际执行了 Java/native ReblockGVCF 同口径基线：GATK 4.6.2.0 Java 为
约 4.596 s、RSS 320,452 KB，native OpenMP 为约 0.0037 s、RSS 7,096 KB，
两者均输出 15 条记录且参数/输入 hash 一致；由于 workload 未达到 1,024 条记录，
脚本明确写入 `speedup_claim_allowed=false`，这些数字不能作为生产加速承诺。

ReblockGVCF 的直接替换表面本轮又补齐了 GATK 正式选项
`--format-annotations-to-remove` 与 `--do-qual-score-approximation`；原有
native 拼写仍保留。同时 `--floor-blocks` 的 reference-block `PL`/`MIN_DP`
移除语义已按 Java oracle 对齐，OpenMP/Serial verifier 均通过，dispatcher 仍为
142 checks。`--drop-low-quals` 也已改为在合并前丢弃 GQ0/低于阈值的 block，真实
chr17 Java/native 记录与字段 exact，manifest 暴露 dropped block 计数。
`AS_QUALapprox` 的 `Number=A` 序列也与 GATK 一致保留 REF 空槽（`|alt|...`），
正式 QUAL approximation 的 Java/native oracle 在两个 backend 均 exact。
concrete variant 的 drop-mode 也已覆盖 GATK 的两个分支：低置信 ALT 调用被丢弃，
而 PL-derived GQ≥30 且最佳调用为 hom-ref 的站点投影到 REF/<NON_REF>，保留最佳
ALT 的 PL/GQ、DP/MIN_DP、END 并清理 stale INFO/FORMAT；该高置信 hom-ref 转 block
已加入独立真实 Java/native oracle，OpenMP/Serial 均 exact。

本轮新增的 Mutect2 语义门禁也已重放：OpenMP/Serial pinned GATK 4.6.2.0
oracle 均为 6/6 site-set、GT/DP/AD exact；TLOD 最大绝对差仍为 3.500717。
同名双端 read 的非覆盖 mate 已作为 `fragment_only` 进入分组似然，但不会进入独立
AD/DP/annotation 行；这提高了 fragment-first 兼容性，但尚未消除完整 assembly/
posterior 的 release-specific 差异。

本轮新增的 FilterMutectCalls release-specific 默认参数 guard 也已通过双 backend：
GATK 4.6.2.0 的 `--min-median-mapping-quality` 原始 `-1` sentinel 在正常模式由
getter 解析成有效 30、microbial mode 解析成 20；正确 `MMQ=Number=R` 的
`REF,ALT=60,29` 调用得到 Java/native 一致 `map_qual` 行为。这只证明参数默认值
和该 hard-filter 边界，不能证明完整 somatic ErrorProbabilities joint posterior。

本轮还补齐了 HC/Mutect2 CLI 默认 `WellformedReadFilter` 的两个可验证边界：
`CigarContainsNoNOperator` 与 `HasReadGroup`。两个标志通过统一 Kokkos mask 执行，
并在 telemetry/manifest 暴露；双 backend API smoke、HC/Mutect2 contract 与 GATK
oracle 均通过。共享 kernel API 默认仍保持向后兼容的宽松策略。

本次 Mutect2/FilterMutectCalls 增量补齐 orientation `ROQ` writer 边界：对 pinned
GATK 4.6.2.0 canonical/reverse-complement prior、`F1R2/F2R1 Number=R` 和 MNV
fixture，GATK/native `ROQ=[93,93,1]` 与 orientation FILTER pattern exact。
native 已将零/确定性 orientation artifact posterior 的旧 `1000/0` 编码改为
GATK `QualityUtils.errorProbToQual()` 的 `[1,93]` bounded phred；OpenMP/Serial
定向 oracle 均通过。该证据仅覆盖 quality serialization，不扩展为完整
FilterMutectCalls joint posterior 或 Mutect2 全模型等价。

本轮随后收紧了 OA/XM 的一个容易被 flat-buffer 隐藏的边界：HTSlib Host reader 与
`ReadBatch`/`ReadFilterInput` 现在分别保存 tag presence，Kokkos
`NonChimericOriginalAlignmentReadFilter` 因而能区分“present-but-empty”和“missing”。
新增 API smoke 覆盖空值/缺失值组合；OpenMP/Serial 的 HTS reader、HC、Mutect2、
PairHMM 定向回归和 pinned GATK oracle 均保持通过，未改变现有 6/6 site-set 与
GT/DP/AD exact 结果。

GenotypeGVCFs 随后按证据范围晋级为 `contract-compatible`：本地 VCF/GVCF、
Fast-GATK `gendb://` sparse workspace、索引输出、stream-by-locus、重复区间集合，
以及单样本/多样本/多等位/InbreedingCoeff 四个 GATK 4.6.2.0 oracle 均通过。
`--include-non-variant-sites` oracle 共输出 1,001 个站点，其中 998 个为真正的
REF-only（ALT=`.`）且 GT/DP/RGQ/INFO-DP 全部 exact；其余 3 个是具体变异，不再被
误报为三条 REF-only 差异。opaque TileDB、release-specific cohort calibration、
全基因组 provenance/raw bit identity 仍保持显式 GATK fallback。

GenotypeGVCFs 的真实 Java/native oracle 现在额外运行
`--use-posteriors-to-calculate-qual`（`--gp-qual`）：对 HaplotypeCaller 产生的、
仅含 PL/GQ 而没有 FORMAT/GP 的 gVCF，Java 与 native 三条 variant row 均整行
exact，QUAL 保持默认 cohort 结果，native manifest 的 posterior kernel calls
与 samples 均为 0。这固定了 GATK `hasPosteriors` 条件，避免把 PL 错当 GP posterior
而重写 QUAL；只有显式 posterior-assignment 输入才进入 native posterior kernel。

同一 GenotypeGVCFs oracle 还覆盖了 `--use-new-qual-calculator false` 的 legacy
AF-calculator 边界：GATK 不运行新 EM kernel，但仍会把 HaplotypeCaller 输入中
`MLEAC/MLEAF` 的 Number=A 值从 concrete ALT+`<NON_REF>` 投影到最终 ALT 列表。
native 现在在 ALT union 和移除 `<NON_REF>` 两个阶段均通过 Kokkos
allele-field remap 完成该投影；OpenMP/Serial 各 3 条记录的 Java/native 行与
semantic header 均 exact，且 legacy manifest 的 cohort AF kernel calls 为 0。

ReblockGVCF 的 GATK 约束也已用独立 oracle 固化：Java 只接受同一样本的非重叠
输入 shard，多样本 VCF 在初始化阶段拒绝；新增的
`verify_reblock_gatk_shards.py` 对两个 `.idx` shard 的 block/ALT/PL/AD/GT/DP/GQ
核心字段在 OpenMP/Serial 均 1/1 exact。native 额外的 sample-major 多样本 kernel
contract 仍是 prototype，不计入 GATK 直接替换能力。

同一 oracle 还暴露并锁定了一个 RGQ 边界：native 现在正确解析
`--rgq-threshold-to-no-call VALUE` 的 separated/inline 两种写法，并在不带
`--drop-low-quals` 时也按 PL[0] 阈值转换低质量 deletion 为 GQ0 reference block；
双 backend contract 与真实 GATK shard oracle 均通过。

随后补齐了 deletion annotation 清理：低质量 deletion 进入 GQ0 reference block
时不再保留 stale `FORMAT/AD`，而是写出 `GT=0/0`、`PL=0,0,0`、`GQ=0`、
`DP=MIN_DP` 并保留原 deletion 的 `END` 覆盖范围；该断言已接入现有双 backend
Reblock contract，避免下游 CombineGVCFs/GenotypeGVCFs 继续消费过期等位深度。

本轮并行收口的新增边界已进入当前 122 项矩阵：ValidateVariants 的 Boolean/互斥选项与 strict/warn
fail-closed 由 GATK oracle 覆盖；GenomicsDBImport 的真实 TileDB workspace publication、输入索引、
GATK GenotypeGVCFs reopen 与 legacy bridge export 由双 oracle 覆盖；LeftAlignAndTrimVariants 的多等位/符号 ALT/
任意倍性字段重映射由 contract+GATK oracle 覆盖；VariantEval 的 CountVariants aggregate
与 Sample GATKReport 由 Java oracle 逐字节覆盖；CollectReadCounts 的 HDF5
本次 oracle 进一步固定 CountVariants 的 genotype-aware ignore-AC0：含 FORMAT/GT 的
hom-ref/all-no-call ALT 记录计为 reference loci，sites-only 记录保持 site-level variant，
并覆盖 nCalledLoci、类型计数、singleton 与 insertion/deletion ratio。
本轮新增 `--keep-ac0`/`-keep-ac0` bounded slice：genotyped AC=0 SNP/INDEL 在
CountVariants aggregate 中恢复为 site-level variant；pinned GATK 4.6.2.0
aggregate report 逐字节一致，bare/true/false 与短别名由 native-only checks 固定，
不宣称其它 evaluator/stratifier 表已达 Java parity。
SimpleCountCollection sample/dictionary/interval/count metadata 由 Java round-trip 覆盖；
GatherBQSRReports 按首个 Java
report 的 `Arguments`/`quantizing_levels` 合并并对非默认 quantizing-levels 做五表
GATK oracle；FilterMutectCalls 在完整 AD/F1R2 过滤后支持
`--sites-only-vcf-output` 的 8 列 writer 边界，并固定 MMQ raw `-1` sentinel 到
normal=30/microbial=20 的 getter 语义。CollectReadCounts 默认输出 GATK
HDF5SimpleCountCollection，显式 `--format TSV` 才走文本路径，read-start/filter/
interval-merging 的 Java oracle 通过。随后 `CollectReadCounts` 还补齐了 GATK
`-XL/--exclude-intervals`：native 在计数 target（而非整条 read）上执行 subtract，
输出分裂后的区间并记录 excluded telemetry；pinned GATK 4.6.2.0 oracle 覆盖该输出和
`-imr` 别名。Java CNV validation 对 `ALL` 和非零 `-ip/-ixp` 均拒绝，native 同样
fail-closed。这些边界均使用统一 OutputManifest/telemetry；本 slice 的
OpenMP/Serial 定向 contract 与 pinned oracle 均通过，完整矩阵状态以本审计的最新快照为准。

## 可重放的核对命令

以下命令不修改输入数据；其中 `ctest -N` 只列出测试，不运行完整套件：

```bash
cd /home/turing-agents/Documents/fast-gatk

python3 - <<'PY'
import collections, json
p = "fastgatk-native/dispatcher/tool_registry.json"
d = json.load(open(p, encoding="utf-8"))
entries = d.get("tools", d)
if isinstance(entries, dict):
    entries = list(entries.values())
print("entries", len(entries))
print("status", collections.Counter(e.get("status") for e in entries))
print("native_binary", sum(bool(e.get("native_binary")) for e in entries))
print("fallback_command", sum(bool(e.get("fallback_command")) for e in entries))
PY

third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest \
  --test-dir fastgatk-native/build -N
third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest \
  --test-dir fastgatk-native/build-serial -N

python3 fastgatk-native/dispatcher/verify_dispatcher.py
python3 fastgatk-native/scripts/compute_progress_score.py
python3 fastgatk-native/scripts/verify_fixture_digests.py
python3 fastgatk-native/scripts/verify_kokkos_backend_matrix.py
# 128 records/shard is a micro regression smoke; it is not a speedup claim.
FASTGATK_GENOTYPE_BINARY=fastgatk-native/build/fastgatk-genotype-gvcf \
FASTGATK_HC_BINARY=fastgatk-native/build/fastgatk-hc-call \
python3 fastgatk-native/scripts/benchmark_genotype_gvcf.py \
  --records-per-shard 4096 --repetitions 1
third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest \
  --test-dir fastgatk-native/build --output-on-failure \
  -R 'fastgatk-runtime-smoke|fastgatk-runtime-pipeline-smoke|fastgatk-kernels-api-smoke|fastgatk-flow-pairhmm-gatk-oracle'
third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest \
  --test-dir fastgatk-native/build-serial --output-on-failure \
  -R 'fastgatk-runtime-smoke|fastgatk-runtime-pipeline-smoke|fastgatk-kernels-api-smoke|fastgatk-flow-pairhmm-gatk-oracle'
```

本次实际执行结果为：registry `48 / (31, 1, 16)`；dispatcher `142 checks`；
OpenMP 与 Serial 最近一次完整 CTest 均 `125/125`（OpenMP 766.81 s、Serial 1024.16 s），新增 ValidateVariants、GenomicsDBImport/bridge、LeftAlignAndTrimVariants、VariantEval、HDF5 metadata、VariantFiltration/MarkDuplicates Java oracles、BQSR CRAM/BI-BD indel、ApplyVQSR、SortSam、HC min-pruning/soft-clip、VariantRecalibrator model、FilterMutectCalls GermlineFilter、VariantsToTable molten/split、SelectVariants filtered-GT、CollectAllelicCounts read-filter、GatherBQSRReports/FilterMutectCalls/CollectReadCounts GATK oracles 与 scheduler retry 后
GenotypeGVCFs 定向回归均 `6/6`（含
无 GP 输入的 `--gp-qual` Java/native exact guardrail），其中 `ctest -R fastgatk-reblock-gvcf`
在两个 backend 均 `4/4`，真实 GATK oracle 覆盖默认、floor/drop、低质量变体丢弃、
高置信 hom-ref 转 block 与 QUALapprox；resource-limit gate 均 `1/1`，核心 HC/Flow/indel 定向测试均
`4/4`，进度评分与 CountBases 均 `2/2`，fixture digest 均 `1/1`；backend matrix
`pass` 且 benchmark checksums 相等。文档中的 82/82、83/83、84/84 等仅是历史快照，
不能覆盖当前 125 项清单。

此外直接执行 `fastgatk-native/scripts/verify_pipeline_local.sh` 时，当前环境的
Nextflow 26.04.6 已实际跑通 `nextflow_smoke.nf`、`nextflow_scatter_gather.nf`
和通用 `nextflow_gatk_compat.nf --tool CountReads`；本地 SLURM wrapper、scatter/
gather、resume 和 manifest/index 检查均通过。这证明现有 DSL2/dispatcher 表面可用，
但不替代真实 SLURM allocation、多节点故障重试和大规模吞吐门禁。

## 为什么“54%”目前不可量化复现

主计划按运行时/I/O、Kokkos 核心、P1 主流程、长尾工具、GPU/集群生产化加权，这
比按行数或命令数合理。现在权重、门禁分数和分数快照已经提交在
`progress_score.json`，所以 `0.5385 -> 54%` 可以由 CI 重放；但 fixture digest、
每个工具的独立门禁证据和“完成”是否要求生产发布仍不是该全局快照的自动推导项。
历史文档中出现过 50%、51%、53%、54% 是日期快照，不是可相加的指标。

下一步应把现有工作流评分细化到每个 registry 工具，并继续使用固定五项门禁；取值
只能是 `0/0.25/0.5/0.75/1`：

1. `api_cli`：参数别名、unknown-option fail-closed、registry/fallback 可重放；
2. `oracle`：与 GATK/Picard/HTSlib 的工具级结果对照；严格模式需要 raw-bit/逐
   字段，近似模式必须声明 tolerance；
3. `format_sidecar`：VCF/BAM/CRAM、header、index、BGZF/Tabix、stats/F1R2 等
   sidecar 和失败原子性；
4. `resource_io`：cgroup/SLURM 内存、byte queue、spill、索引 seek、重试和峰值
   RSS/IO 证据；
5. `e2e_perf`：真实文件边界、Nextflow/SLURM 路径和固定输入下的 prepare/kernel/
   encode/total wall、RSS、I/O；

定义模块分数：

```text
module_score = 0.20*api_cli
             + 0.30*oracle
             + 0.15*format_sidecar
             + 0.15*resource_io
             + 0.20*e2e_perf
```

建议的工作流权重（合计 100%）是：runtime/I/O/compat 20%、Kokkos 核心 20%、
HC/Mutect2/BQSR 主流程 30%、联合分型/CNV/VQSR/Picard 等长尾 20%、GPU/真实
SLURM 多节点 10%。全局分数为工作流分数的加权和；任一工具只有五项都有证据且
没有未声明 fallback 时才能把 `oracle=1`，不能用同一条 API smoke 给整个工具包
记满分。当前 `score_is_expert_estimate=true` 应继续保留；该门禁验证算术和 registry
快照一致，不把专家赋值伪装成自动质量结论。

同时建议把 registry 的单一 `status` 拆成四个正交字段：
`cli_compatible`、`format_compatible`、`algorithm_equivalent`、
`production_ready`。现在 HC 的 `contract-compatible` 与其 notes 中“尚未
bit-identical/biologically equivalent”并存，容易被调度器或用户误读为可以无条件
替换 GATK。

## 当前剩余边界

### 已有较强证据、适合复制模式的部分

- `HostBatch -> KernelPlan.prepare -> Kokkos Views -> execute -> collect` 已在
  PairHMM、SW、graph、activity、BQSR、genotype、somatic/reference-confidence
  等 API 中使用；OpenMP/Serial checksum 一致。
- HTSlib CIGAR/flags/tags、bounded reader/pipeline、OutputManifest、VCF/Tabix 和
  多个 GATK 4.6.2.0 fixture 已有定向 oracle。
- PairHMM strict CPU 路径、flow-PairHMM fixture、HC 的若干 ploidy/RCM/annotation
  边界和 GenotypeGVCFs 的若干字段已有证据。

### 仍不能作为“直接替换”的部分

- **HC**：完整 GATK read threading/SeqGraph/AssemblyRegion 跨区 ownership、
  pruning/traversal、soft-clip pileup、完整 haplotype posterior 和更大真实 corpus；
  当前 registry notes 仍明确不是 bit-identical/biologically equivalent。
- **Mutect2/FilterMutectCalls**：release-specific somatic posterior、完整
  read-orientation/contamination/filter learning、joint calibration 和 sidecar
  provenance；当前 TLOD 差异是实质算法边界，不是 writer 格式问题。
- **BQSR**：大规模/长读、known-sites 全边界、object-store/cloud 语义和全量
  ApplyBQSR 端到端压测；小 fixture 表格一致不足以推出生产吞吐。
- **其它 42 个 registry prototype**：即使有 contract/oracle 脚本，也必须逐工具
  完成真实文件、错误/索引/sidecar、资源压力、Nextflow/SLURM 和性能门禁后再升级。
- **GPU/多节点**：当前只证明 Kokkos host execution spaces。CUDA/HIP/SYCL
  correctness、device OOM/backpressure、跨节点 shard/gather/retry 仍没有本次
  审计可以确认的生产证据。

## 下一轮最有价值的工作

按“降低替换风险/提高可测量性”的收益排序：

1. **把进度评分扩展到逐工具并接入发布门禁**：现有快照已经给出工作流分数；下一步
   给 48 个 registry 条目填写五项分数、证据路径、fixture SHA-256、当前 status 和
   fallback reason。这样能把工具级完成度变成可审计数字，也会立刻暴露“有 binary
   但无 oracle”的条目。
2. **做 HC 真正的工具级 oracle corpus**：围绕复杂 indel、跨 AssemblyRegion、
   soft-clip、重复 k-mer、polyploid/multi-ALT 和 gVCF block，逐条比较 records、
   INFO/FORMAT、header/index、manifest；随后才做更大规模 E2E benchmark。当前 HC
   是当前复杂主流程中已被 registry 提升的入口，最值得把“contract-compatible”
   推进为 `algorithm_equivalent` 的候选；CountBasesInReference、CountReads、FlagStat
   和 IndexFeatureFile 的提升主要是固定 reference/read-walker/index 契约，不应与
   HC 的完整算法门槛混为一谈。
3. **补 Mutect2 的 TLOD/posterior parity**：固定 tumor/normal、PON、orientation
   和 contamination corpus，分离 activity、assembly、PairHMM、posterior、writer
   五阶段，避免把 TLOD 差异误归因于 I/O。FilterMutectCalls 随后复用同一证据。
4. **建立真实文件边界性能基线**：运行
   `fastgatk-native/scripts/benchmark_end_to_end.py`，同时保存 native 与 pinned
   GATK 的 decode/prepare/kernel/encode/write、wall、RSS、filesystem I/O、输出
   digest；当前脚本已覆盖 HC、region/contig streaming、Mutect2 和 FilterMutectCalls。
   只有这样才能回答“比 Java 快多少”，而不是只看 kernel p50。
5. **再做硬件矩阵**：先保持同一 fixture 在 Serial/OpenMP 的 checksum 和 strict
   oracle，再在实际 CUDA/HIP/SYCL 设备上运行相同 schema 7 benchmark；显式记录
   execution space、device memory、queue/backpressure 和 fallback，禁止将 host
   模拟结果写成 GPU 证据。
6. **最后做真实 SLURM/Nextflow 压测**：多 shard、resume、失败重试、scratch
   满/内存紧缩、索引 seek 和 gather；收集每个 shard 的 manifest/telemetry，确认
   结果顺序、sidecar 原子提交和 scheduler allocation 与本地 smoke 不同的边界。

## 本审计的边界

本次已重跑完整双 backend CTest（共 125 项），OpenMP/Serial 均 `125/125` 通过；
此外重放了 registry `compute_progress_score`、dispatcher 142 checks、双 backend
核心定向门禁、fixture digest 1/1、resource-limit 1/1，以及 backend matrix。
文档中的 82/82、83、85、86、88、90、91、93、96、98、101、103、110、113/113
等数字只是历史基线。没有因为测试条目多就上调算法完成度，也没有因为 48 个
native binary 路径存在就取消显式 Java fallback。

2026-09-01 bounded VQSR slice：ApplyVQSR 新增真实 GATK 4.6.2.0 AS multi-ALT
语义 oracle。Java one-ALT recal records 的 scalar `culprit` 现在按
contig/POS/REF/ALT join，并投影为输入记录的 Number=A
`AS_VQSLOD`/`AS_FilterStatus`/`AS_culprit`；site FILTER 继续按 applicable ALT
的最宽松结果决定。`verify_apply_vqsr_gatk_oracle.py`、CTest
`fastgatk-apply-vqsr-gatk-oracle` 和 `verify_all.sh` 已接入，OpenMP/Serial
定向运行均 pass。此项只收敛 VCF/AS writer 与 per-ALT join 边界，不改变
ApplyVQSR/VariantRecalibrator 的 prototype 分数；完整 VQSR model、
recalibration-table、cloud/tranche provenance 仍显式 fallback。

2026-09-01 bounded VQSR model-artifact slice：VariantRecalibrator Host 现可读取并
按 annotation 名称重排 pinned GATK 4.6.2.0 GATKReport v1.1 的均值/标准差、
正负 Gaussian PMix、means/covariance 表，复用 Kokkos scoring；pinned Java
oracle 比较 VQSLOD、坐标/记录数、稳定的正资源标签和模型表维度；GATK 模型侧
负位点选择未序列化，因此作为非 bit-identical metadata 边界。scalar recal VCF
writer 对齐 GATK 的 `N/<VQSR>`、`END`、小写 `culprit`、空 QUAL/FILTER 与
`--sites-only-vcf-output` 8-column site-only 形状。新增
`fastgatk-variant-recalibrator-gatk-model-oracle` CTest，并在 OpenMP/Serial
均通过；不改变总进度分数。Java raw-bit VBEM、完整 recalibration-table/
provenance 与完整 AS model parity 仍显式 fallback。

2026-09-01 bounded FilterMutectCalls GermlineFilter slice：native 新增
`GermlineFilter.germlineProbability` 的 Kokkos 后验路径，Host 从 NLOD/POPAF 与
tumor AD 取证据，并以 GATK `QualityUtils.errorProbToQual` 的 [1,93] 边界输出
Number=A `GERMQ`。POPAF 缺失时不伪造注释（fail-closed），可选 `PGERMLINE` 不被
冒充 GERMQ；极低/极高 population-AF 两个 corner 与 interior posterior
materialization 由 pinned GATK 4.6.2.0 oracle 在 OpenMP/Serial 双 backend 验证。
新增 CTest `fastgatk-filter-mutect-calls-germline-oracle` 并接入 `verify_all.sh` 与
dispatcher notes；此项只覆盖可审计 posterior/output 边界，完整 release-calibrated
joint FilterMutectCalls learning 仍显式 fallback，不改变总进度分数。

2026-09-01 并行收口的另外三个边界：GenomicsDBImport native sparse workspace
保留显式 `-V`/sample-map 输入顺序并稳定去重，重开/chr2-only 查询验证 sample
header 与 missing genotype；HaplotypeCaller `--dont-use-soft-clipped-bases` 只
屏蔽 terminal/internal soft-clip 派生 evidence，四个位点 pinned GATK oracle 的
core FORMAT/稳定 annotation 在双 backend 通过；CombineGVCFs 新增
`--sites-only-vcf-output`，在完整 allele/DP/FORMAT merge 后以独立 writer header
和 `bcf_subset_format` 发布 8 列 site-only VCF，reference block 不额外提升 DP，
GATK 4.6.2.0 shape/payload oracle 双 backend exact。三项均不改变 native sparse
非 TileDB、HC 完整 assembly 与 CombineGVCFs joint-genotyping 仍显式 fallback 的范围。

2026-09-01 最后一轮并行边界：VariantsToTable 的 `--moltenize`、`-ASF/-ASGF`
顺序消费、biallelic `Number=R` 投影与 `-EMD` presence 语义均由 pinned GATK
4.6.2.0 oracle 固定，OpenMP/Serial 双 backend 通过；CalculateContamination
对所有位点被 `MIN_COVERAGE` 过滤的合法输入输出 `0.0/1.0` 与 header-only
segmentation，并支持 `-matched/-segments`；GatherTranches 强制 `--mode`、支持
`-tranche`，version-6 header-only scatter shard 按零行 contributor 合并。由于
新增 VariantsToTable 专用 oracle，当前 CTest 清单为 129 项；新增切片仍不改变
registry prototype/adapter 状态，也不代表完整算法 bit-identical。
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

2026-09-01 最新双 backend 回归：BQSR 的 BI/BD per-base indel quality report
（`Quantized`、`RecalTable0/1/2`、截短 tag fail-closed）、SelectVariants 的
`--set-filtered-gt-to-nocall`（FT 过滤后 allele compaction 和 AC/AN/AF 重算），
以及 CollectAllelicCounts 的 GATK 默认 read-filter 集合和
`-DF/-RF/--disable-tool-default-read-filters` 顺序均已由 pinned GATK 4.6.2.0
oracle 固定，OpenMP/Serial 双 backend 通过。完整 CTest 清单 125/125 通过
（OpenMP 766.81 s、Serial 1024.16 s）。这些证据仍只覆盖声明的边界，registry
状态与 54% 加权完成度不变；完整 BQSR、SelectVariants、CNV 模型、真实 GPU/多节点
和云端生产语义仍未达到直接替换门槛。

2026-09-01 GetPileupSummaries inverted filter slice：GATK 4.6.2.0 help 提供
`-XRF/--inverted-read-filter`，而 native 之前只将长名 fail-closed、短名视为 unknown。
现在针对当前 14 个 LocusWalker built-in predicates 完整执行 `InvertedReadFilter`，且顺序为
默认/普通 `-RF` mask 后取反；不把 inverse 解释为禁用 default。采用 Java 可接受的
`--disable-tool-default-read-filters true -XRF MappingQualityReadFilter
--minimum-mapping-quality 50` low-MAPQ table fixture，OpenMP/Serial 均与 Java 输出逐字节
一致，manifest telemetry、registry、oracle/verify_all 已同步。该项只收敛明确的 read-filter
I/O 边界，完整 LocusWalker/cloud 语义仍保持 fallback。

2026-09-01 CountReads/FlagStat inverted filter slice：pinned GATK 4.6.2.0 的
`FlagStat -XRF MappingQualityReadFilter --minimum-mapping-quality 30
--maximum-mapping-quality 50` 在 chr17 fixture 输出 475 条记录及 12 个确定 counters，
而此前 native 对 `-XRF` 返回 unknown option。共享 `read_metrics_tool` 现将
`-XRF/--inverted-read-filter` 映射为已支持 predicate 的 GATK `InvertedReadFilter`
语义，仍先应用默认和普通 `-RF` filters；参数化 filter 的依赖检查也识别该反向形式。
OpenMP/Serial direct Java boundary 均 exact，existing `verify_flag_stat.py`、dispatcher
registry、manifest telemetry 及 verify_all 接入均已更新；这只扩大成熟的 ReadWalker
filter surface，不改变完整 read-filter/plugin、cloud 或全 GATK release parity 的 fallback。

2026-09-01 BaseRecalibrator read-filter slice：native 真实实现 GATK 4.6.2.0
七个默认 read-filter 的独立 predicate，并按 Java plugin descriptor 的
“默认集合（或 `--disable-tool-default-read-filters` 清空）→`-DF` 移除→`-RF` 添加”
顺序解析；未知及 inverted filter 名称 fail-closed。新增
`fastgatk-bqsr-read-filter-gatk-oracle`，以含 MAPQ 0/255、secondary、duplicate、
vendor-fail、supplementary 的合成 SAM，在 default、禁用 MAPQ-zero、清空默认后显式
启用 MQ-zero、显式过滤 supplementary 四个场景比较 GATK 4.6.2.0 的
`Quantized`/`RecalTable0/1/2`，双 backend 各 `1/1` 通过（有效 read 2/3/6/1）。
该 bounded 证据只覆盖 BQSR filter/plugin 参数边界，不改变 54% 全局完成度，完整
BQSR 模型、大样本、cloud 与真实 GPU/多节点仍未完成。

2026-09-01 ApplyBQSR report alias slice：pinned GATK 4.6.2.0 help 明确
`-bqsr` 是 `--bqsr-recal-file` 的短别名；此前 native 仅接受长拼写并对 `-bqsr`
返回 unknown option。native 现已接入同一 report parser，新增
`fastgatk-apply-bqsr-alias-gatk-oracle`，用 Java BaseRecalibrator report 比较
Java/native 长短拼写的 493 条 `PrintReads` 记录；OpenMP/Serial 各 `1/1` 通过。
该 bounded slice 只证明 CLI/report 选择边界，不改变 54% 加权完成度，也不宣称完整
BQSR 模型、cloud 或 MD5 等价。

2026-09-01 CallCopyRatioSegments degenerate-statistics slice：native 现在保留 GATK
4.6.2.0 SimpleCopyRatioCaller 对空/单例 copy-neutral 集合的 IEEE-754 语义：统计
结果为 NaN，非 neutral segments 的比较全部落入 neutral；不再将单例标准差误设
为 0，也不把 outlier 过滤为空回退到原集合。`-O` 默认自动生成 `<output stem>.igv.seg`，
`--legacy-output` 仍可覆盖路径；manifest 对非有限统计写 JSON null。pinned GATK
called/legacy bytes oracle 在 OpenMP/Serial 均通过；full ModelSegments/Cloud semantics
remain explicit fallback。

2026-09-01 CallCopyRatioSegments non-finite input slice：GATK 4.6.2.0 的
`CopyRatioSegment` TSV decoder 接受 `NaN`、`Infinity`、`-Infinity` 的
`MEAN_LOG2_COPY_RATIO`，并由 SimpleCopyRatioCaller 通过 IEEE-754 `pow`/统计路径
传播，输出 neutral call，同时保留非有限值文本。native 现不再提前拒绝这些行，
called 与 `.igv.seg` 均与 pinned Java bytes exact；新增
`fastgatk-call-copy-ratio-segments-nonfinite-gatk-oracle`，OpenMP/Serial 均通过，
manifest 的 NaN 统计保持 JSON null。该 bounded 输入/计算边界不代表完整
ModelSegments 概率分段或 cloud parity。

2026-09-01 CallCopyRatioSegments compensated-sum slice：native 的长度加权均值与
方差累加现在复现 Java 17 `DoubleStream.sum()` 的 `computeFinalSum`，即返回高位
累加器减去负补偿字，并保留同号 infinity fallback。一个只使用精确 2 的幂 copy
ratio 的 pinned 边界夹具证明旧实现会把 GATK 的 neutral call 翻成 amplification；
`fastgatk-call-copy-ratio-segments-compensated-sum-gatk-oracle` 对 called/IGV 两份输出
执行 byte-exact 比较并在 OpenMP/Serial 通过。该 bounded 算术修复不扩大完整
ModelSegments posterior、cloud 或 bit-identical sampler 声明。

新增 SelectVariants FILTER exclusion bounded oracle：`--exclude-filtered` 与
`--exclude-filtered-variants` 支持 GATK-compatible optional Boolean，非-PASS
FILTER 在 true 模式排除，PASS/`.` 保留，false 模式全保留且多标签按 GATK 规范
排序；OutputManifest 提供选择状态与过滤计数。pinned GATK 4.6.2.0 oracle、原有
SelectVariants contract 均在 OpenMP/Serial 通过；registry 仍为 prototype，完整
ID/pedigree/random/cloud/JEXL 语义未完成。

2026-09-01 GenotypeGVCFs dense-reference writer slice：native 在
`--gatk-compatible-annotations --include-non-variant-sites` 下对单态参考位点执行
GATK cleanup（MIN_DP→DP、GQ→RGQ、移除 PL/MIN_DP），aggregate 与
`--stream-by-locus` 均与 pinned GATK 4.6.2.0 的 101 条数据行逐行 exact；OpenMP/Serial
定向 oracle 均通过。CTest 清单随后为 129 项；默认 native 诊断档仍保留 RCQ/RCP/GQ/PL，
因此该证据只覆盖兼容 writer profile，不代表完整联合分型 bit-identical。
### GetPileupSummaries interval-exclusion boundary (2026-09-01)

Implemented and oracle-verified GATK `-XL/--exclude-intervals`, `-ip`, and `-ixp` semantics for native GetPileupSummaries, including padded exclusion and separate-token Boolean default-filter forms. OpenMP/Serial outputs match pinned GATK 4.6.2.0 byte-for-byte for the covered fixtures; broader LocusWalker/cloud semantics remain fallback.

### ModelSegments supplied partition boundary (2026-09-01)

Native `ModelSegments --segments` now consumes a validated Picard interval
list, bypasses kernel segmentation, and preserves supplied segment coordinates.
Copy-ratio points are assigned by GATK's midpoint rule. The new pinned GATK
4.6.2.0 oracle verifies exact boundaries, per-segment point counts, and `.cr.seg`
rows for OpenMP and Serial. This improves workflow compatibility for a prior
multisample segmentation hand-off; full release-specific MCMC/model/report
parity remains explicit fallback.

### DenoiseReadCounts HDF5 PoN interval-identity boundary (2026-09-01)

Native DenoiseReadCounts now rejects a case whose SimpleCountCollection
interval list is not exactly equal to the HDF5 PoN original interval list
(contig/start/end and order), matching GATK 4.6.2.0 before filtered-panel
subsetting. The new pinned Java oracle covers valid acceptance and
subset/superset/coordinate-mismatch rejection; native and Java agree on all
cases in both OpenMP and Serial. This closes a high-impact PoN denominator
compatibility boundary without claiming complete CNV model or Spark/cloud
parity.

2026-09-01 GenotypeGVCFs output-allele subset slice：native 在共享 Kokkos cohort
AF posterior 后应用 GATK GenotypingEngine 的 per-ALT standard-confidence pruning，
默认阈值 30（`--standard-min-confidence-threshold-for-calling` / `--stand-call-conf`），
并在 writer 前重映射 Number=A/Number=G 的 annotation、PL、AD、GT。新增
`fastgatk-genotype-gvcf-assignment-gatk-oracle`，两样本三 ALT unsupported-sibling
fixture 与 pinned GATK 4.6.2.0 的 `PREFER_PLS`、`USE_PLS_TO_ASSIGN` 均整行 exact；
OpenMP/Serial 各 1/1 通过。该切片只证明 output-allele subset 边界，完整 cohort
calibration 与 raw bit identity 仍为显式边界。

同一 `fastgatk-genotype-gvcf-assignment-gatk-oracle` 还固定了 GATK
`GenotypeGVCFsEngine.createMinimalArgs` 的继承参数行为：即使命令行传入
`--genotype-assignment-method`/`--gam`，Java 最终仍强制 `PREFER_PLS`。native
`--gatk-compatible-annotations` profile 不再额外写 GP/PG 或改变 GQ；诊断
profile 的 richer Kokkos assignment 仍仅用于实验回归。该 alias/row 边界在
OpenMP/Serial 均通过，不能外推为完整 cohort bit identity。

VariantFiltration 新增 `--invalidate-previous-filters` optional-Boolean：只清除
site FILTER、保留 genotype FT，再执行新规则；pinned GATK 4.6.2.0 Java oracle
在两个 backend 均通过。

### 2026-09-01 完整双 backend 回归收口

使用 bundled CTest 对当前 134 项清单执行 `--output-on-failure`：OpenMP
`134/134`（2184.97 s），Serial `134/134`（2146.46 s），零失败。该结果覆盖本轮
新增 GenotypeGVCFs assignment、VariantFiltration invalidation、DenoiseReadCounts
PoN interval-identity 以及既有全部 contract/oracle/smoke/benchmark；它证明当前
源码在两个 CPU execution space 的回归闭合，不改变 54% 的算法加权完成度，也不替代
真实 CUDA/HIP/SYCL GPU、多节点 SLURM/Nextflow 压测和完整 Java 生物学模型等价。

### 2026-09-01 ModelSegments supplied-partition 增量

ModelSegments 现真正执行 GATK advanced `--segments`：Picard interval-list
按 sequence dictionary 顺序且不得重叠，copy-ratio 以 interval midpoint 归属，
并跳过 kernel segmentation。pinned GATK 4.6.2.0
`fastgatk-model-segments-input-segments-gatk-oracle` 在 OpenMP/Serial 均验证
边界、点计数和 `.cr.seg` 行一致；native 另补齐四个 `.modelBegin/.modelFinal`
`.cr.param/.af.param` GATK-shaped 参数 sidecar。该增量闭合直接替换的文件/分区
契约，但参数值仍是 bounded deterministic summary，不提升完整 CNV MCMC 等价。
标准 MCMC/smoothing 参数名现可解析并记录，但 Java 的多轮 credible-interval
smoothing 仍未执行，manifest 以 `smoothing_applied=false` 明确该边界。

### 2026-09-01 ModelSegments CopyRatioModeller 条件模型增量

native probabilistic ModelSegments 现对原始 copy-ratio 点执行一个可复现的
Gibbs-conditional responsibility prepass：按照 GATK `CopyRatioSegmentedData`
的 interval-midpoint 归属，将点放入 Kokkos Views；`RangePolicy` 计算点级
outlier responsibility、每段条件均值、全局 variance 与 Beta(5,95) outlier
probability。固定顺序 Host reduction 只汇总全局标量，条件均值与参数用于
bounded chain 的初始化/步长，并写入 manifest 与 `.modelBegin/.modelFinal.cr.param`
的全局摘要。新增
`fastgatk-model-segments-copy-ratio-conditionals-gatk-oracle`，在含离群点的
fixture 上同时运行 pinned GATK 4.6.2.0 与 native，验证 partition、point count、
GATK-shaped parameter schema、有限参数和主模态保持；OpenMP/Serial 均通过，既有
ModelSegments contract 也保持通过。该切片不宣称 Java slice-sampler 随机抽样、
latent-indicator MCMC、credible-interval smoothing 或 posterior raw-bit identity。

同一 fixture 现附带真实 `REF_COUNT/ALT_COUNT`，native 保留聚合 count evidence 并执行
`AlleleFractionModeller` 形状的 deterministic binomial responsibility conditional：
每段 minor fraction 与 mean-bias/bias-variance/outlier 全局参数均经 Kokkos
`RangePolicy` 更新，再用于 MAF posterior 初始化/步长；`.af.param` 的 GATK-shaped
字段来自该条件状态而非固定常数。CTest #126
`fastgatk-model-segments-allele-fraction-conditionals-gatk-oracle` 与 #125 copy-ratio
oracle 在 OpenMP/Serial 各自通过；Java Gamma-bias 边缘化、slice-sampler 随机抽样、
完整 latent state 与 credible-interval smoothing 仍是明确 fallback。

### 2026-09-01 GenotypeGVCFs multi-sample reference-confidence depth

发现并修正多样本 dense reference block 的真实数据流缺口：合并后的
FORMAT/DP 已是 sample-major 全 cohort 数组，但 INFO/DP 仍可能沿用第一个
source block 的值。native 现在在 GATK-compatible writer 边界按合并后的
FORMAT/DP 求 cohort 总深度；pinned GATK 4.6.2.0 两样本 oracle 覆盖 101
行，aggregate 与 `--stream-by-locus` 均整行 exact，OpenMP/Serial 各通过。
新增 CTest `fastgatk-genotype-gvcf-multisample-reference-confidence-gatk-oracle`。
该切片闭合的是多样本 reference-confidence annotation 语义，不代表完整
joint-model calibration 或 raw bit identity。

CombineGVCFs PL-less reference-confidence slice (2026-09-01)：native now
tracks whether any merged call has FORMAT/PL and does not synthesize missing
likelihood vectors for an all-PL-less cohort. For pure REF/<NON_REF> blocks,
the writer removes promoted INFO/DP while retaining END, matching pinned GATK
4.6.2.0. `fastgatk-combine-gvcfs-plless-gatk-oracle` compares the two-sample
compressed VCF rows exactly on OpenMP and Serial；PL-bearing paths remain
covered by the existing CombineGVCFs oracle。

### 2026-09-01 当前版本完整双 backend 回归

新增多样本 reference-confidence writer 后，初次整套回归发现既有具体变异记录的
INFO/DP 被误套 cohort FORMAT/DP 求和规则（GenotypeGVCFs #60/#61/#64）。native
现仅对 dense reference-only 记录（`allele_count==1`）执行该修正，具体变异仍保留
GATK variant annotation engine 的 DP 语义；受影响的 8 项 GenotypeGVCFs 测试在
OpenMP/Serial 重跑均 8/8。随后当前版本按 #1–#59、#60–#67、#68–#140 分段执行，
两套后端均合计 140/140 通过（OpenMP 2213.42 s、Serial 2266.39 s），未留下失败
测试进程。该回归闭合的是当前 CPU host/backend 门禁，不提升 registry prototype
为完整 Java 算法或 GPU/多节点生产等价。

### 2026-09-01 CreateReadCountPanelOfNormals degenerate-SVD boundary

GATK 4.6.2.0 rejects a multi-sample PoN when eigensamples were requested but
the standardized panel has no singular value above `EPSILON`. Native now
performs the same pre-write guard, so duplicated/identical normals cannot
produce an apparently valid HDF5 PoN that silently disables denoising. The
pinned `fastgatk-create-read-count-panel-of-normals-degenerate-svd-gatk-oracle`
checks Java/native rejection and absent partial output on OpenMP and Serial;
rank-zero and single-sample cases retain GATK's valid no-SVD behavior.

### 2026-09-01 新增边界回归收口

在此前当前版本 140/140 CPU 分段回归之外，本轮新增的
`fastgatk-combine-gvcfs-plless-gatk-oracle`、
`fastgatk-filter-mutect-calls-orientation-joint-oracle` 与
`fastgatk-create-read-count-panel-of-normals-degenerate-svd-gatk-oracle`，以及
GenotypeGVCFs maximum-ALT、ModelSegments smoothing、Mutect2 normal-evidence、HC AssemblyRegion
boundary、HC RCM PL-range、Mutect2 global-mismapping-rate、FilterMutectCalls normal-artifact 与
DenoiseReadCounts integer-input、PoN v7 sample-metadata 共十二项 pinned oracle，OpenMP/Serial 均各自通过；其中 PL-less oracle 还在首次 CTest 中暴露并修复了
`--stream-merge` 对纯 reference block 多写 `INFO/DP` 的 writer 回归。故当前
配置清单为 152 项，但“完整连续 CPU 基线”仍明确为 140/140，新增十二项以独立
oracle 证据计入，不虚构 152/152 的整套连续运行结果。

### 2026-09-01 Mutect2 global-mismapping-rate boundary

Pinned GATK 4.6.2.0 rejects `--phred-scaled-global-read-mismapping-rate 0` inside
`AlleleLikelihoods.normalizeLikelihoods`; negative values disable the cap and
positive values enable it. Native previously let zero reach the Kokkos
normalization helper as a zero-width cap, which would flatten per-read
likelihoods and change TLOD. `calling_pipeline.cpp` now rejects zero before
assembly while retaining the negative-disabled path. The new
`fastgatk-mutect2-mismapping-rate-boundary` oracle runs pinned Java and native
with the same fixture and passes on both OpenMP and Serial (CTest #93); it is a
bounded direct-replacement guard, not a claim that the remaining release-specific
Mutect2 posterior is bit-identical. The broad chr17 oracle remains 6/6 sites,
GT/DP/AD exact, with shared-site TLOD max absolute delta 3.500717.

The same boundary now preserves the disabled-cap representation end to end:
negative Q is passed as GATK's `Double.NEGATIVE_INFINITY` sentinel rather than a
large finite approximation. The shared Kokkos `normalize_likelihoods_kokkos`
kernel skips the cap in this mode while still producing best-read telemetry, so a
valid PairHMM `-Infinity` cell remains `-Infinity`. The bounded oracle now runs
both GATK 4.6.2.0 and native with Q=-1, and directly invokes the pinned Java
`AlleleLikelihoods` API on `[0.0,-Infinity,-3.0]`; OpenMP/Serial CTest #93 and the
Kokkos API smoke test both pass. This closes the disabled-cap representation
boundary but does not reduce the broad fixture's remaining TLOD delta, which is
still attributed to upstream assembly/PairHMM evidence inputs.

### 2026-09-01 Mutect2 tumor/normal evidence-width 修复

进一步审计 pinned GATK 4.6.2.0 `SomaticGenotypingEngine` 后确认，tumor 与
normal 的 `AlleleLikelihoods` 是独立 evidence 集合，低覆盖 normal 不应因其
fragment 数不同于 tumor 而被当作全 missing。native `calculate_posterior` 现在
保留 normal 自身 candidate-major 行宽，Kokkos posterior 对 normal evidence
单独遍历并按 unit 只计一次。新增
`fastgatk-somatic-posterior-normal-count-gatk-oracle`，固定 4 个 tumor/2 个
normal evidence units；OpenMP/Serial 均通过（CTest #150，定向运行）。这只收口
tumor/normal evidence shape，不宣称 release-specific posterior calibration 或
完整 FilterMutectCalls bit-identical。

### 2026-09-01 ModelSegments credible-interval smoothing

Native ModelSegments now implements the GATK adjacent-segment smoothing
decision when smoothing controls are explicitly requested. A Kokkos
`RangePolicy` evaluates same-contig candidates using posterior median
differences and either segment's 10--90% credible width; merged summaries use
GATK's inverse-variance normal approximation, followed by a deterministic
final refit on the merged partition. The pinned
`fastgatk-model-segments-smoothing-gatk-oracle` verifies merge and strict
non-merge boundaries against GATK 4.6.2.0 on both CPU backends. Intermediate
refits between multiple smoothing rounds remain an explicit bounded follow-up.

同日端到端资源基线（`benchmark_end_to_end.py`）在 chr17 fixture 上测得 native
HC 0.0812 s / 15,156 KiB RSS，Java GATK HC 4.058 s / 416,100 KiB RSS；由于
两边仍有 prototype-specific 参数差异，该结果只作为资源利用率观察，不作为最终
性能结论或生产 speedup 承诺。

### 2026-09-01 GenotypeGVCFs maximum-ALT joint-genotyping boundary

新增 `--max-alternate-alleles`（GATK 默认 6）的真实联合分型切片：Kokkos
likelihood-score kernel 按每个 sample 的最低 PL genotype 及其相对 hom-ref
距离排序 proper ALT，稳定保留原始顺序，并在 AF 计算前投影 PL/AD/GT 与
Number=A 字段。被删除 ALT 的原 genotype 变为 no-call，GQ 缺失，投影 PL
按 GATK 规则归一化到最小值 0。Pinned GATK 4.6.2.0 七 ALT→六 ALT oracle
在 aggregate 与 `--stream-by-locus` 两条路径、OpenMP/Serial 均 Java/native
数据行 exact；CTest 已注册为
`fastgatk-genotype-gvcf-max-alternate-alleles-gatk-oracle`。完整
release-specific cohort calibration/raw bit identity 仍为独立边界。

### 2026-09-01 GenotypeGVCFs discovered-ALT annotation boundary

新增可验证的 `--annotate-with-num-discovered-alleles` Boolean：native 在
`--max-alternate-alleles` 之前通过联合分型 writer 记录 concrete ALT 的
发现数，排除 `<NON_REF>`，纯 REF/reference-confidence block 不写 NDA。
该语义在现有 pinned GATK 4.6.2.0 七 ALT→六 ALT oracle 中固定为
`INFO/NDA=7`，并检查 GATK-compatible INFO 顺序、manifest option/telemetry；
aggregate 与 `--stream-by-locus` 两条路径的 OpenMP/Serial 均通过。该增量
复用 `fastgatk-genotype-gvcf-max-alternate-alleles-gatk-oracle`，不重复注册
CTest；完整 joint posterior calibration 仍为显式边界。

### 2026-09-01 FilterMutectCalls normal-artifact threshold boundary

Pinned GATK 4.6.2.0 对照发现，`NormalArtifactFilter` 的 NALOD posterior
是 ARTIFACT error probability，必须与 FilterMutectCalls 的 effective error
threshold 比较，而不是“只要大于 0 就过滤”。native 已将该判断接入
Kokkos 后验路径；`NALOD=-3`、normal AF=0.02、tumour AF=0.10、
`CONSTANT=0.1` 的 positive-but-below-threshold fixture（预期 posterior
`0.0089197225`）在 Java/native 均为 `PASS`。新增
`fastgatk-filter-mutect-calls-normal-artifact-oracle`，OpenMP/Serial 定向
CTest 均通过；完整 release-specific empirical clustering 与跨 filter
posterior calibration 仍为 fallback。

### 2026-09-01 DenoiseReadCounts integer COUNT boundary

审计 GATK 4.6.2.0 `SimpleCountCollection` 后确认 TSV `COUNT` 由
`DataLine.getInt()` 解码，而不是 double。native DenoiseReadCounts 现在在
进入 fractional-coverage/median/log2 Kokkos kernel 前拒绝 decimal、科学计数法、
负数及 signed-32-bit 溢出；合法整数 fixture 的 standardized values 与 Java
一致。新增 pinned oracle
`fastgatk-denoise-read-counts-integer-input-gatk-oracle`，OpenMP/Serial 均通过，
并已加入 CMake/`verify_all.sh`/dispatcher registry。`--normalization-target-coverage`
并非该版本 DenoiseReadCounts 参数，继续保持显式 fallback。

### 2026-09-01 逐工具 progress-score 审计

现有 `compute_progress_score.py` 已在不改变 global score
(`0.5385`/54%) 的前提下输出 `tool_audits`，覆盖 registry 的全部 48 项。
每项固定输出 registry `status`、五项机械门禁分数
(`api_cli`/`oracle`/`format_sidecar`/`resource_io`/`e2e_perf`)、证据脚本路径、
绑定 pinned fixture manifest 与证据脚本 bytes 的 tool-specific SHA-256
digest，以及 registry notes 提取的 fallback reason。证据映射采用显式表，
新增或删除 registry 条目、证据文件缺失、registry 字段不完整或 digest
manifest 缺失都会 fail-closed；逐工具分数只是发布审计视图，不反向修改
原有专家加权 workflow 分数。`fastgatk-progress-score-contract` 会执行这些
48 项完整校验。

### 2026-09-01 CreateReadCountPanelOfNormals PoN sample-path alignment

审计 GATK 4.6.2.0 的 HDF5-SVD PoN v7 后确认，原始与过滤后的
`sample_filenames` 都是输入 count 文件的绝对路径，并且分别与原始/保留
sample row 对齐，不是 `@RG SM` 标识。native preprocess 现在直接接收 CLI
输入的绝对路径；因此过滤 sample 后 `/original_data/sample_filenames` 和
`/panel/sample_filenames` 仍保持正确 row provenance。新增
`verify_create_read_count_panel_of_normals_sample_metadata_gatk_oracle.py`，
与 pinned GATK 4.6.2.0 writer 对照并在 OpenMP/Serial 均通过；完整
Spark/cloud 及 release-specific HDF5 byte parity 仍显式 fallback。

### 2026-09-01 ModelSegments allele-fraction initialization

ModelSegments 的 count-backed allele-fraction conditional 路径已由简单
`min(REF,ALT)/total` 初值推进到 GATK `AlleleFractionInitializer` 语义：
Kokkos `RangePolicy` 计算 `I(0.5, ALT+1, REF+1)` 的 alt-minor responsibility，
按 segment 聚合 responsibility-weighted minor reads，并应用 `+1/+2` flat-prior
pseudocount。manifest 记录 `allele_fraction_initial_segment_means`；独立
`verify_model_segments_allele_fraction_initialization_gatk_oracle.py` 用 pinned
GATK 4.6.2.0 锁定 `(2,1),(1,2)` 的 `29/64`，同时检查 Java/native partition/schema，
OpenMP/Serial 均通过。完整 slice sampler、latent Gamma bias 和 release-specific
MCMC posterior 仍为显式 fallback。

### 2026-09-01 GenotypeGVCFs GP/PG input compatibility

GenotypeGVCFs 现在在 cohort/annotation 消费完成后按 GATK writer 边界清理输入
`FORMAT/GP`、`FORMAT/PG` 的 stale vectors，保留 diagnostic profile 的显式输出；
`PREFER_PLS` 与 `<NON_REF>` allele-subset 路径不再把输入 GP/PG 泄漏到兼容输出。
新增 `verify_genotype_gvcf_gp_input_gatk_oracle.py`，以双样本 fixture 检查
`--gp-qual` gate、aggregate 与 `--stream-by-locus` 行级输出；OpenMP/Serial
GenotypeGVCFs 回归各 11/11 通过，CTest 清单新增
`fastgatk-genotype-gvcf-gp-input-gatk-oracle`。完整 release-specific joint
posterior calibration 仍为显式 fallback。

### 2026-09-01 GenotypeGVCFs STARTS_IN writer boundary

审计 pinned GATK 4.6.2.0 后确认
`--only-output-calls-starting-in-intervals` 只改变最终 VCF writer：输入 traversal 仍需
保留与 `-L` span-overlap 的记录。native 新实现因此不会提前丢弃区间外起点的长 deletion，
而在 Kokkos cohort/joint 计算和 Host annotation 完成后，才按 finalized POS 与标准化
include intervals 做 STARTS_IN 过滤。自包含 fixture 中 Java 默认输出 POS `[10,12]`，
启用该开关后输出 `[12]`；native aggregate 与 `--stream-by-locus` 的完整数据行均 exact。
新 CTest 在 OpenMP/Serial 分别 1/1 通过，manifest 的 mode、skipped 和实际 output count
均被断言。两构建 `ctest -N` 当前各列出 **204** 项；该证据不提升 release-specific
cohort calibration、raw-bit provenance 或 GPU production 语义的 fallback 范围。

本轮新增 ApplyVQSR 默认 cutoff 边界：pinned GATK 4.6.2.0 证明仅提供
`--tranches-file` 而没有 `--truth-sensitivity-filter-level` 时，不应把第一条
tranche 当作请求级别；Java 回退到有效 `VQSLOD >= 0.0`，低分记录使用
`LOW_VQSLOD`。native 已修正，新增
`fastgatk-apply-vqsr-default-cutoff-gatk-oracle`，OpenMP/Serial 均通过。显式
truth-sensitivity 仍进入 tranche walking；VQSR VBEM 全模型收敛仍是 fallback。

本轮 VariantRecalibrator VBEM 主模型增量：native 已对齐 GATK
`GaussianMixtureModel` 的 Java-Random 初始 mu/Sigma、K-means 空簇重启、
`I/200` empirical covariance、zero empirical mean、Wishart shrinkage（mean
outer product 只乘一次 `shrinkageFactor`）、post-E-step pMix 收敛以及
standard-deviation threshold 排除。扩展 pinned
`verify_variant_recalibrator_vbem_gatk_oracle.py` 为双高斯门禁：Java 4.6.2.0
正负模型 PMix/mean/covariance 最大绝对差 `5.23e-9`，六条 VQSLOD 最大差
`5e-5`；修复前 squared-factor 路径的最大 VQSLOD 差为 `1.22744`。OpenMP/Serial
的 CTest `fastgatk-variant-recalibrator-vbem-gatk-oracle` 均通过。Java
raw-bit VBEM、完整 resource calibration 与 recalibration-table provenance
仍明确为后续 fallback，不将该 bounded slice 宣称为全模型完成。

本轮补齐 CallCopyRatioSegments 输入分区校验：native 按 sequence dictionary
顺序并拒绝同 contig 重叠，pinned GATK 4.6.2.0 interval-validation oracle
在 OpenMP/Serial 均通过。该边界不代表完整 CNV/ModelSegments 等价。

### 2026-09-01 Mutect2 PairHMM/assembly residual audit

针对 pinned chr17 broad fixture 的 `17:69368` TLOD 最大绝对差
`3.500717`，先用 GATK 4.6.2.0 的 `--pair-hmm-results-file` 导出真实
read/haplotype/quality 请求，再由新增的
`fastgatk-pairhmm-results-gatk-oracle` 通过 `compute_kokkos_bucketed()` 重放。
GATK debug writer 的 706 行在 OpenMP 与 Serial 均通过，Kokkos 与其六位小数
科学计数法输出的最大差为 `4.222168e-6`（固定序列化边界 `5e-6`），因此当前
TLOD residual 不是 PairHMM recurrence 或 FP32/FP64 参数造成的：数值层已经
在同一请求输入上收敛。

同一 broad run 的 native stats 现在显式记录 `pairhmm_haplotypes`、graph
haplotype considered/pruned 和 `pairhmm_graph_snp_posterior_pairs`。固定 fixture
上 GATK 实际 PairHMM 请求为 706 行，native 为 724 对、同为 18 个 local
haplotype，request-set 差为 18；native graph haplotype 为 0。Mutect2 broad
oracle 因此严格保留该 assembly/request-set 差异与 TLOD 差异，未放宽 VCF 断言。
下一步应对齐 ReadThreadingAssembler/AssemblyRegion trimming 的逐 region
haplotype ownership，再用该 PairHMM replay oracle 复验；在此之前强行调 TLOD
或 posterior 阈值是不安全的。

### CombineGVCFs interval/reference-block boundary (2026-09-01)

Pinned GATK 4.6.2.0 probe shows that `-L` selection is not simple span
overlap for reference-confidence blocks: if an interval ends inside a
`REF/<NON_REF>` block, Java warns and drops the block; if the block `END` is
contained by the interval, Java emits the original full block. Native now
shares this predicate between aggregate and `--stream-merge`. The new
`fastgatk-combine-gvcfs-interval-refblock-gatk-oracle` compares both cases and
is required on OpenMP and Serial.

### ApplyVQSR sites-only writer boundary (2026-09-01)

GATK 4.6.2.0's `--sites-only-vcf-output` is a writer-only optional Boolean:
the full input FORMAT/sample payload remains available while VQSLOD and FILTER
decisions are made, then the published VCF contains only the eight site
columns. Native ApplyVQSR now duplicates the writer header and strips FORMAT at
the final write boundary. The pinned
`fastgatk-apply-vqsr-sites-only-gatk-oracle` compares Java/native rows and
manifest telemetry on OpenMP and Serial.

### 2026-09-01 FilterMutectCalls contamination-table sample fallback

审计 pinned GATK 4.6.2.0 `ContaminationFilter` 后确认，
`--contamination-table` 是按 tumor sample 名称查找的 per-sample map；表中不匹配
当前 tumor 的行不能被折叠为全局最大 contamination，未匹配样本应使用
`--contamination-estimate` 的 fallback（默认 0.0）。native 原先把所有表行的最大值
作为全局默认，可能把 OTHER 的 0.90 错施给 TUMOR；现已修正为保留 sample map、
contamination table 默认 0.0，并在显式 contamination estimate 存在时仍读取表，
让匹配行覆盖 fallback。pinned `verify_filter_mutect_calls_contamination_oracle.py`
新增 OTHER/TUMOR 不匹配 fixture，验证 Java/native contamination FILTER pattern，
OpenMP/Serial 均通过；完整 release-specific FilterMutectCalls joint calibration
仍保持显式 fallback。

### 2026-09-01 Mutect2 `--force-active` activity boundary

Pinned GATK 4.6.2.0 on the chr17 fixture keeps all 22 IGV profile rows and
their first four columns unchanged for `--force-active=false|true`, while each
size row changes from inactive `-1.00000` to active `1.00000`. Native previously
rejected the option; the shared Kokkos ActivityProfile API now preserves the
same boundaries and materializes every segment when enabled. The dedicated
`fastgatk-mutect2-force-active-gatk-oracle` checks Java/native profile state,
GATK's unchanged fixture VCF, native candidate-site-set stability, and
stats/manifest telemetry; OpenMP and Serial both pass. This is an activity/
AssemblyRegion ownership boundary, not a claim of complete somatic posterior
bit identity.
### 2026-09-01 Mutect2 `--normal-lod`

Audited GATK 4.6.2.0 `SomaticGenotypingEngine`: matched-normal ALT normal
log10 odds gates candidate emission, with default threshold 2.2. Native now
parses the value, applies it after Kokkos likelihood materialization, and
records it in stats/manifest. The pinned synthetic oracle (tumor 10/20 ALT,
normal 3/20 ALT at 17:69005) is exact at the stable site boundary: `-100`
emits `17:69005 T>C`, while `0` suppresses it, for Java and native on both
OpenMP and Serial. Floating posterior text is not claimed bit-identical.

### 2026-09-01 ModelSegments default KernelSegmenter boundary

审计发现 GATK 4.6.2.0 在未提供 `--segments` 时始终使用
`KernelSegmenter`，而 native 旧路径默认误用阈值分段。现已将无
`--change-point-threshold` 的单样本路径切换为 Kokkos KernelSegmenter；该
native 扩展显式给出 `--change-point-threshold` 时仍保留阈值模式，避免旧的
TSV 契约失效。线性 copy-ratio-only kernel 同时改用线性外积核 SVD 的精确
rank-one 特征 `x -> x`，不再把每个观测乘以随机 anchor 而改变 cost/penalty。
Pinned `fastgatk-model-segments-default-kernel-gatk-oracle` 用 step、gradient、
blip 多段 fixture 与 GATK 4.6.2.0 比对边界，OpenMP/Serial 均通过；Gaussian/
多维 SVD、完整 Java MCMC posterior/raw-bit 仍是明确 fallback。

> 2026-09-01 GenotypeGVCFs spanning-deletion ownership：pinned GATK 4.6.2.0
> 多样本 fixture 证明，若 `*` 没有 concrete deletion span 支持，Java 在
> joint AF/QUAL 后移除该 ALT，并将使用它的 sample 设为 `./.`、GQ=0、归一化
> PL；有真实 deletion span 时则保留 `*`。native 现在跨 shard 收集 deletion
> spans，在 aggregate 与 `--stream-by-locus` 最终 output-allele subset 复现该
> 有界规则；新增 `fastgatk-genotype-gvcf-spanning-deletion-gatk-oracle`，两种
> case 在 OpenMP/Serial 均 exact，manifest 记录 orphan 计数。完整 release-specific
> joint posterior calibration 仍保持 explicit fallback。

### DepthOfCoverage default read-filter boundary (2026-09-01)

Pinned GATK 4.6.2.0 `DepthOfCoverage.getDefaultReadFilters()` resolves
`WellformedReadFilter`, `NotDuplicateReadFilter`,
`NotSecondaryAlignmentReadFilter` and `MappedReadFilter`; QC-fail and
supplementary flags are not default exclusions. Native previously dropped
QC-fail records and silently counted the M islands of a CIGAR `N` read. The
native count path now retains QC-fail/supplementary reads and rejects a whole
record containing reference-skip, matching GATK's
`CigarContainsNoNOperator`. The new
`fastgatk-depth-of-coverage-read-filter-gatk-oracle` compares the locus table
byte-for-byte on OpenMP and Serial and checks `reads_seen=2`,
`reads_used=1`, `reads_filtered=1`.

### DepthOfCoverage deletion-site boundary (2026-09-01)

Pinned GATK 4.6.2.0 confirms the `--include-deletions` plus
`--ignore-deletion-sites` interaction: ignore suppresses deletion-only depth
and base-count contributions. Native now parses this optional boolean,
applies the effective setting before Kokkos projection, and records requested
versus effective deletion counting in the manifest. The new
`fastgatk-depth-of-coverage-ignore-deletion-sites-gatk-oracle` is exact on
both OpenMP and Serial; gene/cloud partitions remain fallback.

### Nextflow/SLURM direct-replacement boundary (2026-09-01)

审计发现通用 `nextflow_gatk_compat.nf` 的 output declaration 固定为
`result*`，因此自定义 `--output custom-counts.txt` 虽然 native 成功，仍会被
Nextflow 判定为缺失输出。现已改为按 caller output stem 声明并由本地真实
Nextflow 26.04.6 verifier 回归。`nextflow.config` 新增 container/cpus/memory
与 `slurm_gpu` profile；`slurm_smoke.sh` 增加 GPU count/type 到 `--gpus` 的
受控映射及显式 CUDA_VISIBLE_DEVICES 传递，fake-sbatch contract 通过。
这只证明本地 deterministic scheduler boundary；真实容器、CUDA/Kokkos
device execution、跨节点 SLURM/preemption 仍未宣称通过。

### ModelSegments combined allele-fraction segmentation audit (2026-09-01)

Java source 与 pinned GATK 4.6.2.0 共同证明 combined mode 使用 interval 内
第一个 allelic site 的 oriented ALT fraction。修复前 native 对全部 site 求平均
并折叠为 MAF，在 0.1→0.9 orientation 翻转但 MAF 恒 0.1 的 fixture 上错误只
输出一段；修复后 OpenMP/Serial 均与 Java 的两段边界一致。Manifest 固定
`first-oriented-alt-per-copy-ratio-interval` 语义；per-locus AF posterior count/
MCMC 数值仍未声明 exact。

### ApplyVQSR filter-control optional Boolean（2026-09-01）

pinned GATK 4.6.2.0 fixture 证明 `--ignore-all-filters` 与
`--exclude-filtered` 都是 optional Boolean。native 先前把它们当纯 flag，且
`ignore-all-filters` 未进入 pre-filtered 重评分 gate；现已修正并以 bare、true、
false 三种调用比较 FILTER/VQSLOD 行。该 scoped oracle 不改变 ApplyVQSR
prototype/fallback 总体状态。
