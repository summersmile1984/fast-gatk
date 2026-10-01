# fast-gatk 项目范围声明（Scope Statement）

> 生成日期：2026-09-17
> 来源：grep 实证 + 主计划 `FAST_GATK_EXECUTION_PLAN.md` §3.2 / §3.3 / §4.2 / §5 / §6 / §10
> 主计划引用："在保持 GATK 命令、文件和结果契约的前提下，用 **C++ Host + Kokkos Kernel + HTSlib** 实现可在 CPU/GPU/异构集群运行的 native toolkit，并能被现有 **GATK + SLURM + Nextflow** 直接调用。"

本文件是项目范围的**权威声明**。所有"我们是否要实现 X"的问题，以本文档为准。
如需修改范围，必须先改本文档并提交对应的代码/CI 变更。

---

## 1. 范围目标（IN-SCOPE）

| 维度 | 内容 | 代码位置 |
|---|---|---|
| Host 实现 | C++20 + HTSlib（无 Rust、无 C ABI 热路径） | `fastgatk-native/src/` + `fastgatk-core/` + `fastgatk-runtime/` |
| Kernel 实现 | Kokkos 5.2.0（`Kokkos::View` / `MDRangePolicy` / `Experimental::simd`） | `fastgatk-kernels/src/` |
| Kokkos backend | SERIAL / OPENMP / CUDA / HIP / SYCL（含 `*_OPENMP` 组合） | `fastgatk-native/CMakeLists.txt:18-22` |
| 调度编排 | Nextflow DSL2 + SLURM scatter/gather + 自带 `slurm_smoke.sh` 包装 | `fastgatk-native/workflow/` |
| 文件 IO | 本地 BAM/CRAM/VCF/BCF/GVCF/TSV/HDF5（通过 HTSlib） | `fastgatk-native/src/*_tool.cpp` |
| CLI 兼容 | `gatk <Tool> ...` 形式的 dispatcher，CLI 兼容、未知参数 fail-closed | `fastgatk-native/dispatcher/` |
| GATK 4.6.2.0 工具子集 | 48 个 contract-compatible + 1 个 adapter（详见 §3） | `fastgatk-native/dispatcher/tool_registry.json` |
| 应用域 | 基因组测序（短读 WGS/WES、RNA-seq、靶向 panel） | 由 HaplotypeCaller / Mutect2 / BQSR / GenomicsDB 链路决定 |

---

## 2. 明确**不在**范围（OUT-OF-SCOPE）

以下项目**不实现**，且不打算实现。任何"实现 X"的提议需先修改本文档。

### 2.1 GATK Spark 工具家族（68 个 `*Spark.java`）

`gatk-source/src/main/java/org/broadinstitute/hellbender/tools/**/*Spark.java` 共 68 个，全部不实现：

**Caller Spark（11 个）：**
HaplotypeCallerSpark、BQSRPipelineSpark、ReadsPipelineSpark、PrintReadsSpark、PrintVariantsSpark、PileupSpark、BwaAndMarkDuplicatesPipelineSpark、CountReadsSpark、CountBasesSpark、CountVariantsSpark、FlagStatSpark、SortSamSpark、RevertSamSpark、MarkDuplicatesSpark

**BQSR Spark（4 个）：**
BaseRecalibratorSpark、BaseRecalibratorSparkFn、ApplyBQSRSpark、ApplyBQSRSparkFn

**CNV/SV Spark（12 个）：**
CollectAllelicCountsSpark、StructuralVariationDiscoveryPipelineSpark、FindBreakpointEvidenceSpark、ExtractSVEvidenceSpark、CalcMetadataSpark、CpxVariantReInterpreterSpark、DiscoverVariantsFromContigAlignmentsSAMSpark、SvDiscoverFromLocalAssemblyContigAlignmentsSpark、ExtractOriginalAlignmentRecordsByNameSpark、FindBadGenomicKmersSpark、CompareDuplicatesSpark

**PathSeq Spark（7 个）：**
PathSeqPipelineSpark、PathSeqScoreSpark、PathSeqFilterSpark、PathSeqBwaSpark、PSBwaAlignerSpark、PSPairedUnpairedSplitterSpark、ContainsKmerReadFilterSpark

**BWA Spark（2 个）：**
BwaSpark、BwaSparkEngine

**Metrics Spark（8 个）：**
InsertSizeMetricsCollectorSpark、QualityYieldMetricsCollectorSpark、MeanQualityByCycleSpark、QualityScoreDistributionSpark、CollectBaseDistributionByCycleSpark、CollectInsertSizeMetricsSpark、CollectQualityYieldMetricsSpark、CollectMultipleMetricsSpark、MetricsCollectorSpark、MetricsCollectorSparkTool

**Hopscotch & Sparkifier（8 个）：**
HopscotchCollectionSpark、HopscotchMapSpark、HopscotchMultiMapSpark、HopscotchSetSpark、HopscotchUniqueMultiMapSpark、ReadFilterSparkifier、ReadTransformerSparkifier、DownsampleableSparkReadShard

**Example Spark（11 个，纯示例）：**
ExampleAssemblyRegionWalkerSpark、ExampleCollectMultiMetricsSpark、ExampleCollectSingleMetricsSpark、ExampleIntervalWalkerSpark、ExampleLocusWalkerSpark、ExampleMultiMetricsCollectorSpark、ExampleReadWalkerWithReferenceSpark、ExampleReadWalkerWithVariantsSpark、ExampleSingleMetricsCollectorSpark、ExampleVariantWalkerSpark

**工具型辅助（1 个）：**
ParallelCopyGCSDirectoryIntoHDFSSpark

### 2.2 Spark 运行时基础设施

| 项目 | 状态 | 证据 |
|---|---|---|
| `SparkContext` / `JavaSparkContext` | ❌ 无 | `grep -rE 'SparkContext\|JavaSparkContext' fastgatk-{native,kernels,runtime,core}/src/` → 0 命中 |
| RDD / Dataset<Row> | ❌ 无 | 同上 |
| Spark on YARN/K8s 提交器 | ❌ 无 | `third_party/` 不包含 spark 依赖 |
| Spark serializer / Kryo 注册 | ❌ 无 | 同上 |
| Hadoop FileSystem API | ❌ 无 | `grep -rE 'hadoop\|hdfs://' fastgatk-*/src/` → 0 命中 |
| Dataproc / EMR / Glue launcher | ❌ 无 | 无对应实现 |
| `third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-spark.jar` | ⚠️ vendored 不 link | GATK 4.6.2.0 release tarball 自带 Spark 变体 jar，仅作 oracle 比对参考存在；我们仅 link `gatk-package-4.6.2.0-local.jar`，从不 link 该 `-spark.jar`。`fastgatk-scope-contract` 显式豁免此文件 |

### 2.3 云对象存储

| 协议 | 状态 | 证据 |
|---|---|---|
| `s3://` | ❌ 拒绝 | `fastgatk-native/src/genotype_gvcf_tool.cpp:352-354`、`genomicsdb_import_tool.cpp:573-574` —— 命中后抛 `BACKEND_UNAVAILABLE` |
| `gs://` | ❌ 拒绝 | 同上 |
| `https://` / `http://` | ❌ 拒绝 | 同上 |
| `azure://` / `wasb://` | ❌ 拒绝 | 同上（未实现 + 拒绝） |
| 工具级 cloud fallback | 19 个工具显式列为 fallback | `tool_registry.json` 中 `cloud-*` 共 19 处（详见 §3） |

### 2.4 与基因组测序无关的 GATK 工具族

| 工具族 | 状态 | 主计划 §3.3 行 |
|---|---|---|
| PathSeq（病原体分类） | ❌ 不实现 | "P3 fallback" |
| BWA-MEM（read alignment） | ❌ 不实现（外部依赖） | "P3 fallback" |
| StructuralVariation 全链路 | ❌ 不实现 | "P3 fallback" |
| Funcotator / FuncotateSegments（癌症注释） | ❌ 不实现 | "P3 CPU/cache-first，无 native" |
| MuTect / M1（旧体细胞） | ❌ 不实现 | 只做 Mutect2 |
| DRAGSTR（STR 建模） | ❌ 不实现 | "P3 CPU/fallback" |
| GenomicsDB Java workspace（TileDB） | ⚠️ Adapter only | `GenomicsDBImport` 是唯一 `status: "adapter"`，产物仍走 fastgatk-inputs.tsv |

---

## 3. 工具分类（48 个 + 1 adapter）的"Spark/cloud 是否在范围"

每个工具的 `fallback_boundaries` 字段列出**实际不实现**的范围。下列汇总：

### 3.1 Spark 显式 fallback（4 个）
| 工具 | Spark fallback |
|---|---|
| BaseRecalibrator | `spark-distributed-bqsr` |
| ApplyBQSR | `spark-distributed` |
| CollectF1R2Counts | `spark-distributed-collection` |
| VariantRecalibrator | `spark-distributed-recalibration` |

### 3.2 Cloud 显式 fallback（19 个，按工具）
- `cloud-scatter-gather-orchestration`：GatherBQSRReports、GatherTranches
- `cloud-report-rendering`：AnalyzeCovariates
- `cloud-filtering-stats`：FilterMutectCalls
- `cloud-em-distribution`：LearnReadOrientationModel
- `cloud-indexed-feature-semantics`：CollectAllelicCounts
- `cloud-inputs`：DepthOfCoverage、VariantFiltration
- `cloud-scatter-segmentation`：ModelSegments
- `cloud-scatter-application`：ApplyVQSR
- `cloud-object-store-staging`：ReblockGVCF、SortSam、MarkDuplicates、CombineGVCFs
- `cloud-remote-feature-inputs`：VariantEval
- `cloud-uri-staging`：GatherPileupSummaries
- `cloud-report-bit-format`：CalculateContamination
- `cloud-and-dictionary-edge-cases`：GetPileupSummaries
- `cloud-remote-inputs`：ValidateVariants
- `interval-list-cloud-semantics`：LeftAlignAndTrimVariants

### 3.3 不含 Spark/cloud 显式 fallback（仍 out-of-scope 视情况而定）
- 16 个工具只有 `compatibility_scope` 未补齐（HaplotypeCaller、AnnotateIntervals、CountBasesInReference、CompareReferences、CheckReferenceCompatibility、FastaReferenceMaker、FastaAlternateReferenceMaker、ShiftFasta、IndexFeatureFile、CountReads、FlagStat、SplitIntervals、FilterIntervals、PreprocessIntervals、GatherVcfs、GenomicsDBImport）—— 这些是**审计未完成**，不是范围问题。

---

## 4. 范围内 - 但**仅在调度层**

### 4.1 Nextflow + SLURM 编排（**纯 shell wrapper**，不是 native 实现）
`fastgatk-native/workflow/` 下 7 个文件：
- `scatter_gather_smoke.sh`（interval 切分 + 并发跑 native binary + gather）
- `slurm_smoke.sh`（检测 `SLURM_CPUS_PER_TASK` 等环境，检测到已分配作业就不嵌套 `sbatch`）
- `nextflow.config`、`nextflow_scatter_gather.nf`、`nextflow_cluster_smoke.nf`、`nextflow_gatk_compat.nf`、`nextflow_smoke.nf`

→ 全部是 **bash / Nextflow DSL2**，**不引入 Spark 运行时**。
→ WP-15（主计划 §4.2）的范围就是这一层。

### 4.2 GPU/CPU 异构（Kokkos 后端，不是 Spark）
CMake 已配置 SERIAL / OPENMP / CUDA / HIP / SYCL（含 `*_OPENMP` 组合）：
- 当前 `build/` 实际启用：`Kokkos_ENABLE_OPENMP=ON, Kokkos_ENABLE_SERIAL=ON`，其余 OFF
- GPU 后端走 Kokkos native（`pairhmm_kokkos.cpp`、`smith_waterman_kokkos.cpp`、`filter_model_kokkos.cpp`），不是 Spark on GPU

### 4.3 Dispatcher 对 Spark 参数的 fail-closed 行为
- `--gatk-config-file`（可携带 Spark/cloud 默认值）→ fail-closed：`fastgatk-native/dispatcher/fastgatk.py:963-975`
- `--`（GATK Spark tool/Spark argument 分割符）→ fail-closed：`fastgatk-native/dispatcher/fastgatk.py:1014-1023`

---

## 5. 范围变更流程

任何对 §2（OUT-OF-SCOPE）内容的添加或删除必须：
1. 修改本文件 `SCOPE.md` 对应章节
2. 在 `tool_registry.json` 对应工具的 `fallback_boundaries` 加/删条目
3. 同步 `FAST_GATK_EXECUTION_PLAN.md` §3.2 / §3.3
4. 通过 `progress_score.json` 的 `workflow_weights.gpu_slurm_cluster` 重新加权
5. 提交 PR 并由至少一位维护者 review

---

## 6. 验证命令

下列命令任何人都可以本地复现，验证本文件 §2 的每一项：

```bash
# 2.1 Spark 工具清单 (60 = 50 实际工具 + 10 Example)
find gatk-source/src/main/java -name "*Spark*.java" -path "*/tools/*" | wc -l  # 期望 60
# 非 example 的真正工具
find gatk-source/src/main/java -name "*Spark*.java" -path "*/tools/*" ! -name "Example*" | wc -l  # 期望 50
# Example 类
find gatk-source/src/main/java -name "Example*Spark*.java" -path "*/tools/*" | wc -l  # 期望 10

# 2.2 Spark 运行时
grep -rE 'SparkContext|JavaSparkContext|RDD|Dataset<' fastgatk-{native,kernels,runtime,core}/src/ | wc -l  # 期望 0
grep -rE 'hadoop|hdfs://' fastgatk-{native,kernels,runtime,core}/src/ | wc -l  # 期望 0
ls third_party/ | grep -iE 'spark|hadoop|hdfs' | wc -l  # 期望 0

# 2.3 云对象存储 - 应只在拒绝逻辑中出现
grep -rnE 'rfind\("s3://"|rfind\("gs://"|rfind\("https?://"' fastgatk-native/src/  # 期望 4 行（全在 genotype_gvcf_tool.cpp / genomicsdb_import_tool.cpp 的拒绝分支）

# 2.3 工具级 cloud fallback
grep -c '"cloud' fastgatk-native/dispatcher/tool_registry.json  # 期望 19
grep -c '"spark' fastgatk-native/dispatcher/tool_registry.json  # 期望 4

# 4.3 Dispatcher fail-closed
grep -nE 'UNSUPPORTED_PARAMETER|spark' fastgatk-native/dispatcher/fastgatk.py | head -10

# 完整契约门禁 - 上方 grep 的 CI 强制版本（CTest 名: fastgatk-scope-contract）
python3 fastgatk-native/scripts/verify_scope_contract.py
```

如果任何一条不符，本文件必须更新，且对应修正 `fastgatk-native/scripts/verify_scope_contract.py` 中固化期望值。
