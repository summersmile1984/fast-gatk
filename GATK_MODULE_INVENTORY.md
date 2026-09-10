# GATK 4.6.2.0 源码模块清单与替换状态

各模块的 C++ Host、Kokkos kernel、内存/IO、验收门槛和实施依赖见
[MODULE_IMPLEMENTATION_PLANS.md](</home/turing-agents/Documents/fast-gatk/MODULE_IMPLEMENTATION_PLANS.md>)。

本清单由工作区中的 `src/main/java` 自动统计后人工复核，目标是避免只覆盖 HaplotypeCaller 主路径。状态含义：

- **native-priority**：应优先实现 Kokkos/HTSlib native 路径。
- **adapter**：保留成熟外部后端（GenomicsDB、Picard、HDF5 等），增加资源和兼容适配。
- **hybrid**：控制流/复杂对象保留 CPU，数值批处理使用 Kokkos。
- **fallback-first**：先由 Java GATK/Spark/Picard 执行，待 profile 后再移植。

## 顶层源码区域

| 区域 | Java 文件 | 约行数 | 命令行类 | 建议状态 | Kokkos 方案 |
|---|---:|---:|---:|---|---|
| `cmdline` | 44 | 3,451 | 0 | native-priority | C++/C ABI dispatcher、参数解析、错误码、版本/帮助契约；不要把 JVM 选项传给 device kernel |
| `engine` | 131 | 18,375 | 0 | hybrid | interval/shard/halo、reads/reference/features traversal；CPU I/O + Kokkos batch boundary |
| `exceptions` | 3 | 597 | 0 | native-priority | 统一错误分类：bad input、resource exhausted、backend unavailable、output contract failure |
| `metrics` | 19 | 1,783 | 0 | native-priority | 线程安全计数器、阶段时间、RSS/device/scratch/I/O telemetry |
| `tools` | 972 | 187,758 | 210 | 按子包分层 | 见下表；未完成命令必须 fallback |
| `transformers` | 15 | 1,323 | 0 | hybrid | read/variant transform 先 CPU；可批量且无状态的过滤再用 Kokkos |
| `utils` | 401 | 78,005 | 0 | 按子包分层 | 见下表 |

## `utils` 子模块

| 子包 | 文件/行数 | 主要职责 | 实现状态 |
|---|---:|---|---|
| `pairhmm` | 13 / 2,672 | regular/flow PairHMM、quality model | native-priority；ragged SoA、length bucket、team/vector/scratch |
| `smithwaterman` | 6 / 709 | SW 对齐和 CIGAR | native-priority；batched wavefront DP，保留 overhang/tie-break |
| `alignment`, `baq`, `clipping` | 5 / 1,681 | read alignment/BAQ/clipping | hybrid；短 read CPU SIMD，批量 BAQ 可选 device |
| `activityprofile`, `locusiterator`, `downsampling` | 24 / 3,783 | active region、pileup、downsample | hybrid；按 contig streaming，热点计数/过滤 Kokkos |
| `read`, `codecs`, `variant`, `fasta`, `reference` | 116 / 24,879 | read object、SAM/VCF codec、reference | adapter/native-priority；HTSlib/BGZF/libdeflate，C++ flat records |
| `io`, `nio`, `gcs`, `runtime` | 28 / 5,335 | 本地/云 I/O、运行时 | native-priority；bounded queues、prefetch、spill、cgroup/scratch probe |
| `recalibration` | 18 / 4,016 | BQSR covariates/tables | native-priority；SoA + thread-local integer reductions |
| `genotyper`, `haplotype`, `pileup` | 19 / 5,265 | likelihood/genotype/haplotype utilities | hybrid；复杂集合 CPU，小矩阵/归约 Kokkos |
| `mcmc`, `svd`, `clustering` | 22 / 1,358 | CNV/VQSR 数值运算 | hybrid；Kokkos Kernels dense math，控制流 CPU |
| `collections`, `iterators`, `functional`, `param` | 26 / 3,214 | 通用数据结构/迭代器 | CPU；native 侧避免 Java 风格对象图，使用 arena/flat vector |
| `spark` | 3 / 484 | Spark 工具辅助 | fallback-first；Nextflow+SLURM 代替 Spark 任务编排 |
| `bwa`, `dragstr`, `illumina`, `fragments` | 13 / 1,945 | 专项 read/STR 工具 | fallback-first；按 profile 选择独立 kernel |
| `report`, `tsv`, `text`, `help`, `logging`, `config`, `python`, `R` | 33 / 6,282 | 报告、文本、配置、脚本桥接 | CPU/fallback；不是 GPU 优先项 |

## 工具命令包

| 包 | 命令类数 | 代表命令 | 状态与方案 |
|---|---:|---|---|
| `tools/walkers/haplotypecaller` | 2 annotated commands / 87 files | HaplotypeCaller、RampedHaplotypeCaller | hybrid；host graph + Kokkos k-mer/SW/PairHMM |
| `tools/walkers/mutect` | 8 annotated commands / 57 files | Mutect2、FilterMutectCalls、GetPileupSummaries、LearnReadOrientationModel | hybrid；复用 HC/Mutect Kokkos kernel，VCF 语义与 sidecar Host；FilterMutect scalar boundary 走统一 KernelPlan |
| `tools/walkers/bqsr` | 4 | BaseRecalibrator、ApplyBQSR、GatherBQSRReports、AnalyzeCovariates | native-priority；计数 reduction + streaming transform |
| `tools/walkers/root` | 9 annotated commands / 12 direct files | GenotypeGVCFs、CombineGVCFs、ReblockGVCF、SelectVariants 等 | hybrid；locus batch math，writer/annotation CPU |
| `tools/walkers/vqsr` | 8 annotated commands / 32 files | VariantRecalibrator、ApplyVQSR、GatherTranches | hybrid；VBEM/GMM 用 Kokkos Kernels |
| `tools/walkers/annotator` | 1 annotated command / 83 files | VariantAnnotator 和各 annotation plugin | CPU/cache-first；批量数值 annotation 可后移植 |
| `tools/walkers/variantutils` | 8 annotated commands / 10 files | LeftAlignAndTrimVariants、VariantsToTable 等 | CPU/HTSlib adapter；保持 VCF semantics |
| `tools/walkers/varianteval` | 1 annotated command / 56 files | VariantEval 及 evaluator/stratification | CPU；复杂插件和报告逻辑不适合第一阶段 GPU |
| `tools/walkers/sv` | 7 annotated commands / 9 files | SVAnnotate、SVCluster、SVConcordance、SVStratify | fallback-first；证据计数/聚类后续选择性加速 |
| `tools/walkers/coverage`、`fasta`、`filters`、`validation` 等 | 其余 walkers | DepthOfCoverage、CountReads、FlagStat、SplitIntervals、CountBasesInReference、FastaReferenceMaker、FastaAlternateReferenceMaker、ShiftFasta、IndexFeatureFile、Validate* 等 | CPU/HTSlib；固定字节/过滤 predicate 可批量化但收益低 |
| `tools/spark` | 40 | MarkDuplicatesSpark、SortSamSpark、BwaSpark、PathSeq* | fallback-first；native 不复刻 Spark API，改由 SLURM scatter 或单节点 pipeline |
| `tools/copynumber` | 19 | CollectReadCounts、DenoiseReadCounts、ModelSegments、GermlineCNVCaller | hybrid；counts/SVD/dense math Kokkos，HDF5/模型控制 CPU |
| `tools/funcotator` | 4 | Funcotator、FuncotateSegments、FilterFuncotations | CPU/cache-first；数据源查询和字符串处理不适合 GPU |
| `tools/genomicsdb` | 1 | GenomicsDBImport | adapter；复用 GenomicsDB/TileDB，自适应 batch/tmp/fd |
| `tools/sv` | 4 | CondenseDepthEvidence、PrintSVEvidence 等 | fallback-first；先保持输出契约 |
| `tools/reference` | 2 | CheckReferenceCompatibility、CompareReferences | CPU/I/O；reference cache 和并行 checksum |
| `tools/dragstr` | 2 | CalibrateDragstrModel、ComposeSTRTableFile | CPU/fallback；后续可批量化模型计算 |

## 外部依赖边界

GATK 的构建脚本声明了以下不在 `src/main/java` 中完整实现的能力：

| 依赖 | 兼容策略 |
|---|---|
| htsjdk | 用 HTSlib C API 做 native I/O；保留 SAM/BAM/CRAM/VCF/BCF/index 语义 |
| Picard | Picard 命令先以 Java fallback；MarkDuplicates/SortSam 另做 external-memory native engine |
| Intel GKL | 作为 CPU PairHMM/SW/DEFLATE baseline；native Kokkos 结果必须与其做 oracle 比较 |
| GenomicsDB/TileDB | 第一阶段保留 backend；只替换 adapter、分区、batch、tmp 和资源监控 |
| Spark/Hadoop | 不在 Kokkos runtime 内重建；Nextflow+SLURM 负责作业级散列与聚合 |
| Google Cloud NIO/GCS connector | 先由 HTSlib/云 URI adapter + prefetch 支持；缓存计入内存预算 |
| HDF5/Parquet/Snappy | CNV/报告路径保留格式；只有 dense math 进入 Kokkos |

## Native registry 建议

```text
tool_registry:
  HaplotypeCaller: native(hc_v1), fallback(gatk)
  Mutect2: native(mutect_v1), fallback(gatk)
  BaseRecalibrator: native(bqsr_v1), fallback(gatk)
  ApplyBQSR: native(apply_bqsr_v1), fallback(gatk)
  GenotypeGVCFs: adapter(genomicsdb + native_genotype_v1), fallback(gatk)
  GenomicsDBImport: adapter(genomicsdb), fallback(gatk)
  MarkDuplicates: fallback(picard), native(md_external_v1 when enabled)
  all other tools: fallback(gatk) until a tool-specific contract is registered
```

每个 registry entry 应包含：支持的参数集合、输出/索引/sidecar 列表、deterministic/fast 模式、可用 backend、最小测试集、资源估算器和失败回退策略。未知参数不能静默忽略。
