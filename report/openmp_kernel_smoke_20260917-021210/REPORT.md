# OpenMP Kernel Smoke Verification

生成时间：20260917-021210（Asia/Taipei）
工作树：/home/turing-agents/Documents/fast-gatk
构建目录：/home/turing-agents/Documents/fast-gatk/fastgatk-native/build（OpenMP 后端，OMP_NUM_THREADS=4）
方法：对每个工具跑最小 invocation，解析 manifest 中的 ``telemetry.execution_space`` 或 ``*.kernel_*_execution_space`` 字段。

## 汇总

- ✅ OpenMP 验证（``telemetry.execution_space == OpenMP``）：1 个工具
- ✅ OpenMP kernel 命中（manifest 含 ``*_execution_space=OpenMP``）：4 个工具（含上）
- ⚪ rc=0 但无 manifest（HTSlib/Host 路径，无 Kokkos kernel）：31 个工具
- 🟡 Serial/Host-only（registry 声明）：0 个工具
- ⚠️ 跳过（fixture 缺失）：13 个工具
- ❌ 失败：0 个工具

## 详细结果

| Tool | rc | wall(s) | execution_space | kernel_spaces | status |
| --- | ---: | ---: | --- | --- | --- |
| HaplotypeCaller | 0 | 0.019 | OpenMP | execution_space=OpenMP, pairhmm_normalization_execution_space=OpenMP, pairhmm_marginalization_execut | ✅ OpenMP (manifest) |
| BaseRecalibrator | 0 | 0.612 | — | kernel_execution_space=OpenMP, covariate_kernel_execution_space=OpenMP | ✅ OpenMP kernels (2) |
| ApplyBQSR | 0 | 0.031 | — | — | ⚪ HTSlib/Host path |
| Mutect2 | 0 | 0.085 | — | somatic_likelihood_execution_space=OpenMP, somatic_posterior_execution_space=OpenMP, error_correctio | ✅ OpenMP kernels (2) |
| FilterMutectCalls | 0 | 0.013 | — | — | ⚪ HTSlib/Host path |
| CollectAllelicCounts | 0 | 0.006 | — | — | ⚪ HTSlib/Host path |
| CollectF1R2Counts | 0 | 0.012 | — | — | ⚪ HTSlib/Host path |
| CollectReadCounts | 0 | 0.009 | — | — | ⚪ HTSlib/Host path |
| CombineGVCFs | - | - | - | - | ⚠️ skipped: missing fixture dependency |
| GenotypeGVCFs | 0 | 0.009 | — | genotype_kernel_execution_space=OpenMP, pl_remap_kernel_execution_space=OpenMP, allele_count_kernel_ | ✅ OpenMP kernels (7) |
| ReblockGVCF | 0 | 0.008 | — | — | ⚪ HTSlib/Host path |
| SelectVariants | 0 | 0.009 | — | — | ⚪ HTSlib/Host path |
| VariantFiltration | 0 | 0.008 | — | — | ⚪ HTSlib/Host path |
| LeftAlignAndTrimVariants | 0 | 0.008 | — | — | ⚪ HTSlib/Host path |
| ValidateVariants | 0 | 0.003 | — | — | ⚪ HTSlib/Host path |
| VariantsToTable | 0 | 0.003 | — | — | ⚪ HTSlib/Host path |
| VariantEval | 0 | 0.003 | — | — | ⚪ HTSlib/Host path |
| GatherVcfs | 0 | 0.008 | — | — | ⚪ HTSlib/Host path |
| GatherTranches | 0 | 0.001 | — | — | ⚪ HTSlib/Host path |
| SortSam | 0 | 0.010 | — | — | ⚪ HTSlib/Host path |
| MarkDuplicates | 0 | 0.010 | — | — | ⚪ HTSlib/Host path |
| GatherBQSRReports | 0 | 0.005 | — | — | ⚪ HTSlib/Host path |
| AnalyzeCovariates | 0 | 0.009 | — | — | ⚪ HTSlib/Host path |
| CountReads | 0 | 0.003 | — | — | ⚪ HTSlib/Host path |
| FlagStat | 0 | 0.003 | — | — | ⚪ HTSlib/Host path |
| SplitIntervals | 0 | 0.003 | — | — | ⚪ HTSlib/Host path |
| FilterIntervals | - | - | - | - | ⚠️ skipped: missing fixture dependency |
| PreprocessIntervals | 0 | 0.002 | — | — | ⚪ HTSlib/Host path |
| AnnotateIntervals | 0 | 0.002 | — | — | ⚪ HTSlib/Host path |
| CountBasesInReference | 0 | 0.005 | — | — | ⚪ HTSlib/Host path |
| CompareReferences | 0 | 0.002 | — | — | ⚪ HTSlib/Host path |
| CheckReferenceCompatibility | 0 | 0.002 | — | — | ⚪ HTSlib/Host path |
| FastaReferenceMaker | 0 | 0.018 | — | — | ⚪ HTSlib/Host path |
| FastaAlternateReferenceMaker | 0 | 0.020 | — | — | ⚪ HTSlib/Host path |
| ShiftFasta | 0 | 0.018 | — | — | ⚪ HTSlib/Host path |
| IndexFeatureFile | 0 | 0.005 | — | — | ⚪ HTSlib/Host path |
| GetPileupSummaries | - | - | - | - | ⚠️ skipped: missing fixture dependency |
| CalculateContamination | - | - | - | - | ⚠️ skipped: missing fixture dependency |
| LearnReadOrientationModel | - | - | - | - | ⚠️ skipped: missing fixture dependency |
| DepthOfCoverage | 0 | 0.003 | — | — | ⚪ HTSlib/Host path |
| ApplyVQSR | - | - | - | - | ⚠️ skipped: missing fixture dependency |
| VariantRecalibrator | - | - | - | - | ⚠️ skipped: missing fixture dependency |
| GenomicsDBImport | - | - | - | - | ⚠️ skipped: missing fixture dependency |
| GatherPileupSummaries | - | - | - | - | ⚠️ skipped: missing fixture dependency |
| CallCopyRatioSegments | - | - | - | - | ⚠️ skipped: missing fixture dependency |
| ModelSegments | - | - | - | - | ⚠️ skipped: missing fixture dependency |
| DenoiseReadCounts | - | - | - | - | ⚠️ skipped: missing fixture dependency |
| CreateReadCountPanelOfNormals | - | - | - | - | ⚠️ skipped: missing fixture dependency |

详细 JSON：`/home/turing-agents/Documents/fast-gatk/report/openmp_kernel_smoke_20260917-021210/summary.json`