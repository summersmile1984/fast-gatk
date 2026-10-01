# OpenMP Kernel Smoke Verification

生成时间：20260917-020930（Asia/Taipei）
工作树：/home/turing-agents/Documents/fast-gatk
构建目录：/home/turing-agents/Documents/fast-gatk/fastgatk-native/build（OpenMP 后端，OMP_NUM_THREADS=4）
方法：对每个工具跑最小 invocation，解析 manifest 中的 ``telemetry.execution_space``。

## 汇总

- ✅ OpenMP kernel 验证：35 个工具
- 🟡 Serial/Host-only（registry 声明，非 OpenMP）：0 个工具
- ⚠️ 跳过（fixture 缺失）：12 个工具
- ❌ 失败：1 个工具

## 详细结果

| Tool | rc | wall(s) | execution_space | kernel_spaces | status |
| --- | ---: | ---: | --- | --- | --- |
| HaplotypeCaller | 0 | 0.027 | OpenMP | execution_space=OpenMP, pairhmm_normalization_execution_space=OpenMP, pairhmm_marginalization_execution_space=OpenMP, pairhmm_uncertainty_execution_space=OpenMP, pairhmm_execution_space=OpenMP, sw_execution_space=OpenMP, graph_haplotype_sw_execution_space=OpenMP, graph_execution_space=OpenMP, activity_execution_space=OpenMP, read_filter_execution_space=OpenMP, reference_confidence_execution_space=OpenMP, call_confidence_execution_space=OpenMP | ✅ OpenMP |
| BaseRecalibrator | 0 | 0.616 | — | kernel_execution_space=OpenMP, covariate_kernel_execution_space=OpenMP | ⚠️ no manifest |
| ApplyBQSR | 0 | 0.026 | — | — | ⚠️ no manifest |
| Mutect2 | 0 | 0.079 | — | somatic_likelihood_execution_space=OpenMP, somatic_posterior_execution_space=OpenMP | ⚠️ no manifest |
| FilterMutectCalls | 0 | 0.010 | — | — | ⚠️ no manifest |
| CollectAllelicCounts | 0 | 0.004 | — | — | ⚠️ no manifest |
| CollectF1R2Counts | 0 | 0.011 | — | — | ⚠️ no manifest |
| CollectReadCounts | 0 | 0.009 | — | — | ⚠️ no manifest |
| CombineGVCFs | - | - | - | - | ⚠️ missing fixture dependency |
| GenotypeGVCFs | 0 | 0.009 | — | genotype_kernel_execution_space=OpenMP, pl_remap_kernel_execution_space=OpenMP, allele_count_kernel_execution_space=OpenMP, cross_sample_reference_execution_space=OpenMP, allele_field_remap_kernel_execution_space=OpenMP, cohort_af_kernel_execution_space=OpenMP, genotype_prior_kernel_execution_space=OpenMP | ⚠️ no manifest |
| ReblockGVCF | 0 | 0.008 | — | — | ⚠️ no manifest |
| SelectVariants | 0 | 0.009 | — | — | ⚠️ no manifest |
| VariantFiltration | 0 | 0.009 | — | — | ⚠️ no manifest |
| LeftAlignAndTrimVariants | 0 | 0.009 | — | — | ⚠️ no manifest |
| ValidateVariants | 0 | 0.002 | — | — | ⚠️ no manifest |
| VariantsToTable | 0 | 0.003 | — | — | ⚠️ no manifest |
| VariantEval | 0 | 0.003 | — | — | ⚠️ no manifest |
| GatherVcfs | 0 | 0.008 | — | — | ⚠️ no manifest |
| GatherTranches | 0 | 0.003 | — | — | ⚠️ no manifest |
| SortSam | 0 | 0.010 | — | — | ⚠️ no manifest |
| MarkDuplicates | 0 | 0.011 | — | — | ⚠️ no manifest |
| GatherBQSRReports | 0 | 0.005 | — | — | ⚠️ no manifest |
| AnalyzeCovariates | 0 | 0.008 | — | — | ⚠️ no manifest |
| CountReads | 0 | 0.002 | — | — | ⚠️ no manifest |
| FlagStat | 0 | 0.002 | — | — | ⚠️ no manifest |
| SplitIntervals | 0 | 0.002 | — | — | ⚠️ no manifest |
| FilterIntervals | - | - | - | - | ⚠️ missing fixture dependency |
| PreprocessIntervals | 0 | 0.002 | — | — | ⚠️ no manifest |
| AnnotateIntervals | 0 | 0.002 | — | — | ⚠️ no manifest |
| CountBasesInReference | 0 | 0.005 | — | — | ⚠️ no manifest |
| CompareReferences | 0 | 0.002 | — | — | ⚠️ no manifest |
| CheckReferenceCompatibility | 0 | 0.001 | — | — | ⚠️ no manifest |
| FastaReferenceMaker | 0 | 0.016 | — | — | ⚠️ no manifest |
| FastaAlternateReferenceMaker | 0 | 0.017 | — | — | ⚠️ no manifest |
| ShiftFasta | 0 | 0.017 | — | — | ⚠️ no manifest |
| IndexFeatureFile | 0 | 0.005 | — | — | ⚠️ no manifest |
| GetPileupSummaries | - | - | - | - | ⚠️ missing fixture dependency |
| CalculateContamination | - | - | - | - | ⚠️ missing fixture dependency |
| LearnReadOrientationModel | - | - | - | - | ⚠️ missing fixture dependency |
| DepthOfCoverage | 0 | 0.003 | — | — | ⚠️ no manifest |
| ApplyVQSR | - | - | - | - | ⚠️ missing fixture dependency |
| VariantRecalibrator | 2 | 0.946 | — | — | ❌ rc=2 |
| GenomicsDBImport | - | - | - | - | ⚠️ missing fixture dependency |
| GatherPileupSummaries | - | - | - | - | ⚠️ missing fixture dependency |
| CallCopyRatioSegments | - | - | - | - | ⚠️ missing fixture dependency |
| ModelSegments | - | - | - | - | ⚠️ missing fixture dependency |
| DenoiseReadCounts | - | - | - | - | ⚠️ missing fixture dependency |
| CreateReadCountPanelOfNormals | - | - | - | - | ⚠️ missing fixture dependency |

详细 JSON：`/home/turing-agents/Documents/fast-gatk/report/openmp_kernel_smoke_20260917-020930/summary.json`