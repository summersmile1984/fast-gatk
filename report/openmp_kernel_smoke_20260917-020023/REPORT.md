# OpenMP Kernel Smoke Verification

生成时间：20260917-020023（Asia/Taipei）
工作树：/home/turing-agents/Documents/fast-gatk
构建目录：/home/turing-agents/Documents/fast-gatk/fastgatk-native/build（OpenMP 后端，OMP_NUM_THREADS=4）
方法：对每个工具跑最小 invocation，解析 manifest 中的 ``telemetry.execution_space``。

## 汇总

- ✅ OpenMP kernel 验证：28 个工具
- 🟡 Serial/Host-only（registry 声明，非 OpenMP）：0 个工具
- ⚠️ 跳过（fixture 缺失）：5 个工具
- ❌ 失败：15 个工具

## 详细结果

| Tool | rc | wall(s) | execution_space | kernel_spaces | status |
| --- | ---: | ---: | --- | --- | --- |
| HaplotypeCaller | 0 | 0.021 | OpenMP | execution_space=OpenMP, pairhmm_normalization_execution_space=OpenMP, pairhmm_marginalization_execution_space=OpenMP, pairhmm_uncertainty_execution_space=OpenMP, pairhmm_execution_space=OpenMP, sw_execution_space=OpenMP, graph_haplotype_sw_execution_space=OpenMP, graph_execution_space=OpenMP, activity_execution_space=OpenMP, read_filter_execution_space=OpenMP, reference_confidence_execution_space=OpenMP, call_confidence_execution_space=OpenMP | ✅ OpenMP |
| BaseRecalibrator | 0 | 0.641 | — | kernel_execution_space=OpenMP, covariate_kernel_execution_space=OpenMP | ⚠️ no manifest |
| ApplyBQSR | 2 | 0.009 | — | — | ❌ rc=2 |
| Mutect2 | 0 | 0.092 | — | somatic_likelihood_execution_space=OpenMP, somatic_posterior_execution_space=OpenMP | ⚠️ no manifest |
| FilterMutectCalls | 0 | 0.012 | — | — | ⚠️ no manifest |
| CollectAllelicCounts | 0 | 0.005 | — | — | ⚠️ no manifest |
| CollectF1R2Counts | 0 | 0.015 | — | — | ⚠️ no manifest |
| CollectReadCounts | 0 | 0.011 | — | — | ⚠️ no manifest |
| CombineGVCFs | 2 | 0.005 | — | — | ❌ rc=2 |
| GenotypeGVCFs | 0 | 0.012 | — | genotype_kernel_execution_space=OpenMP, pl_remap_kernel_execution_space=OpenMP, allele_count_kernel_execution_space=OpenMP, cross_sample_reference_execution_space=OpenMP, allele_field_remap_kernel_execution_space=OpenMP, cohort_af_kernel_execution_space=OpenMP, genotype_prior_kernel_execution_space=OpenMP | ⚠️ no manifest |
| ReblockGVCF | 0 | 0.011 | — | — | ⚠️ no manifest |
| SelectVariants | 0 | 0.015 | — | — | ⚠️ no manifest |
| VariantFiltration | 0 | 0.009 | — | — | ⚠️ no manifest |
| LeftAlignAndTrimVariants | 0 | 0.009 | — | — | ⚠️ no manifest |
| ValidateVariants | 0 | 0.003 | — | — | ⚠️ no manifest |
| VariantsToTable | 0 | 0.004 | — | — | ⚠️ no manifest |
| VariantEval | 0 | 0.004 | — | — | ⚠️ no manifest |
| GatherVcfs | 0 | 0.009 | — | — | ⚠️ no manifest |
| GatherTranches | 2 | 0.002 | — | — | ❌ rc=2 |
| SortSam | 0 | 0.012 | — | — | ⚠️ no manifest |
| MarkDuplicates | 0 | 0.013 | — | — | ⚠️ no manifest |
| GatherBQSRReports | 0 | 0.008 | — | — | ⚠️ no manifest |
| AnalyzeCovariates | 2 | 0.007 | — | — | ❌ rc=2 |
| CountReads | 0 | 0.004 | — | — | ⚠️ no manifest |
| FlagStat | 0 | 0.003 | — | — | ⚠️ no manifest |
| SplitIntervals | 0 | 0.003 | — | — | ⚠️ no manifest |
| FilterIntervals | 2 | 0.008 | — | — | ❌ rc=2 |
| PreprocessIntervals | 2 | 0.004 | — | — | ❌ rc=2 |
| AnnotateIntervals | 2 | 0.003 | — | — | ❌ rc=2 |
| CountBasesInReference | 0 | 0.008 | — | — | ⚠️ no manifest |
| CompareReferences | 0 | 0.002 | — | — | ⚠️ no manifest |
| CheckReferenceCompatibility | 2 | 0.003 | — | — | ❌ rc=2 |
| FastaReferenceMaker | 0 | 0.017 | — | — | ⚠️ no manifest |
| FastaAlternateReferenceMaker | 0 | 0.019 | — | — | ⚠️ no manifest |
| ShiftFasta | 2 | 0.002 | — | — | ❌ rc=2 |
| IndexFeatureFile | 0 | 0.006 | — | — | ⚠️ no manifest |
| GetPileupSummaries | 2 | 0.004 | — | — | ❌ rc=2 |
| CalculateContamination | 2 | 0.005 | — | — | ❌ rc=2 |
| LearnReadOrientationModel | 2 | 0.004 | — | — | ❌ rc=2 |
| DepthOfCoverage | 0 | 0.008 | — | — | ⚠️ no manifest |
| ApplyVQSR | 2 | 0.002 | — | — | ❌ rc=2 |
| VariantRecalibrator | 2 | 1.080 | — | — | ❌ rc=2 |
| GenomicsDBImport | 2 | 0.003 | — | — | ❌ rc=2 |
| GatherPileupSummaries | - | - | - | - | ⚠️ missing fixture dependency |
| CallCopyRatioSegments | - | - | - | - | ⚠️ missing fixture dependency |
| ModelSegments | - | - | - | - | ⚠️ missing fixture dependency |
| DenoiseReadCounts | - | - | - | - | ⚠️ missing fixture dependency |
| CreateReadCountPanelOfNormals | - | - | - | - | ⚠️ missing fixture dependency |

详细 JSON：`/home/turing-agents/Documents/fast-gatk/report/openmp_kernel_smoke_20260917-020023/summary.json`