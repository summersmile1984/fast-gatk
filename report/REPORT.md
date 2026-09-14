# fast-gatk 总对比报表

生成方式：`aggregate_report.py`（准确率解析自最新回归日志；
效率/资源来自 `benchmark_*.py`）

| 工具 | 准确率(通过/总数) | native wall(s) | java wall(s) | 资源占用 |
| --- | ---: | ---: | ---: | --- |
| AnalyzeCovariates | 6/6 | 0.0142 | — | — |
| AnnotateIntervals | 2/2 | 0.0368 | — | — |
| ApplyBQSR | 2/2 | — | — | — |
| ApplyVQSR | 12/12 | 0.0360 | — | — |
| BaseRecalibrator | 18/18 | — | — | — |
| CalculateContamination | 4/4 | 0.4975 | — | — |
| CallCopyRatioSegments | 10/10 | 0.0924 | — | — |
| CheckReferenceCompatibility | 2/2 | 0.0019 | — | — |
| CollectAllelicCounts | 4/4 | 0.0393 | — | — |
| CollectF1R2Counts | 2/2 | 0.0145 | — | — |
| CollectReadCounts | 4/4 | 0.0443 | — | — |
| CombineGVCFs | 8/8 | 0.0431 | — | — |
| CompareReferences | 2/2 | 0.0238 | — | — |
| CountBasesInReference | 2/2 | 0.0363 | — | — |
| CountReads | 2/2 | — | — | — |
| CreateReadCountPanelOfNormals | 6/6 | 0.0649 | — | — |
| DenoiseReadCounts | 8/8 | 0.0659 | — | — |
| DepthOfCoverage | 8/8 | — | — | — |
| FastaAlternateReferenceMaker | 2/2 | — | — | — |
| FastaReferenceMaker | 2/2 | — | — | — |
| FilterIntervals | 2/2 | 0.0370 | — | — |
| FilterMutectCalls | 34/34 | fail | — | 7860 KB RSS |
| FlagStat | 2/2 | — | — | — |
| GatherBQSRReports | 2/2 | 0.0263 | — | — |
| GatherPileupSummaries | 4/4 | — | — | — |
| GatherTranches | 2/2 | 0.0040 | — | — |
| GatherVcfs | 4/4 | 0.0363 | — | — |
| GenomicsDBImport | 16/16 | — | — | — |
| GenotypeGVCFs | 50/50 | 1.3176 | — | — |
| GetPileupSummaries | 4/4 | 0.0185 | — | — |
| HaplotypeCaller | 116/116 | 0.0733 | 4.0986 | 17476 KB RSS |
| IndexFeatureFile | 2/2 | 0.0054 | — | — |
| LearnReadOrientationModel | 4/4 | 0.0230 | — | — |
| LeftAlignAndTrimVariants | 8/8 | 0.1017 | — | — |
| MarkDuplicates | 8/8 | 0.0286 | — | — |
| ModelSegments | 22/22 | 0.0741 | — | — |
| Mutect2 | 64/64 | 0.1318 | — | 21476 KB RSS |
| PreprocessIntervals | 2/2 | 0.0161 | — | — |
| ReblockGVCF | 12/12 | 0.0325 | — | — |
| SelectVariants | 12/12 | 0.2179 | — | — |
| ShiftFasta | 2/2 | 0.0090 | — | — |
| SortSam | 8/8 | 0.0039 | — | — |
| SplitIntervals | 2/2 | 0.0089 | — | — |
| ValidateVariants | 6/6 | 0.0581 | — | — |
| VariantEval | 8/8 | 0.0066 | — | — |
| VariantFiltration | 12/12 | 0.0284 | — | — |
| VariantRecalibrator | 16/16 | 0.0577 | — | — |
| VariantsToTable | 4/4 | 0.0289 | — | — |

> 注：
> 1. 准确率 = 该工具所有 oracle/契约测试在双后端日志中的通过数/总数（失败标红）。
> 2. native wall = benchmark 脚本报告的 p50/wall 墙钟；java wall 仅在 benchmark 同时跑
>    Java 时才有值（多数 benchmark 只测 native，无 Java 对比）。
> 3. 资源占用目前只有 `benchmark_end_to_end.py`（HC/Mutect2）报告 max RSS + FS 读写；
>    其余工具的逐工具 RSS 未采集（benchmark 只测墙钟）。
> 4. 若要「逐工具 native vs java 的效率 + RSS 对比」，需扩展各 `benchmark_*.py` 统一
>    加 Java 计时与 `/usr/bin/time` 包裹。