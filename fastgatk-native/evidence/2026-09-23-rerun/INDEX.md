# Verify-script re-run index — 2026-09-23

- Tools covered: **48** (excluding `genomicsdb-export` (native-only bridge, no verify script))
- Scripts total: **272**
- Passed: **272**    Failed: **0**    Skipped: **0**
- Total elapsed: **10677.704s**

Per-tool summaries (cross-ref with `fastgatk-native/docs/cli-alignment.md` *Bit-identical / bounded-parity evidence* table):

| Tool | Native binary | Scripts | Passed | Failed | Skipped | Elapsed (s) | Report |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `AnalyzeCovariates` | `fastgatk-analyze-covariates` | 2 | 2 | 0 | 0 | 30.294 | [analyze-covariates-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/analyze-covariates/analyze-covariates-rerun-20260923.md) |
| `AnnotateIntervals` | `fastgatk-annotate-intervals` | 1 | 1 | 0 | 0 | 61.289 | [annotate-intervals-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/annotate-intervals/annotate-intervals-rerun-20260923.md) |
| `ApplyBQSR` | `fastgatk-apply-bqsr` | 1 | 1 | 0 | 0 | 50.253 | [apply-bqsr-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/apply-bqsr/apply-bqsr-rerun-20260923.md) |
| `ApplyVQSR` | `fastgatk-apply-vqsr` | 6 | 6 | 0 | 0 | 268.403 | [apply-vqsr-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/apply-vqsr/apply-vqsr-rerun-20260923.md) |
| `BaseRecalibrator` | `fastgatk-bqsr` | 9 | 9 | 0 | 0 | 576.811 | [bqsr-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/bqsr/bqsr-rerun-20260923.md) |
| `CalculateContamination` | `fastgatk-calculate-contamination` | 2 | 2 | 0 | 0 | 34.792 | [calculate-contamination-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/calculate-contamination/calculate-contamination-rerun-20260923.md) |
| `CallCopyRatioSegments` | `fastgatk-call-copy-ratio-segments` | 5 | 5 | 0 | 0 | 53.955 | [call-copy-ratio-segments-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/call-copy-ratio-segments/call-copy-ratio-segments-rerun-20260923.md) |
| `CheckReferenceCompatibility` | `fastgatk-check-reference-compatibility` | 1 | 1 | 0 | 0 | 5.77 | [check-reference-compatibility-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/check-reference-compatibility/check-reference-compatibility-rerun-20260923.md) |
| `CollectAllelicCounts` | `fastgatk-collect-allelic-counts` | 2 | 2 | 0 | 0 | 46.938 | [collect-allelic-counts-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/collect-allelic-counts/collect-allelic-counts-rerun-20260923.md) |
| `CollectF1R2Counts` | `fastgatk-collect-f1r2-counts` | 1 | 1 | 0 | 0 | 12.183 | [collect-f1r2-counts-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/collect-f1r2-counts/collect-f1r2-counts-rerun-20260923.md) |
| `CollectReadCounts` | `fastgatk-collect-read-counts` | 3 | 3 | 0 | 0 | 68.002 | [collect-read-counts-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/collect-read-counts/collect-read-counts-rerun-20260923.md) |
| `CombineGVCFs` | `fastgatk-combine-gvcfs` | 4 | 4 | 0 | 0 | 118.562 | [combine-gvcfs-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/combine-gvcfs/combine-gvcfs-rerun-20260923.md) |
| `CompareReferences` | `fastgatk-compare-references` | 1 | 1 | 0 | 0 | 10.465 | [compare-references-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/compare-references/compare-references-rerun-20260923.md) |
| `CountBasesInReference` | `fastgatk-count-bases-in-reference` | 1 | 1 | 0 | 0 | 8.472 | [count-bases-in-reference-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/count-bases-in-reference/count-bases-in-reference-rerun-20260923.md) |
| `CountReads` | `fastgatk-count-reads` | 2 | 2 | 0 | 0 | 153.796 | [count-reads-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/count-reads/count-reads-rerun-20260923.md) |
| `CreateReadCountPanelOfNormals` | `fastgatk-create-read-count-panel-of-normals` | 4 | 4 | 0 | 0 | 63.448 | [create-read-count-panel-of-normals-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/create-read-count-panel-of-normals/create-read-count-panel-of-normals-rerun-20260923.md) |
| `DenoiseReadCounts` | `fastgatk-denoise-read-counts` | 4 | 4 | 0 | 0 | 61.985 | [denoise-read-counts-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/denoise-read-counts/denoise-read-counts-rerun-20260923.md) |
| `DepthOfCoverage` | `fastgatk-depth-of-coverage` | 4 | 4 | 0 | 0 | 131.323 | [depth-of-coverage-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/depth-of-coverage/depth-of-coverage-rerun-20260923.md) |
| `FastaAlternateReferenceMaker` | `fastgatk-fasta-alternate-reference-maker` | 1 | 1 | 0 | 0 | 13.385 | [fasta-alternate-reference-maker-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/fasta-alternate-reference-maker/fasta-alternate-reference-maker-rerun-20260923.md) |
| `FastaReferenceMaker` | `fastgatk-fasta-reference-maker` | 1 | 1 | 0 | 0 | 68.202 | [fasta-reference-maker-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/fasta-reference-maker/fasta-reference-maker-rerun-20260923.md) |
| `FilterIntervals` | `fastgatk-filter-intervals` | 1 | 1 | 0 | 0 | 20.663 | [filter-intervals-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/filter-intervals/filter-intervals-rerun-20260923.md) |
| `FilterMutectCalls` | `fastgatk-filter-mutect-calls` | 10 | 10 | 0 | 0 | 206.064 | [filter-mutect-calls-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/filter-mutect-calls/filter-mutect-calls-rerun-20260923.md) |
| `FlagStat` | `fastgatk-flag-stat` | 1 | 1 | 0 | 0 | 111.86 | [flag-stat-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/flag-stat/flag-stat-rerun-20260923.md) |
| `GatherBQSRReports` | `fastgatk-gather-bqsr-reports` | 1 | 1 | 0 | 0 | 61.666 | [gather-bqsr-reports-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/gather-bqsr-reports/gather-bqsr-reports-rerun-20260923.md) |
| `GatherPileupSummaries` | `fastgatk-gather-pileup-summaries` | 2 | 2 | 0 | 0 | 15.816 | [gather-pileup-summaries-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/gather-pileup-summaries/gather-pileup-summaries-rerun-20260923.md) |
| `GatherTranches` | `fastgatk-gather-tranches` | 1 | 1 | 0 | 0 | 16.108 | [gather-tranches-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/gather-tranches/gather-tranches-rerun-20260923.md) |
| `GatherVcfs` | `fastgatk-gather-vcfs` | 2 | 2 | 0 | 0 | 38.56 | [gather-vcfs-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/gather-vcfs/gather-vcfs-rerun-20260923.md) |
| `GenomicsDBImport` | `fastgatk-genomicsdb-import` | 7 | 7 | 0 | 0 | 128.648 | [genomicsdb-import-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/genomicsdb-import/genomicsdb-import-rerun-20260923.md) |
| `GenotypeGVCFs` | `fastgatk-genotype-gvcf` | 25 | 25 | 0 | 0 | 2138.171 | [genotype-gvcf-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/genotype-gvcf/genotype-gvcf-rerun-20260923.md) |
| `GetPileupSummaries` | `fastgatk-get-pileup-summaries` | 2 | 2 | 0 | 0 | 81.738 | [get-pileup-summaries-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/get-pileup-summaries/get-pileup-summaries-rerun-20260923.md) |
| `HaplotypeCaller` | `fastgatk-hc-call` | 61 | 61 | 0 | 0 | 2623.018 | [hc-call-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/hc-call/hc-call-rerun-20260923.md) |
| `IndexFeatureFile` | `fastgatk-index-feature-file` | 1 | 1 | 0 | 0 | 105.244 | [index-feature-file-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/index-feature-file/index-feature-file-rerun-20260923.md) |
| `LearnReadOrientationModel` | `fastgatk-learn-read-orientation-model` | 2 | 2 | 0 | 0 | 30.775 | [learn-read-orientation-model-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/learn-read-orientation-model/learn-read-orientation-model-rerun-20260923.md) |
| `LeftAlignAndTrimVariants` | `fastgatk-left-align-trim` | 4 | 4 | 0 | 0 | 199.623 | [left-align-trim-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/left-align-trim/left-align-trim-rerun-20260923.md) |
| `MarkDuplicates` | `fastgatk-mark-duplicates` | 4 | 4 | 0 | 0 | 69.859 | [mark-duplicates-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/mark-duplicates/mark-duplicates-rerun-20260923.md) |
| `ModelSegments` | `fastgatk-model-segments` | 10 | 10 | 0 | 0 | 95.668 | [model-segments-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/model-segments/model-segments-rerun-20260923.md) |
| `Mutect2` | `fastgatk-mutect2` | 38 | 38 | 0 | 0 | 1161.277 | [mutect2-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/mutect2/mutect2-rerun-20260923.md) |
| `PreprocessIntervals` | `fastgatk-preprocess-intervals` | 1 | 1 | 0 | 0 | 27.16 | [preprocess-intervals-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/preprocess-intervals/preprocess-intervals-rerun-20260923.md) |
| `ReblockGVCF` | `fastgatk-reblock-gvcf` | 6 | 6 | 0 | 0 | 172.245 | [reblock-gvcf-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/reblock-gvcf/reblock-gvcf-rerun-20260923.md) |
| `SelectVariants` | `fastgatk-select-variants` | 6 | 6 | 0 | 0 | 423.558 | [select-variants-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/select-variants/select-variants-rerun-20260923.md) |
| `ShiftFasta` | `fastgatk-shift-fasta` | 1 | 1 | 0 | 0 | 5.814 | [shift-fasta-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/shift-fasta/shift-fasta-rerun-20260923.md) |
| `SortSam` | `fastgatk-sort-sam` | 4 | 4 | 0 | 0 | 41.845 | [sort-sam-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/sort-sam/sort-sam-rerun-20260923.md) |
| `SplitIntervals` | `fastgatk-split-intervals` | 1 | 1 | 0 | 0 | 53.65 | [split-intervals-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/split-intervals/split-intervals-rerun-20260923.md) |
| `ValidateVariants` | `fastgatk-validate-variants` | 3 | 3 | 0 | 0 | 141.143 | [validate-variants-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/validate-variants/validate-variants-rerun-20260923.md) |
| `VariantEval` | `fastgatk-variant-eval` | 4 | 4 | 0 | 0 | 127.712 | [variant-eval-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/variant-eval/variant-eval-rerun-20260923.md) |
| `VariantFiltration` | `fastgatk-variant-filtration` | 6 | 6 | 0 | 0 | 377.549 | [variant-filtration-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/variant-filtration/variant-filtration-rerun-20260923.md) |
| `VariantRecalibrator` | `fastgatk-variant-recalibrator` | 10 | 10 | 0 | 0 | 238.426 | [variant-recalibrator-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/variant-recalibrator/variant-recalibrator-rerun-20260923.md) |
| `VariantsToTable` | `fastgatk-variants-to-table` | 3 | 3 | 0 | 0 | 125.516 | [variants-to-table-rerun-20260923.md](fastgatk-native/evidence/2026-09-23-rerun/variants-to-table/variants-to-table-rerun-20260923.md) |

## Bridge binary with no verify script

- `fastgatk-genomicsdb-export` — (native-only bridge, no verify script); cross-reference the docs.

## GATK-alignment status

Run `aggregate_alignment_status.py` after a re-run to cross-reference
this index with the *Bit-identical / bounded-parity evidence* table
in `fastgatk-native/docs/cli-alignment.md`.  It writes
`ALIGNMENT_STATUS.md` next to this file.

## Native-vs-Java resource comparison

Run `compare_native_vs_java.py` to classify every descendant process
of the rerun into `native` (comm starts with `fastgatk-`) vs `java`
(comm == `java`) buckets, and emit `NATIVE_VS_JAVA.md` next to this
file.  Anomalies (native slower than Java, or native RSS exceeding
Java's) are auto-flagged in that report.

## Workflow documentation

Full workflow, schema details, and an anomaly follow-up playbook
live in `fastgatk-native/docs/regression-evidence.md`.

## Reproduction

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --repo .
```

```bash
python3 fastgatk-native/scripts/verify_rerun_report.py --repo .
```

```bash
python3 fastgatk-native/scripts/aggregate_alignment_status.py --repo .
```

```bash
python3 fastgatk-native/scripts/compare_native_vs_java.py --repo .
```
