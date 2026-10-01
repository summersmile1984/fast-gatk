# Verify-script re-run index — 2026-09-28

- Tools covered: **50** (excluding `genomicsdb-export` (native-only bridge, no verify script))
- Scripts total: **293**
- Passed: **292**    Failed: **1**    Skipped: **0**
- Total elapsed: **10421.234s**

Per-tool summaries (cross-ref with `fastgatk-native/docs/cli-alignment.md` *Bit-identical / bounded-parity evidence* table):

| Tool | Native binary | Scripts | Passed | Failed | Skipped | Elapsed (s) | Report |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `AnalyzeCovariates` | `fastgatk-analyze-covariates` | 2 | 2 | 0 | 0 | 35.147 | [analyze-covariates-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/analyze-covariates/analyze-covariates-rerun-20260928.md) |
| `AnnotateIntervals` | `fastgatk-annotate-intervals` | 1 | 1 | 0 | 0 | 43.664 | [annotate-intervals-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/annotate-intervals/annotate-intervals-rerun-20260928.md) |
| `ApplyBQSR` | `fastgatk-apply-bqsr` | 1 | 1 | 0 | 0 | 46.614 | [apply-bqsr-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/apply-bqsr/apply-bqsr-rerun-20260928.md) |
| `ApplyVQSR` | `fastgatk-apply-vqsr` | 6 | 6 | 0 | 0 | 187.317 | [apply-vqsr-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/apply-vqsr/apply-vqsr-rerun-20260928.md) |
| `BaseRecalibrator` | `fastgatk-bqsr` | 9 | 9 | 0 | 0 | 432.52 | [bqsr-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/bqsr/bqsr-rerun-20260928.md) |
| `CalculateContamination` | `fastgatk-calculate-contamination` | 2 | 2 | 0 | 0 | 31.306 | [calculate-contamination-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/calculate-contamination/calculate-contamination-rerun-20260928.md) |
| `CallCopyRatioSegments` | `fastgatk-call-copy-ratio-segments` | 5 | 5 | 0 | 0 | 40.602 | [call-copy-ratio-segments-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/call-copy-ratio-segments/call-copy-ratio-segments-rerun-20260928.md) |
| `CheckReferenceCompatibility` | `fastgatk-check-reference-compatibility` | 1 | 1 | 0 | 0 | 3.818 | [check-reference-compatibility-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/check-reference-compatibility/check-reference-compatibility-rerun-20260928.md) |
| `CollectAllelicCounts` | `fastgatk-collect-allelic-counts` | 2 | 2 | 0 | 0 | 38.462 | [collect-allelic-counts-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/collect-allelic-counts/collect-allelic-counts-rerun-20260928.md) |
| `CollectF1R2Counts` | `fastgatk-collect-f1r2-counts` | 1 | 1 | 0 | 0 | 7.804 | [collect-f1r2-counts-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/collect-f1r2-counts/collect-f1r2-counts-rerun-20260928.md) |
| `CollectReadCounts` | `fastgatk-collect-read-counts` | 3 | 3 | 0 | 0 | 45.712 | [collect-read-counts-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/collect-read-counts/collect-read-counts-rerun-20260928.md) |
| `CombineGVCFs` | `fastgatk-combine-gvcfs` | 4 | 4 | 0 | 0 | 87.367 | [combine-gvcfs-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/combine-gvcfs/combine-gvcfs-rerun-20260928.md) |
| `CompareReferences` | `fastgatk-compare-references` | 1 | 1 | 0 | 0 | 3.69 | [compare-references-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/compare-references/compare-references-rerun-20260928.md) |
| `CountBasesInReference` | `fastgatk-count-bases-in-reference` | 1 | 1 | 0 | 0 | 5.141 | [count-bases-in-reference-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/count-bases-in-reference/count-bases-in-reference-rerun-20260928.md) |
| `CountReads` | `fastgatk-count-reads` | 2 | 2 | 0 | 0 | 110.991 | [count-reads-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/count-reads/count-reads-rerun-20260928.md) |
| `CreateReadCountPanelOfNormals` | `fastgatk-create-read-count-panel-of-normals` | 4 | 4 | 0 | 0 | 29.326 | [create-read-count-panel-of-normals-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/create-read-count-panel-of-normals/create-read-count-panel-of-normals-rerun-20260928.md) |
| `DenoiseReadCounts` | `fastgatk-denoise-read-counts` | 4 | 4 | 0 | 0 | 48.196 | [denoise-read-counts-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/denoise-read-counts/denoise-read-counts-rerun-20260928.md) |
| `DepthOfCoverage` | `fastgatk-depth-of-coverage` | 4 | 4 | 0 | 0 | 87.494 | [depth-of-coverage-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/depth-of-coverage/depth-of-coverage-rerun-20260928.md) |
| `FastaAlternateReferenceMaker` | `fastgatk-fasta-alternate-reference-maker` | 1 | 1 | 0 | 0 | 8.111 | [fasta-alternate-reference-maker-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/fasta-alternate-reference-maker/fasta-alternate-reference-maker-rerun-20260928.md) |
| `FastaReferenceMaker` | `fastgatk-fasta-reference-maker` | 1 | 1 | 0 | 0 | 49.216 | [fasta-reference-maker-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/fasta-reference-maker/fasta-reference-maker-rerun-20260928.md) |
| `FilterIntervals` | `fastgatk-filter-intervals` | 1 | 1 | 0 | 0 | 14.057 | [filter-intervals-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/filter-intervals/filter-intervals-rerun-20260928.md) |
| `FilterMutectCalls` | `fastgatk-filter-mutect-calls` | 13 | 13 | 0 | 0 | 1192.591 | [filter-mutect-calls-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/filter-mutect-calls/filter-mutect-calls-rerun-20260928.md) |
| `FlagStat` | `fastgatk-flag-stat` | 1 | 1 | 0 | 0 | 89.862 | [flag-stat-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/flag-stat/flag-stat-rerun-20260928.md) |
| `Funcotator` | `fastgatk-funcotator` | 1 | 1 | 0 | 0 | 12.664 | [funcotator-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/funcotator/funcotator-rerun-20260928.md) |
| `GatherBQSRReports` | `fastgatk-gather-bqsr-reports` | 1 | 1 | 0 | 0 | 32.948 | [gather-bqsr-reports-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/gather-bqsr-reports/gather-bqsr-reports-rerun-20260928.md) |
| `GatherPileupSummaries` | `fastgatk-gather-pileup-summaries` | 2 | 2 | 0 | 0 | 11.147 | [gather-pileup-summaries-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/gather-pileup-summaries/gather-pileup-summaries-rerun-20260928.md) |
| `GatherTranches` | `fastgatk-gather-tranches` | 1 | 1 | 0 | 0 | 8.023 | [gather-tranches-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/gather-tranches/gather-tranches-rerun-20260928.md) |
| `GatherVcfs` | `fastgatk-gather-vcfs` | 2 | 2 | 0 | 0 | 21.317 | [gather-vcfs-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/gather-vcfs/gather-vcfs-rerun-20260928.md) |
| `GenomicsDBImport` | `fastgatk-genomicsdb-import` | 7 | 7 | 0 | 0 | 84.031 | [genomicsdb-import-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/genomicsdb-import/genomicsdb-import-rerun-20260928.md) |
| `GenotypeGVCFs` | `fastgatk-genotype-gvcf` | 25 | 25 | 0 | 0 | 2038.662 | [genotype-gvcf-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/genotype-gvcf/genotype-gvcf-rerun-20260928.md) |
| `GetPileupSummaries` | `fastgatk-get-pileup-summaries` | 2 | 2 | 0 | 0 | 45.021 | [get-pileup-summaries-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/get-pileup-summaries/get-pileup-summaries-rerun-20260928.md) |
| `HaplotypeCaller` | `fastgatk-hc-call` | 74 | 74 | 0 | 0 | 2369.056 | [hc-call-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/hc-call/hc-call-rerun-20260928.md) |
| `IndexFeatureFile` | `fastgatk-index-feature-file` | 1 | 1 | 0 | 0 | 31.923 | [index-feature-file-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/index-feature-file/index-feature-file-rerun-20260928.md) |
| `LearnReadOrientationModel` | `fastgatk-learn-read-orientation-model` | 2 | 2 | 0 | 0 | 20.372 | [learn-read-orientation-model-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/learn-read-orientation-model/learn-read-orientation-model-rerun-20260928.md) |
| `LeftAlignAndTrimVariants` | `fastgatk-left-align-trim` | 4 | 4 | 0 | 0 | 164.203 | [left-align-trim-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/left-align-trim/left-align-trim-rerun-20260928.md) |
| `MarkDuplicates` | `fastgatk-mark-duplicates` | 4 | 4 | 0 | 0 | 53.373 | [mark-duplicates-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/mark-duplicates/mark-duplicates-rerun-20260928.md) |
| `ModelSegments` | `fastgatk-model-segments` | 10 | 10 | 0 | 0 | 67.903 | [model-segments-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/model-segments/model-segments-rerun-20260928.md) |
| `Mutect2` | `fastgatk-mutect2` | 40 | 40 | 0 | 0 | 1512.051 | [mutect2-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/mutect2/mutect2-rerun-20260928.md) |
| `PreprocessIntervals` | `fastgatk-preprocess-intervals` | 1 | 1 | 0 | 0 | 18.099 | [preprocess-intervals-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/preprocess-intervals/preprocess-intervals-rerun-20260928.md) |
| `ReblockGVCF` | `fastgatk-reblock-gvcf` | 6 | 5 | 1 | 0 | 124.169 | [reblock-gvcf-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/reblock-gvcf/reblock-gvcf-rerun-20260928.md) |
| `SelectVariants` | `fastgatk-select-variants` | 6 | 6 | 0 | 0 | 309.704 | [select-variants-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/select-variants/select-variants-rerun-20260928.md) |
| `ShiftFasta` | `fastgatk-shift-fasta` | 1 | 1 | 0 | 0 | 3.669 | [shift-fasta-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/shift-fasta/shift-fasta-rerun-20260928.md) |
| `SortSam` | `fastgatk-sort-sam` | 4 | 4 | 0 | 0 | 27.479 | [sort-sam-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/sort-sam/sort-sam-rerun-20260928.md) |
| `SplitIntervals` | `fastgatk-split-intervals` | 1 | 1 | 0 | 0 | 34.829 | [split-intervals-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/split-intervals/split-intervals-rerun-20260928.md) |
| `ValidateVariants` | `fastgatk-validate-variants` | 3 | 3 | 0 | 0 | 111.394 | [validate-variants-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/validate-variants/validate-variants-rerun-20260928.md) |
| `VariantAnnotator` | `fastgatk-variant-annotator` | 2 | 2 | 0 | 0 | 22.795 | [variant-annotator-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/variant-annotator/variant-annotator-rerun-20260928.md) |
| `VariantEval` | `fastgatk-variant-eval` | 4 | 4 | 0 | 0 | 86.756 | [variant-eval-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/variant-eval/variant-eval-rerun-20260928.md) |
| `VariantFiltration` | `fastgatk-variant-filtration` | 6 | 6 | 0 | 0 | 324.97 | [variant-filtration-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/variant-filtration/variant-filtration-rerun-20260928.md) |
| `VariantRecalibrator` | `fastgatk-variant-recalibrator` | 10 | 10 | 0 | 0 | 161.204 | [variant-recalibrator-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/variant-recalibrator/variant-recalibrator-rerun-20260928.md) |
| `VariantsToTable` | `fastgatk-variants-to-table` | 3 | 3 | 0 | 0 | 84.464 | [variants-to-table-rerun-20260928.md](fastgatk-native/evidence/2026-09-28-rerun/variants-to-table/variants-to-table-rerun-20260928.md) |

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
