# CLI Alignment — fastgatk-native ↔ GATK 4.6.2.0

**Slug:** `cli-alignment-doc`

## How to regenerate

```
python3 fastgatk-native/scripts/generate_cli_alignment_doc.py
```

The script reads only:

- `fastgatk-native/CMakeLists.txt` — binary list (`add_executable(fastgatk-*)`).
- `fastgatk-native/src/<tool>.cpp` — argv-parser block (regex on `is_option(argument, "--name")` and `argument == "--name"`).
- `fastgatk-native/include/fastgatk/cli/gatk_cli_catalog.hpp` — per-tool catalog the native `fastgatk::cli::consume(...)` function consults.
- `third_party/gatk-package/gatk-4.6.2.0/gatkdoc/<class>.json` — GATK 4.6.2.0 `arguments[]` schema for the matching Java tool.

Data sources are pinned: **GATK 4.6.2.0** (Java + JSON), **Kokkos 5.2.0** for native.  This document reflects the current `CMakeLists.txt` and `src/*.cpp`; it does not make claims about numerical/algorithmic parity, only the CLI surface.

## Status enum

- `accepted` — native argv parser handles the option literal.
- `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is silently consumed by `fastgatk::cli::consume(...)`.
- `unsupported` — GATK exposes the option but the native tool does not handle it (will likely fail or be ignored).
- `native_only` — native tool accepts the option literal but GATK 4.6.2.0 has no such argument (custom native switch).

Total entries: **51** (production binaries; smoke/test/oracle binaries live in `fastgatk-native/tests/` and are intentionally excluded).

## Overview

| Native binary | GATK Java class | GATK reference JSON | Source cpp | Count of `accepted` | Status |
| --- | --- | --- | --- | --- | --- |
| `fastgatk-genomicsdb-export` | `GenomicsDBExport (native-only bridge binary)` | `<none>` | `src/genomicsdb_export_bridge.cpp` | 0 | native-only bridge |
| `fastgatk-hc-call` | `HaplotypeCaller` | `org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.json` | `src/hc_call.cpp` | 91 | `accepted`=91, `unsupported`=17, `ignored_via_catalog`=66 |
| `fastgatk-bqsr` | `BaseRecalibrator` | `org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.json` | `src/bqsr_tool.cpp` | 19 | `accepted`=19, `unsupported`=26, `ignored_via_catalog`=94 |
| `fastgatk-apply-bqsr` | `ApplyBQSR` | `org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.json` | `src/bqsr_tool.cpp` | 20 | `accepted`=20, `unsupported`=38, `ignored_via_catalog`=113 |
| `fastgatk-genotype-gvcf` | `GenotypeGVCFs` | `org_broadinstitute_hellbender_tools_walkers_GenotypeGVCFs.json` | `src/genotype_gvcf_tool.cpp` | 24 | `accepted`=24, `unsupported`=55, `ignored_via_catalog`=153 |
| `fastgatk-reblock-gvcf` | `ReblockGVCF` | `org_broadinstitute_hellbender_tools_walkers_variantutils_ReblockGVCF.json` | `src/reblock_gvcf_tool.cpp` | 22 | `accepted`=22, `unsupported`=76, `ignored_via_catalog`=189 |
| `fastgatk-select-variants` | `SelectVariants` | `org_broadinstitute_hellbender_tools_walkers_variantutils_SelectVariants.json` | `src/select_variants_tool.cpp` | 28 | `accepted`=28, `unsupported`=90, `ignored_via_catalog`=235 |
| `fastgatk-gather-vcfs` | `GatherVcfs (Picard)` | `picard_vcf_GatherVcfs.json` | `src/gather_vcfs_tool.cpp` | 10 | `accepted`=10, `unsupported`=99, `ignored_via_catalog`=235 |
| `fastgatk-left-align-trim` | `LeftAlignAndTrimVariants` | `org_broadinstitute_hellbender_tools_walkers_variantutils_LeftAlignAndTrimVariants.json` | `src/left_align_tool.cpp` | 20 | `accepted`=20, `unsupported`=112, `ignored_via_catalog`=251 |
| `fastgatk-funcotator` | `Funcotator` | `org_broadinstitute_hellbender_tools_funcotator_Funcotator.json` | `src/funcotator_tool.cpp` | 10 | `accepted`=10, `unsupported`=163, `ignored_via_catalog`=251 |
| `fastgatk-variant-annotator` | `VariantAnnotator` | `org_broadinstitute_hellbender_tools_walkers_annotator_VariantAnnotator.json` | `src/variant_annotator_tool.cpp` | 12 | `accepted`=12, `unsupported`=209, `ignored_via_catalog`=251 |
| `fastgatk-variant-filtration` | `VariantFiltration` | `org_broadinstitute_hellbender_tools_walkers_filters_VariantFiltration.json` | `src/variant_filtration_tool.cpp` | 27 | `accepted`=27, `unsupported`=222, `ignored_via_catalog`=272 |
| `fastgatk-apply-vqsr` | `ApplyVQSR` | `org_broadinstitute_hellbender_tools_walkers_vqsr_ApplyVQSR.json` | `src/apply_vqsr_tool.cpp` | 13 | `accepted`=13, `unsupported`=236, `ignored_via_catalog`=298 |
| `fastgatk-variant-recalibrator` | `VariantRecalibrator` | `org_broadinstitute_hellbender_tools_walkers_vqsr_VariantRecalibrator.json` | `src/variant_recalibrator_tool.cpp` | 14 | `accepted`=14, `unsupported`=250, `ignored_via_catalog`=344 |
| `fastgatk-sort-sam` | `SortSam (Picard)` | `picard_sam_SortSam.json` | `src/sort_sam_tool.cpp` | 14 | `accepted`=14, `unsupported`=254, `ignored_via_catalog`=344 |
| `fastgatk-mark-duplicates` | `MarkDuplicates (Picard)` | `picard_sam_markduplicates_MarkDuplicates.json` | `src/mark_duplicates_tool.cpp` | 8 | `accepted`=8, `unsupported`=298, `ignored_via_catalog`=344 |
| `fastgatk-mutect2` | `Mutect2` | `org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.json` | `src/mutect2_tool.cpp` | 70 | `accepted`=70, `unsupported`=316, `ignored_via_catalog`=427 |
| `fastgatk-combine-gvcfs` | `CombineGVCFs` | `org_broadinstitute_hellbender_tools_walkers_CombineGVCFs.json` | `src/combine_gvcf_tool.cpp` | 14 | `accepted`=14, `unsupported`=337, `ignored_via_catalog`=454 |
| `fastgatk-filter-mutect-calls` | `FilterMutectCalls` | `org_broadinstitute_hellbender_tools_walkers_mutect_filtering_FilterMutectCalls.json` | `src/filter_mutect_tool.cpp` | 42 | `accepted`=42, `unsupported`=350, `ignored_via_catalog`=474 |
| `fastgatk-genomicsdb-import` | `GenomicsDBImport` | `org_broadinstitute_hellbender_tools_genomicsdb_GenomicsDBImport.json` | `src/genomicsdb_import_tool.cpp` | 11 | `accepted`=11, `unsupported`=366, `ignored_via_catalog`=508 |
| `fastgatk-gather-bqsr-reports` | `GatherBQSRReports` | `org_broadinstitute_hellbender_tools_walkers_bqsr_GatherBQSRReports.json` | `src/gather_bqsr_tool.cpp` | 8 | `accepted`=8, `unsupported`=372, `ignored_via_catalog`=508 |
| `fastgatk-analyze-covariates` | `AnalyzeCovariates` | `org_broadinstitute_hellbender_tools_walkers_bqsr_AnalyzeCovariates.json` | `src/analyze_covariates_tool.cpp` | 8 | `accepted`=8, `unsupported`=378, `ignored_via_catalog`=512 |
| `fastgatk-variants-to-table` | `VariantsToTable` | `org_broadinstitute_hellbender_tools_walkers_variantutils_VariantsToTable.json` | `src/variants_to_table_tool.cpp` | 22 | `accepted`=22, `unsupported`=391, `ignored_via_catalog`=529 |
| `fastgatk-variant-eval` | `VariantEval` | `org_broadinstitute_hellbender_tools_walkers_varianteval_VariantEval.json` | `src/variant_eval_tool.cpp` | 21 | `accepted`=21, `unsupported`=400, `ignored_via_catalog`=566 |
| `fastgatk-validate-variants` | `ValidateVariants` | `org_broadinstitute_hellbender_tools_walkers_variantutils_ValidateVariants.json` | `src/validate_variants_tool.cpp` | 16 | `accepted`=16, `unsupported`=413, `ignored_via_catalog`=586 |
| `fastgatk-get-pileup-summaries` | `GetPileupSummaries` | `org_broadinstitute_hellbender_tools_walkers_contamination_GetPileupSummaries.json` | `src/get_pileup_summaries_tool.cpp` | 20 | `accepted`=20, `unsupported`=421, `ignored_via_catalog`=604 |
| `fastgatk-calculate-contamination` | `CalculateContamination` | `org_broadinstitute_hellbender_tools_walkers_contamination_CalculateContamination.json` | `src/calculate_contamination_tool.cpp` | 9 | `accepted`=9, `unsupported`=427, `ignored_via_catalog`=607 |
| `fastgatk-gather-pileup-summaries` | `GatherPileupSummaries` | `<none>` | `src/gather_pileup_summaries_tool.cpp` | 0 | native-only bridge |
| `fastgatk-learn-read-orientation-model` | `LearnReadOrientationModel` | `org_broadinstitute_hellbender_tools_walkers_readorientation_LearnReadOrientationModel.json` | `src/learn_read_orientation_model_tool.cpp` | 11 | `accepted`=11, `unsupported`=433, `ignored_via_catalog`=607 |
| `fastgatk-denoise-read-counts` | `DenoiseReadCounts` | `org_broadinstitute_hellbender_tools_copynumber_DenoiseReadCounts.json` | `src/denoise_read_counts_tool.cpp` | 3 | `accepted`=3, `unsupported`=439, `ignored_via_catalog`=616 |
| `fastgatk-create-read-count-panel-of-normals` | `CreateReadCountPanelOfNormals` | `org_broadinstitute_hellbender_tools_copynumber_CreateReadCountPanelOfNormals.json` | `src/create_read_count_panel_of_normals_tool.cpp` | 4 | `accepted`=4, `unsupported`=447, `ignored_via_catalog`=631 |
| `fastgatk-call-copy-ratio-segments` | `CallCopyRatioSegments` | `org_broadinstitute_hellbender_tools_copynumber_CallCopyRatioSegments.json` | `src/call_copy_ratio_segments_tool.cpp` | 4 | `accepted`=4, `unsupported`=453, `ignored_via_catalog`=639 |
| `fastgatk-model-segments` | `ModelSegments` | `org_broadinstitute_hellbender_tools_copynumber_ModelSegments.json` | `src/model_segments_tool.cpp` | 3 | `accepted`=3, `unsupported`=459, `ignored_via_catalog`=668 |
| `fastgatk-gather-tranches` | `GatherTranches` | `org_broadinstitute_hellbender_tools_walkers_vqsr_GatherTranches.json` | `src/gather_tranches_tool.cpp` | 5 | `accepted`=5, `unsupported`=465, `ignored_via_catalog`=673 |
| `fastgatk-annotate-intervals` | `AnnotateIntervals` | `org_broadinstitute_hellbender_tools_copynumber_AnnotateIntervals.json` | `src/annotate_intervals_tool.cpp` | 6 | `accepted`=6, `unsupported`=479, `ignored_via_catalog`=698 |
| `fastgatk-count-bases-in-reference` | `CountBasesInReference` | `org_broadinstitute_hellbender_tools_walkers_fasta_CountBasesInReference.json` | `src/count_bases_in_reference_tool.cpp` | 5 | `accepted`=5, `unsupported`=493, `ignored_via_catalog`=721 |
| `fastgatk-compare-references` | `CompareReferences` | `org_broadinstitute_hellbender_tools_reference_CompareReferences.json` | `src/compare_references_tool.cpp` | 10 | `accepted`=10, `unsupported`=508, `ignored_via_catalog`=744 |
| `fastgatk-check-reference-compatibility` | `CheckReferenceCompatibility` | `org_broadinstitute_hellbender_tools_reference_CheckReferenceCompatibility.json` | `src/check_reference_compatibility_tool.cpp` | 6 | `accepted`=6, `unsupported`=522, `ignored_via_catalog`=767 |
| `fastgatk-fasta-reference-maker` | `FastaReferenceMaker` | `org_broadinstitute_hellbender_tools_walkers_fasta_FastaReferenceMaker.json` | `src/fasta_reference_tool.cpp` | 6 | `accepted`=6, `unsupported`=536, `ignored_via_catalog`=790 |
| `fastgatk-fasta-alternate-reference-maker` | `FastaAlternateReferenceMaker` | `org_broadinstitute_hellbender_tools_walkers_fasta_FastaAlternateReferenceMaker.json` | `src/fasta_reference_tool.cpp` | 10 | `accepted`=10, `unsupported`=550, `ignored_via_catalog`=813 |
| `fastgatk-shift-fasta` | `ShiftFasta` | `org_broadinstitute_hellbender_tools_walkers_fasta_ShiftFasta.json` | `src/shift_fasta_tool.cpp` | 8 | `accepted`=8, `unsupported`=565, `ignored_via_catalog`=836 |
| `fastgatk-index-feature-file` | `IndexFeatureFile` | `org_broadinstitute_hellbender_tools_IndexFeatureFile.json` | `src/index_feature_file_tool.cpp` | 4 | `accepted`=4, `unsupported`=571, `ignored_via_catalog`=840 |
| `fastgatk-count-reads` | `CountReads` | `org_broadinstitute_hellbender_tools_CountReads.json` | `src/read_metrics_tool.cpp` | 15 | `accepted`=15, `unsupported`=580, `ignored_via_catalog`=858 |
| `fastgatk-flag-stat` | `FlagStat` | `org_broadinstitute_hellbender_tools_FlagStat.json` | `src/read_metrics_tool.cpp` | 15 | `accepted`=15, `unsupported`=589, `ignored_via_catalog`=876 |
| `fastgatk-split-intervals` | `SplitIntervals` | `org_broadinstitute_hellbender_tools_walkers_SplitIntervals.json` | `src/split_intervals_tool.cpp` | 16 | `accepted`=16, `unsupported`=602, `ignored_via_catalog`=896 |
| `fastgatk-filter-intervals` | `FilterIntervals` | `org_broadinstitute_hellbender_tools_copynumber_FilterIntervals.json` | `src/filter_intervals_tool.cpp` | 10 | `accepted`=10, `unsupported`=608, `ignored_via_catalog`=912 |
| `fastgatk-preprocess-intervals` | `PreprocessIntervals` | `org_broadinstitute_hellbender_tools_copynumber_PreprocessIntervals.json` | `src/preprocess_intervals_tool.cpp` | 6 | `accepted`=6, `unsupported`=622, `ignored_via_catalog`=936 |
| `fastgatk-collect-read-counts` | `CollectReadCounts` | `org_broadinstitute_hellbender_tools_copynumber_CollectReadCounts.json` | `src/collect_read_counts_tool.cpp` | 20 | `accepted`=20, `unsupported`=633, `ignored_via_catalog`=948 |
| `fastgatk-collect-f1r2-counts` | `CollectF1R2Counts` | `<none>` | `src/collect_f1r2_counts_tool.cpp` | 0 | native-only bridge |
| `fastgatk-collect-allelic-counts` | `CollectAllelicCounts` | `org_broadinstitute_hellbender_tools_copynumber_CollectAllelicCounts.json` | `src/collect_allelic_counts_tool.cpp` | 15 | `accepted`=15, `unsupported`=642, `ignored_via_catalog`=968 |
| `fastgatk-depth-of-coverage` | `DepthOfCoverage` | `org_broadinstitute_hellbender_tools_walkers_coverage_DepthOfCoverage.json` | `src/depth_of_coverage_tool.cpp` | 19 | `accepted`=19, `unsupported`=654, `ignored_via_catalog`=999 |


## Bit-identical / bounded-parity evidence

This section cross-references the CLI alignment table with the tool-dispatcher registry and the pinned GATK 4.6.2.0 oracle scripts. It is auto-generated from:

- `fastgatk-native/dispatcher/tool_registry.json` — per-tool status (`contract-compatible` / `adapter`), determinism, fallback policy.
- `fastgatk-native/scripts/verify_*.py` — `byte-identical` and `bit-identical` assertions against pinned GATK output.
- `third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar` — pinned GATK Java runtime invoked by every oracle script.
- `fastgatk-native/tests/pinned_fixture_digests.sha256` — SHA-256 corpus that pins the exact bytes of every test fixture.

### Pinned artifacts

| Artifact | Status | Identifier |
| --- | --- | --- |
| Pinned GATK 4.6.2.0 jar | present (406.0 MiB) | `gatk-package-4.6.2.0-local.jar` |
| SHA-256 of pinned jar | computed | `40b494b1e356931c9143be09c57ba6dfbb8d150167361aac2ceec24089b75667` |
| Pinned fixture digests | present (12 entries) | `pinned_fixture_digests.sha256` |

### Per-binary parity evidence

| Native binary | GATK class | registry status | determinism | oracle scripts | bit-identical scripts | bounded oracles | sample bit-identical script |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `fastgatk-analyze-covariates` | `AnalyzeCovariates` | `contract-compatible` | `strict` | 2 | 1 | 1 | `verify_analyze_covariates.py` |
| `fastgatk-annotate-intervals` | `AnnotateIntervals` | `contract-compatible` | `strict` | 1 | 0 | 0 | `—` |
| `fastgatk-apply-bqsr` | `ApplyBQSR` | `contract-compatible` | `strict` | 1 | 0 | 1 | `—` |
| `fastgatk-apply-vqsr` | `ApplyVQSR` | `contract-compatible` | `strict` | 6 | 0 | 5 | `—` |
| `fastgatk-bqsr` | `BaseRecalibrator` | `contract-compatible` | `strict` | 12 | 0 | 11 | `—` |
| `fastgatk-calculate-contamination` | `CalculateContamination` | `contract-compatible` | `strict` | 2 | 0 | 1 | `—` |
| `fastgatk-call-copy-ratio-segments` | `CallCopyRatioSegments` | `contract-compatible` | `strict` | 5 | 0 | 4 | `—` |
| `fastgatk-check-reference-compatibility` | `CheckReferenceCompatibility` | `contract-compatible` | `strict` | 1 | 0 | 0 | `—` |
| `fastgatk-collect-allelic-counts` | `CollectAllelicCounts` | `contract-compatible` | `strict` | 2 | 0 | 1 | `—` |
| `fastgatk-collect-f1r2-counts` | `CollectF1R2Counts` | `contract-compatible` | `strict` | 1 | 0 | 0 | `—` |
| `fastgatk-collect-read-counts` | `CollectReadCounts` | `contract-compatible` | `strict` | 2 | 0 | 1 | `—` |
| `fastgatk-combine-gvcfs` | `CombineGVCFs` | `contract-compatible` | `strict` | 4 | 1 | 3 | `verify_combine_gvcfs.py` |
| `fastgatk-compare-references` | `CompareReferences` | `contract-compatible` | `strict` | 1 | 0 | 0 | `—` |
| `fastgatk-count-bases-in-reference` | `CountBasesInReference` | `contract-compatible` | `strict` | 1 | 0 | 0 | `—` |
| `fastgatk-count-reads` | `CountReads` | `contract-compatible` | `strict` | 1 | 0 | 0 | `—` |
| `fastgatk-create-read-count-panel-of-normals` | `CreateReadCountPanelOfNormals` | `contract-compatible` | `strict` | 3 | 0 | 2 | `—` |
| `fastgatk-denoise-read-counts` | `DenoiseReadCounts` | `contract-compatible` | `strict` | 4 | 1 | 3 | `verify_denoise_read_counts.py` |
| `fastgatk-depth-of-coverage` | `DepthOfCoverage` | `contract-compatible` | `strict` | 4 | 0 | 2 | `—` |
| `fastgatk-fasta-alternate-reference-maker` | `FastaAlternateReferenceMaker` | `contract-compatible` | `strict` | 1 | 0 | 1 | `—` |
| `fastgatk-fasta-reference-maker` | `FastaReferenceMaker` | `contract-compatible` | `strict` | 1 | 0 | 0 | `—` |
| `fastgatk-filter-intervals` | `FilterIntervals` | `contract-compatible` | `strict` | 1 | 0 | 0 | `—` |
| `fastgatk-filter-mutect-calls` | `FilterMutectCalls` | `contract-compatible` | `strict` | 13 | 1 | 12 | `verify_filter_mutect_contamination_joint_oracle.py` |
| `fastgatk-flag-stat` | `FlagStat` | `contract-compatible` | `strict` | 1 | 0 | 0 | `—` |
| `fastgatk-funcotator` | `Funcotator` | `contract-compatible` | `strict` | 1 | 1 | 1 | `verify_funcotator_vcf_gatk_oracle.py` |
| `fastgatk-gather-bqsr-reports` | `GatherBQSRReports` | `contract-compatible` | `strict` | 0 | 0 | 0 | `—` |
| `fastgatk-gather-pileup-summaries` | `GatherPileupSummaries` | `contract-compatible` | `strict` | 2 | 0 | 1 | `—` |
| `fastgatk-gather-tranches` | `GatherTranches` | `contract-compatible` | `strict` | 1 | 0 | 0 | `—` |
| `fastgatk-gather-vcfs` | `GatherVcfs` | `contract-compatible` | `strict` | 2 | 0 | 1 | `—` |
| `fastgatk-genomicsdb-export` | `GenomicsDBExport` | `—` | `—` | 0 | 0 | 0 | `—` |
| `fastgatk-genomicsdb-import` | `GenomicsDBImport` | `adapter` | `strict` | 7 | 0 | 4 | `—` |
| `fastgatk-genotype-gvcf` | `GenotypeGVCFs` | `contract-compatible` | `strict` | 25 | 9 | 19 | `verify_gatk_genotype_gvcf.py` |
| `fastgatk-get-pileup-summaries` | `GetPileupSummaries` | `contract-compatible` | `strict` | 1 | 0 | 0 | `—` |
| `fastgatk-hc-call` | `HaplotypeCaller` | `contract-compatible` | `strict, fast` | 68 | 28 | 55 | `verify_hc_af_zero_format_gatk_oracle.py` |
| `fastgatk-index-feature-file` | `IndexFeatureFile` | `contract-compatible` | `strict` | 1 | 0 | 0 | `—` |
| `fastgatk-learn-read-orientation-model` | `LearnReadOrientationModel` | `contract-compatible` | `strict` | 2 | 0 | 1 | `—` |
| `fastgatk-left-align-trim` | `LeftAlignAndTrimVariants` | `contract-compatible` | `strict` | 4 | 0 | 3 | `—` |
| `fastgatk-mark-duplicates` | `MarkDuplicates` | `contract-compatible` | `strict` | 4 | 0 | 3 | `—` |
| `fastgatk-model-segments` | `ModelSegments` | `contract-compatible` | `strict` | 10 | 1 | 9 | `verify_model_segments_copy_ratio_conditionals_gatk_oracle.py` |
| `fastgatk-mutect2` | `Mutect2` | `contract-compatible` | `strict` | 33 | 3 | 19 | `verify_mutect2_as_value_byte_equality_oracle.py` |
| `fastgatk-preprocess-intervals` | `PreprocessIntervals` | `contract-compatible` | `strict` | 1 | 0 | 0 | `—` |
| `fastgatk-reblock-gvcf` | `ReblockGVCF` | `contract-compatible` | `strict` | 3 | 2 | 2 | `verify_reblock_gvcf_droplowqual_gatk_oracle.py` |
| `fastgatk-select-variants` | `SelectVariants` | `contract-compatible` | `strict` | 6 | 1 | 5 | `verify_select_variants_refonly_gatk_oracle.py` |
| `fastgatk-shift-fasta` | `ShiftFasta` | `contract-compatible` | `strict` | 1 | 0 | 0 | `—` |
| `fastgatk-sort-sam` | `SortSam` | `contract-compatible` | `strict` | 4 | 1 | 3 | `verify_sort_sam_cli_boundary_gatk_oracle.py` |
| `fastgatk-split-intervals` | `SplitIntervals` | `contract-compatible` | `strict` | 1 | 0 | 0 | `—` |
| `fastgatk-validate-variants` | `ValidateVariants` | `contract-compatible` | `strict` | 3 | 0 | 2 | `—` |
| `fastgatk-variant-annotator` | `VariantAnnotator` | `contract-compatible` | `strict` | 2 | 2 | 2 | `verify_variant_annotator_coverage_gatk_oracle.py` |
| `fastgatk-variant-eval` | `VariantEval` | `contract-compatible` | `strict` | 4 | 0 | 3 | `—` |
| `fastgatk-variant-filtration` | `VariantFiltration` | `contract-compatible` | `strict` | 6 | 2 | 5 | `verify_variant_filtration_asfilterstatus_gatk_oracle.py` |
| `fastgatk-variant-recalibrator` | `VariantRecalibrator` | `contract-compatible` | `strict` | 8 | 3 | 7 | `verify_variant_recalibrator.py` |
| `fastgatk-variants-to-table` | `VariantsToTable` | `contract-compatible` | `strict` | 2 | 1 | 1 | `verify_variants_to_table_gatk_oracle.py` |

**Aggregate:** 49 tools have at least one oracle verify script in `fastgatk-native/scripts/`; 272 oracle scripts in total; 58 of them assert byte/bit-identical output against pinned GATK 4.6.2.0; 194 carry the bounded-oracle naming convention (`*_gatk_oracle.py`).

**Dispatcher status counts:** 49 `contract-compatible`, 1 `adapter`.  All 50 tools declare `determinism: ["strict", "fast"]` and `fallback_policy: "explicit"`.


## Per-tool details

## `fastgatk-genomicsdb-export` ↔ `GenomicsDBExport (native-only bridge binary)`

- Source: `fastgatk-native/src/genomicsdb_export_bridge.cpp`
- GATK reference: `<none>` (native-only bridge binary)
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

No gatkdoc JSON is shipped with this GATK version; the tool is native-only.

## `fastgatk-hc-call` ↔ `HaplotypeCaller`

- Source: `fastgatk-native/src/hc_call.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--active-probability-threshold` | `` | `double` | `0.002` | `no` | `advanced` | 1091 | `accepted` |
| `--adaptive-pruning` | `` | `boolean` | `false` | `no` | `advanced` | 1486 | `accepted` |
| `--adaptive-pruning-initial-error-rate` | `` | `double` | `0.001` | `no` | `advanced` | 1496 | `accepted` |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | 1697 | `accepted` |
| `--all-site-pls` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--allele-informative-reads-overlap-margin` | `` | `int` | `2` | `no` | `advanced` | 1443 | `accepted` |
| `--alleles` | `` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | 993 | `accepted` |
| `--allow-non-unique-kmers-in-ref` | `` | `boolean` | `false` | `no` | `advanced` | 1525 | `accepted` |
| `--annotate-with-num-discovered-alleles` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--annotation` | `-A` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--annotation-group` | `-G` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--annotations-to-exclude` | `-AX` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--apply-bqd` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--apply-frd` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--assembly-region-out` | `` | `String` | `null` | `no` | `optional` | 1051 | `accepted` |
| `--assembly-region-padding` | `` | `int` | `100` | `no` | `optional` | 1100 | `accepted` |
| `--bam-output` | `-bamout` | `String` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--bam-writer-type` | `` | `WriterType` | `CALLED_HAPLOTYPES` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--base-quality-score-threshold` | `` | `byte` | `18` | `no` | `optional` | 1201 | `accepted` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--comparison` | `-comp` | `List[FeatureInput[VariantContext]]` | `[]` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--contamination-fraction-per-sample-file` | `-contamination-file` | `File` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--contamination-fraction-to-filter` | `-contamination` | `double` | `0.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | 1690 | `accepted` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--dbsnp` | `-D` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--debug-assembly` | `-debug` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-cap-base-qualities-to-map-quality` | `` | `boolean` | `false` | `no` | `advanced` | 1209 | `accepted` |
| `--disable-optimizations` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | 1672 | `accepted` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | 1719 | `accepted` |
| `--disable-spanning-event-genotyping` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--disable-symmetric-hmm-normalizing` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--disable-tool-default-annotations` | `-disable-tool-default-annotations` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | 1682 | `accepted` |
| `--do-not-correct-overlapping-quality` | `` | `boolean` | `false` | `no` | `advanced` | 1158 | `accepted` |
| `--do-not-run-physical-phasing` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--dont-increase-kmer-sizes-for-cycles` | `` | `boolean` | `false` | `no` | `advanced` | 1463 | `accepted` |
| `--dont-use-dragstr-pair-hmm-scores` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--dont-use-dragstr-priors` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--dont-use-soft-clipped-bases` | `` | `boolean` | `false` | `no` | `optional` | 1136 | `accepted` |
| `--dragen-378-concordance-mode` | `` | `Boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--dragen-mode` | `` | `Boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--dragstr-het-hom-ratio` | `` | `int` | `2` | `no` | `optional` | — | `ignored_via_catalog` |
| `--dragstr-params-path` | `` | `GATKPath` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--emit-ref-confidence` | `-ERC` | `ReferenceConfidenceMode` | `NONE` | `no` | `advanced` | 1637 | `accepted` |
| `--enable-all-annotations` | `` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--enable-dynamic-read-disqualification-for-genotyping` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | 1026 | `accepted` |
| `--expected-mismatch-rate-for-read-disqualification` | `` | `double` | `0.02` | `no` | `advanced` | 1403 | `accepted` |
| `--floor-blocks` | `` | `boolean` | `false` | `no` | `advanced` | 1660 | `accepted` |
| `--flow-assembly-collapse-partial-mode` | `` | `boolean` | `false` | `no` | `advanced` | 1247 | `accepted` |
| `--flow-disallow-probs-larger-than-call` | `` | `boolean` | `false` | `no` | `advanced` | 1341 | `accepted` |
| `--flow-fill-empty-bins-value` | `` | `double` | `0.001` | `no` | `advanced` | 1309 | `accepted` |
| `--flow-filter-alleles` | `` | `boolean` | `false` | `no` | `advanced` | 1333 | `accepted` |
| `--flow-filter-alleles-qual-threshold` | `` | `float` | `30.0` | `no` | `advanced` | 1319 | `accepted` |
| `--flow-filter-alleles-sor-threshold` | `` | `float` | `3.0` | `no` | `advanced` | 1326 | `accepted` |
| `--flow-filter-lone-alleles` | `` | `boolean` | `false` | `no` | `advanced` | 1337 | `accepted` |
| `--flow-lump-probs` | `` | `boolean` | `false` | `no` | `advanced` | 1301 | `accepted` |
| `--flow-matrix-mods` | `` | `String` | `null` | `no` | `advanced` | 1367 | `accepted` |
| `--flow-mode` | `` | `FlowMode` | `NONE` | `no` | `advanced` | 1251 | `accepted` |
| `--flow-order-for-annotations` | `` | `List[String]` | `[]` | `no` | `optional` | 1377 | `accepted` |
| `--flow-probability-scaling-factor` | `` | `int` | `10` | `no` | `advanced` | 1345 | `accepted` |
| `--flow-quantization-bins` | `` | `int` | `121` | `no` | `advanced` | 1356 | `accepted` |
| `--flow-remove-non-single-base-pair-indels` | `` | `boolean` | `false` | `no` | `advanced` | 1387 | `accepted` |
| `--flow-remove-one-zero-probs` | `` | `boolean` | `false` | `no` | `advanced` | 1391 | `accepted` |
| `--flow-report-insertion-or-deletion` | `` | `boolean` | `false` | `no` | `advanced` | 1395 | `accepted` |
| `--flow-retain-max-n-probs-base-format` | `` | `boolean` | `false` | `no` | `advanced` | 1399 | `accepted` |
| `--flow-symmetric-indel-probs` | `` | `boolean` | `false` | `no` | `advanced` | 1305 | `accepted` |
| `--flow-use-t0-tag` | `` | `boolean` | `false` | `no` | `advanced` | 1297 | `accepted` |
| `--force-active` | `` | `boolean` | `false` | `no` | `advanced` | 1054 | `accepted` |
| `--force-call-filtered-alleles` | `-genotype-filtered-alleles` | `boolean` | `false` | `no` | `advanced` | 997 | `accepted` |
| `--founder-id` | `-founder-id` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--genotype-assignment-method` | `-gam` | `GenotypeAssignmentMethod` | `USE_PLS_TO_ASSIGN` | `no` | `optional` | 1600 | `accepted` |
| `--graph-output` | `-graph` | `String` | `null` | `no` | `optional` | 1009 | `accepted` |
| `--gvcf-gq-bands` | `-GQB` | `List[Integer]` | `[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 70, 80, 90, 99]` | `no` | `advanced` | 1619 | `accepted` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 860 | `accepted` |
| `--heterozygosity` | `` | `Double` | `0.001` | `no` | `optional` | 1583 | `accepted` |
| `--heterozygosity-stdev` | `` | `double` | `0.01` | `no` | `optional` | 1587 | `accepted` |
| `--indel-heterozygosity` | `` | `double` | `1.25E-4` | `no` | `optional` | 1585 | `accepted` |
| `--indel-size-to-eliminate-in-ref-model` | `` | `int` | `10` | `no` | `advanced` | 1614 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `yes` | `required` | 986 | `accepted` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | 1044 | `accepted` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | 1039 | `accepted` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | 1029 | `accepted` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 1019 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--keep-boundary-flows` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--kmer-size` | `` | `List[Integer]` | `[10, 25]` | `no` | `advanced` | 1452 | `accepted` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--likelihood-calculation-engine` | `` | `Implementation` | `PairHMM` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--linked-de-bruijn-graph` | `` | `boolean` | `false` | `no` | `advanced` | 1490 | `accepted` |
| `--mapping-quality-threshold-for-genotyping` | `` | `int` | `20` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--max-alternate-alleles` | `` | `int` | `6` | `no` | `advanced` | 1424 | `accepted` |
| `--max-assembly-region-size` | `` | `int` | `300` | `no` | `optional` | 1106 | `accepted` |
| `--max-effective-depth-adjustment-for-frd` | `` | `int` | `0` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--max-genotype-count` | `` | `int` | `1024` | `no` | `advanced` | 1433 | `accepted` |
| `--max-mnp-distance` | `-mnp-dist` | `int` | `0` | `no` | `advanced` | 1414 | `accepted` |
| `--max-num-haplotypes-in-population` | `` | `int` | `128` | `no` | `advanced` | 1554 | `accepted` |
| `--max-prob-propagation-distance` | `` | `int` | `50` | `no` | `advanced` | 1109 | `accepted` |
| `--max-reads-per-alignment-start` | `` | `int` | `50` | `no` | `optional` | 1169 | `accepted` |
| `--max-unpruned-variants` | `` | `int` | `100` | `no` | `advanced` | 1505 | `accepted` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--min-assembly-region-size` | `` | `int` | `50` | `no` | `optional` | 1103 | `accepted` |
| `--min-base-quality-score` | `-mbq` | `byte` | `10` | `no` | `optional` | 1179 | `accepted` |
| `--min-dangling-branch-length` | `` | `int` | `4` | `no` | `advanced` | 1511 | `accepted` |
| `--min-pruning` | `` | `int` | `2` | `no` | `advanced` | 1474 | `accepted` |
| `--native-pair-hmm-threads` | `` | `int` | `4` | `no` | `optional` | 1077 | `accepted` |
| `--native-pair-hmm-use-double-precision` | `` | `boolean` | `false` | `no` | `optional` | 1711 | `accepted` |
| `--num-pruning-samples` | `` | `int` | `1` | `no` | `advanced` | 1480 | `accepted` |
| `--num-reference-samples-if-no-call` | `` | `int` | `0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--output` | `-O` | `GATKPath` | `null` | `yes` | `required` | 1049 | `accepted` |
| `--output-mode` | `` | `OutputMode` | `EMIT_VARIANTS_ONLY` | `no` | `optional` | — | `ignored_via_catalog` |
| `--pair-hmm-gap-continuation-penalty` | `` | `int` | `10` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--pair-hmm-implementation` | `-pairHMM` | `Implementation` | `FASTEST_AVAILABLE` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--pair-hmm-results-file` | `` | `GATKPath` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--pcr-indel-model` | `` | `PCRErrorModel` | `CONSERVATIVE` | `no` | `advanced` | 1219 | `accepted` |
| `--pedigree` | `-ped` | `GATKPath` | `null` | `no` | `optional` | — | `unsupported` |
| `--phred-scaled-global-read-mismapping-rate` | `` | `int` | `45` | `no` | `advanced` | 1211 | `accepted` |
| `--pileup-detection` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--ploidy-regions` | `-ploidy-regions` | `FeatureInput[NamedFeature]` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--population-callset` | `-population` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--pruning-lod-threshold` | `` | `double` | `2.302585092994046` | `no` | `advanced` | 1499 | `accepted` |
| `--pruning-seeding-lod-threshold` | `` | `double` | `9.210340371976184` | `no` | `advanced` | 1502 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 1726 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | 1666 | `accepted` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--recover-all-dangling-branches` | `` | `boolean` | `false` | `no` | `advanced` | 1527 | `accepted` |
| `--recover-dangling-heads` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 991 | `accepted` |
| `--reference-model-deletion-quality` | `` | `byte` | `30` | `no` | `advanced` | 1189 | `accepted` |
| `--sample-name` | `-ALIAS` | `String` | `null` | `no` | `optional` | 1017 | `accepted` |
| `--sample-ploidy` | `-ploidy` | `int` | `2` | `no` | `optional` | 1576 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | 1729 | `accepted` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | 1704 | `accepted` |
| `--smith-waterman` | `` | `Implementation` | `FASTEST_AVAILABLE` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-dangling-end-gap-extend-penalty` | `` | `int` | `-6` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-dangling-end-gap-open-penalty` | `` | `int` | `-110` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-dangling-end-match-value` | `` | `int` | `25` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-dangling-end-mismatch-penalty` | `` | `int` | `-50` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-haplotype-to-reference-gap-extend-penalty` | `` | `int` | `-11` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-haplotype-to-reference-gap-open-penalty` | `` | `int` | `-260` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-haplotype-to-reference-match-value` | `` | `int` | `200` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-haplotype-to-reference-mismatch-penalty` | `` | `int` | `-150` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-read-to-haplotype-gap-extend-penalty` | `` | `int` | `-5` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-read-to-haplotype-gap-open-penalty` | `` | `int` | `-30` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-read-to-haplotype-match-value` | `` | `int` | `10` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-read-to-haplotype-mismatch-penalty` | `` | `int` | `-15` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--soft-clip-low-quality-ends` | `` | `boolean` | `false` | `no` | `advanced` | 1145 | `accepted` |
| `--standard-min-confidence-threshold-for-calling` | `-stand-call-conf` | `double` | `30.0` | `no` | `optional` | 1590 | `accepted` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | 1731 | `accepted` |
| `--transform-dragen-mapping-quality` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--use-filtered-reads-for-annotations` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-new-qual-calculator` | `-new-qual` | `boolean` | `true` | `no` | `optional` | — | `ignored_via_catalog` |
| `--use-pdhmm` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--use-pdhmm-overlap-optimization` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--use-posteriors-to-calculate-qual` | `-gp-qual` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | 1730 | `accepted` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--genotype-filtered-alleles` | 999 |
| `--annotate-allele-specific` | 1011 |
| `--region` | 1019 |
| `--interval` | 1020 |
| `--output-manifest` | 1061 |
| `--manifest` | 1061 |
| `--telemetry` | 1064 |
| `--batch-records` | 1066 |
| `--stream-by-contig` | 1068 |
| `--stream-by-region` | 1070 |
| `--threads` | 1077 |
| `--min-depth` | 1083 |
| `--min-alt-support` | 1087 |
| `--minimum-mapping-quality` | 1112 |
| `--min-read-length` | 1118 |
| `--max-read-length` | 1126 |
| `--include-duplicates` | 1134 |
| `--use-haplotype-realignment-for-rcm` | 1154 |
| `--do-not-correct-overlapping-base-qualities` | 1160 |
| `--max-reads-per-locus` | 1168 |
| `--downsampling-seed` | 1175 |
| `--min-base-quality` | 1178 |
| `--flow-assembly-collapse-hmer-size` | 1236 |
| `--flow-probability-threshold` | 1259 |
| `--flow-ligation` | 1269 |
| `--flow-quality` | 1279 |
| `--flow-disallow-soft-clipped` | 1289 |
| `--flow-fill-from-read-orientations` | 1293 |
| `--max-candidates` | 1412 |
| `--mnp-dist` | 1414 |
| `--min-kmer-count` | 1468 |
| `--disable-adaptive-pruning` | 1488 |
| `--disable-artificial-haplotype-recovery` | 1492 |
| `--enable-legacy-graph-cycle-detection` | 1494 |
| `--min-dangling-matching-bases` | 1518 |
| `--error-correct-reads` | 1532 |
| `--error-correction-log-odds` | 1536 |
| `--kmer-length-for-read-error-correction` | 1539 |
| `--min-observations-for-kmer-to-be-solid` | 1546 |
| `--max-haplotype-paths` | 1553 |
| `--max-haplotype-depth` | 1561 |
| `--max-haplotype-combination-alleles` | 1566 |
| `--haplotype-pruning-log10` | 1573 |
| `--stand-call-conf` | 1591 |
| `--no-genotype-priors` | 1617 |
| `--gvcf` | 1637 |
| `--java-options` | 1730 |

## `fastgatk-bqsr` ↔ `BaseRecalibrator`

- Source: `fastgatk-native/src/bqsr_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--binary-tag-name` | `` | `String` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--bqsr-baq-gap-open-penalty` | `` | `double` | `40.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | 408 | `accepted` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--default-base-qualities` | `` | `byte` | `-1` | `no` | `optional` | — | `ignored_via_catalog` |
| `--deletions-default-quality` | `` | `byte` | `45` | `no` | `optional` | 397 | `accepted` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | 373 | `accepted` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | 362 | `accepted` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 250 | `accepted` |
| `--indels-context-size` | `-ics` | `int` | `3` | `no` | `optional` | 390 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `yes` | `required` | 293 | `accepted` |
| `--insertions-default-quality` | `` | `byte` | `45` | `no` | `optional` | 394 | `accepted` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | — | `ignored_via_catalog` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 319 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | 379 | `accepted` |
| `--known-sites` | `` | `List[FeatureInput[Feature]]` | `[]` | `yes` | `required` | 297 | `accepted` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--low-quality-tail` | `` | `byte` | `2` | `no` | `optional` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--maximum-cycle-value` | `-max-cycle` | `int` | `500` | `no` | `optional` | 385 | `accepted` |
| `--mismatches-context-size` | `-mcs` | `int` | `2` | `no` | `optional` | 358 | `accepted` |
| `--mismatches-default-quality` | `` | `byte` | `-1` | `no` | `optional` | — | `ignored_via_catalog` |
| `--output` | `-O` | `GATKPath` | `null` | `yes` | `required` | 314 | `accepted` |
| `--preserve-qscores-less-than` | `` | `int` | `6` | `no` | `optional` | 325 | `accepted` |
| `--quantizing-levels` | `` | `int` | `16` | `no` | `optional` | — | `ignored_via_catalog` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 412 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | 368 | `accepted` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 295 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-original-qualities` | `-OQ` | `Boolean` | `false` | `no` | `optional` | 345 | `accepted` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `-bqsr` | 301 |
| `--bqsr-recal-file` | 303 |
| `--recal-file` | 304 |
| `--output-manifest` | 316 |
| `--manifest` | 316 |
| `--region` | 319 |
| `--batch-records` | 323 |
| `--quantize-quals` | 328 |
| `--static-quantized-quals` | 331 |
| `--round-down-quantized` | 334 |
| `--allow-missing-read-group` | 338 |
| `--global-qscore-prior` | 342 |
| `--emit-original-quals` | 351 |
| `--compute-indel-bqsr-tables` | 355 |
| `--indels` | 355 |
| `--max-cycle` | 385 |
| `--checkpoint` | 400 |
| `--resume-checkpoint` | 402 |
| `--checkpoint-every-batches` | 405 |
| `mismatches_context_size` | 1595 |

## `fastgatk-apply-bqsr` ↔ `ApplyBQSR`

- Source: `fastgatk-native/src/bqsr_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--allow-missing-read-group` | `` | `boolean` | `false` | `no` | `optional` | 338 | `accepted` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--bqsr-recal-file` | `-bqsr` | `File` | `null` | `yes` | `required` | 303 | `accepted` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | 408 | `accepted` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | 373 | `accepted` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | 362 | `accepted` |
| `--emit-original-quals` | `` | `boolean` | `false` | `no` | `optional` | 351 | `accepted` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--global-qscore-prior` | `` | `double` | `-1.0` | `no` | `optional` | 342 | `accepted` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 250 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `yes` | `required` | 293 | `accepted` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `unsupported` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | — | `unsupported` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 319 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | 379 | `accepted` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--output` | `-O` | `GATKPath` | `null` | `yes` | `required` | 314 | `accepted` |
| `--preserve-qscores-less-than` | `` | `int` | `6` | `no` | `optional` | 325 | `accepted` |
| `--quantize-quals` | `` | `int` | `0` | `no` | `optional` | 328 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 412 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | 368 | `accepted` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `no` | `optional` | 295 | `accepted` |
| `--round-down-quantized` | `` | `boolean` | `false` | `no` | `advanced` | 334 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--static-quantized-quals` | `` | `List[Integer]` | `[]` | `no` | `advanced` | 331 | `accepted` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-original-qualities` | `-OQ` | `Boolean` | `false` | `no` | `optional` | 345 | `accepted` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--known-sites` | 297 |
| `--recal-file` | 304 |
| `--output-manifest` | 316 |
| `--manifest` | 316 |
| `--region` | 319 |
| `--batch-records` | 323 |
| `--compute-indel-bqsr-tables` | 355 |
| `--indels` | 355 |
| `--mismatches-context-size` | 358 |
| `-mcs` | 358 |
| `--maximum-cycle-value` | 385 |
| `--max-cycle` | 385 |
| `--indels-context-size` | 390 |
| `-ics` | 390 |
| `--insertions-default-quality` | 394 |
| `--deletions-default-quality` | 397 |
| `--checkpoint` | 400 |
| `--resume-checkpoint` | 402 |
| `--checkpoint-every-batches` | 405 |
| `mismatches_context_size` | 1595 |

## `fastgatk-genotype-gvcf` ↔ `GenotypeGVCFs`

- Source: `fastgatk-native/src/genotype_gvcf_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_GenotypeGVCFs.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--allele-fraction-error` | `` | `double` | `0.001` | `no` | `optional` | — | `ignored_via_catalog` |
| `--annotate-with-num-discovered-alleles` | `` | `boolean` | `false` | `no` | `optional` | 663 | `accepted` |
| `--annotation` | `-A` | `List[String]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--annotation-group` | `-G` | `List[String]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--annotations-to-exclude` | `-AX` | `List[String]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--call-genotypes` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | 628 | `accepted` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--dbsnp` | `-D` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | 647 | `accepted` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | 706 | `accepted` |
| `--disable-tool-default-annotations` | `-disable-tool-default-annotations` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--dont-use-dragstr-priors` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--enable-all-annotations` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | 610 | `accepted` |
| `--flow-order-for-annotations` | `` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--force-output-intervals` | `` | `List[String]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--founder-id` | `-founder-id` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--genomicsdb-max-alternate-alleles` | `` | `int` | `50` | `no` | `optional` | — | `ignored_via_catalog` |
| `--genomicsdb-shared-posixfs-optimizations` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--genomicsdb-use-bcf-codec` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--genomicsdb-use-gcs-hdfs-connector` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--genotype-assignment-method` | `-gam` | `GenotypeAssignmentMethod` | `USE_PLS_TO_ASSIGN` | `no` | `optional` | 698 | `accepted` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 564 | `accepted` |
| `--heterozygosity` | `` | `Double` | `0.001` | `no` | `optional` | 632 | `accepted` |
| `--heterozygosity-stdev` | `` | `double` | `0.01` | `no` | `optional` | 638 | `accepted` |
| `--include-non-variant-sites` | `-all-sites` | `boolean` | `false` | `no` | `optional` | 669 | `accepted` |
| `--indel-heterozygosity` | `` | `double` | `1.25E-4` | `no` | `optional` | 635 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--input-is-somatic` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | 613 | `accepted` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 603 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--keep-combined-raw-annotations` | `-keep-combined` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--keep-specific-combined-raw-annotation` | `-keep-specific-combined` | `List[String]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-alternate-alleles` | `` | `int` | `6` | `no` | `advanced` | 660 | `accepted` |
| `--max-genotype-count` | `` | `int` | `1024` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--merge-input-intervals` | `-merge-input-intervals` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--num-reference-samples-if-no-call` | `` | `int` | `0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--only-output-calls-starting-in-intervals` | `` | `boolean` | `false` | `no` | `advanced` | 671 | `accepted` |
| `--output` | `-O` | `GATKPath` | `null` | `yes` | `required` | 623 | `accepted` |
| `--pedigree` | `-ped` | `GATKPath` | `null` | `no` | `optional` | — | `unsupported` |
| `--population-callset` | `-population` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 706 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 601 | `accepted` |
| `--sample-ploidy` | `-ploidy` | `int` | `2` | `no` | `optional` | 649 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--standard-min-confidence-threshold-for-calling` | `-stand-call-conf` | `double` | `30.0` | `no` | `optional` | 641 | `accepted` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--tumor-lod-to-emit` | `-emit-lod` | `double` | `3.5` | `no` | `optional` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-new-qual-calculator` | `-new-qual` | `boolean` | `true` | `no` | `optional` | 684 | `accepted` |
| `--use-posteriors-to-calculate-qual` | `-gp-qual` | `boolean` | `false` | `no` | `optional` | 691 | `accepted` |
| `--variant` | `-V` | `String` | `null` | `yes` | `required` | 599 | `accepted` |
| `--variant-output-filtering` | `` | `Mode` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | 709 | `accepted` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--interval` | 603 |
| `--region` | 604 |
| `--output-manifest` | 625 |
| `--manifest` | 625 |
| `--stand-call-conf` | 642 |
| `--no-genotype-priors` | 667 |
| `--gatk-compatible-annotations` | 677 |
| `--strict-gatk-annotations` | 678 |
| `--stream-by-locus` | 680 |
| `--stream-loci` | 680 |
| `--no-use-new-qual-calculator` | 682 |
| `--new-qual` | 686 |
| `--gp-qual` | 693 |
| `--gam` | 699 |
| `--java-options` | 709 |

## `fastgatk-reblock-gvcf` ↔ `ReblockGVCF`

- Source: `fastgatk-native/src/reblock_gvcf_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_variantutils_ReblockGVCF.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-site-filters-to-genotype` | `` | `boolean` | `false` | `no` | `advanced` | 391 | `accepted` |
| `--allow-missing-hom-ref-data` | `` | `boolean` | `false` | `no` | `advanced` | 378 | `accepted` |
| `--annotate-with-num-discovered-alleles` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--annotation` | `-A` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--annotation-group` | `-G` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--annotations-to-exclude` | `-AX` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--annotations-to-keep` | `` | `List[String]` | `[]` | `no` | `advanced` | 379 | `accepted` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | 392 | `accepted` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--dbsnp` | `-D` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | 396 | `accepted` |
| `--disable-tool-default-annotations` | `-disable-tool-default-annotations` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--do-qual-score-approximation` | `-do-qual-approx` | `boolean` | `false` | `no` | `advanced` | 375 | `accepted` |
| `--dont-use-dragstr-priors` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--drop-low-quals` | `-drop-low-quals` | `boolean` | `false` | `no` | `advanced` | 373 | `accepted` |
| `--enable-all-annotations` | `` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--floor-blocks` | `` | `boolean` | `false` | `no` | `advanced` | 372 | `accepted` |
| `--flow-order-for-annotations` | `` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--format-annotations-to-remove` | `` | `List[String]` | `[]` | `no` | `advanced` | 382 | `accepted` |
| `--founder-id` | `-founder-id` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--genotype-assignment-method` | `-gam` | `GenotypeAssignmentMethod` | `USE_PLS_TO_ASSIGN` | `no` | `optional` | — | `ignored_via_catalog` |
| `--gvcf-gq-bands` | `-GQB` | `List[Integer]` | `[20, 100]` | `no` | `advanced` | 351 | `accepted` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 311 | `accepted` |
| `--heterozygosity` | `` | `Double` | `0.001` | `no` | `optional` | — | `ignored_via_catalog` |
| `--heterozygosity-stdev` | `` | `double` | `0.01` | `no` | `optional` | — | `ignored_via_catalog` |
| `--indel-heterozygosity` | `` | `double` | `1.25E-4` | `no` | `optional` | — | `ignored_via_catalog` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | 337 | `accepted` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 331 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--keep-all-alts` | `` | `boolean` | `false` | `no` | `advanced` | 374 | `accepted` |
| `--keep-site-filters` | `-keep-filters` | `boolean` | `false` | `no` | `advanced` | 390 | `accepted` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-alternate-alleles` | `` | `int` | `6` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--max-genotype-count` | `` | `int` | `1024` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--num-reference-samples-if-no-call` | `` | `int` | `0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--output` | `-O` | `GATKPath` | `null` | `yes` | `required` | 346 | `accepted` |
| `--pedigree` | `-ped` | `GATKPath` | `null` | `no` | `optional` | — | `unsupported` |
| `--population-callset` | `-population` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 396 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 329 | `accepted` |
| `--rgq-threshold-to-no-call` | `-rgq-threshold` | `double` | `0.0` | `no` | `advanced` | 359 | `accepted` |
| `--sample-ploidy` | `-ploidy` | `int` | `2` | `no` | `optional` | — | `ignored_via_catalog` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--standard-min-confidence-threshold-for-calling` | `-stand-call-conf` | `double` | `30.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--tree-score-threshold-to-no-call` | `` | `double` | `0.0` | `no` | `advanced` | 369 | `accepted` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-new-qual-calculator` | `-new-qual` | `boolean` | `true` | `no` | `optional` | — | `ignored_via_catalog` |
| `--use-posteriors-to-calculate-qual` | `-gp-qual` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--variant` | `-V` | `List[GATKPath]` | `[]` | `yes` | `required` | 327 | `accepted` |
| `--variant-output-filtering` | `` | `Mode` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | 398 | `accepted` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--interval` | 332 |
| `--region` | 332 |
| `--output-manifest` | 348 |
| `--manifest` | 348 |
| `--rgq-threshold` | 359 |
| `--do-qual-score-approx` | 376 |
| `--do-qual-approx` | 376 |
| `--annotations-to-remove` | 383 |
| `--keep-filters` | 390 |
| `--java-options` | 398 |

## `fastgatk-select-variants` ↔ `SelectVariants`

- Source: `fastgatk-native/src/select_variants_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_variantutils_SelectVariants.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--apply-jexl-filters-first` | `-jexl-first` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--call-genotypes` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--concordance` | `-conc` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | 211 | `accepted` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | 249 | `accepted` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | 260 | `accepted` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--discordance` | `-disc` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | 213 | `accepted` |
| `--drop-genotype-annotation` | `-DGA` | `List[String]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--drop-info-annotation` | `-DA` | `List[String]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--exclude-filtered` | `` | `boolean` | `false` | `no` | `optional` | 218 | `accepted` |
| `--exclude-ids` | `-xl-ids` | `Set[String]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--exclude-non-variants` | `` | `boolean` | `false` | `no` | `optional` | 226 | `accepted` |
| `--exclude-sample-expressions` | `-xl-se` | `Set[String]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--exclude-sample-name` | `-xl-sn` | `Set[String]` | `[]` | `no` | `optional` | 199 | `accepted` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--genomicsdb-max-alternate-alleles` | `` | `int` | `50` | `no` | `optional` | — | `ignored_via_catalog` |
| `--genomicsdb-shared-posixfs-optimizations` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--genomicsdb-use-bcf-codec` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--genomicsdb-use-gcs-hdfs-connector` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 125 | `accepted` |
| `--ignore-non-ref-in-types` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | 165 | `accepted` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 158 | `accepted` |
| `--invert-mendelian-violation` | `` | `Boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--invertSelect` | `-invert-select` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--keep-ids` | `-ids` | `Set[String]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--keep-original-ac` | `` | `boolean` | `false` | `no` | `optional` | 239 | `accepted` |
| `--keep-original-dp` | `` | `boolean` | `false` | `no` | `optional` | 243 | `accepted` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-filtered-genotypes` | `` | `int` | `2147483647` | `no` | `optional` | — | `ignored_via_catalog` |
| `--max-fraction-filtered-genotypes` | `` | `double` | `1.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--max-indel-size` | `` | `int` | `2147483647` | `no` | `optional` | 189 | `accepted` |
| `--max-nocall-fraction` | `` | `double` | `1.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--max-nocall-number` | `` | `int` | `2147483647` | `no` | `optional` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--mendelian-violation` | `` | `Boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--mendelian-violation-qual-threshold` | `` | `double` | `0.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--min-filtered-genotypes` | `` | `int` | `0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--min-fraction-filtered-genotypes` | `` | `double` | `0.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--min-indel-size` | `` | `int` | `0` | `no` | `optional` | 183 | `accepted` |
| `--output` | `-O` | `GATKPath` | `null` | `yes` | `required` | 154 | `accepted` |
| `--pedigree` | `-ped` | `GATKPath` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--preserve-alleles` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 260 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `no` | `optional` | 156 | `accepted` |
| `--remove-fraction-genotypes` | `` | `double` | `0.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--remove-unused-alternates` | `` | `boolean` | `false` | `no` | `optional` | 248 | `accepted` |
| `--restrict-alleles-to` | `` | `NumberAlleleRestriction` | `ALL` | `no` | `optional` | 181 | `accepted` |
| `--sample-expressions` | `-se` | `Set[String]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--sample-name` | `-sn` | `Set[String]` | `[]` | `no` | `optional` | 196 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--select` | `` | `ArrayList[String]` | `[]` | `no` | `optional` | 202 | `accepted` |
| `--select-genotype-expressions` | `-select-genotype` | `ArrayList[String]` | `[]` | `no` | `optional` | 206 | `accepted` |
| `--select-random-fraction` | `-fraction` | `double` | `1.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--select-type-to-exclude` | `-xl-select-type` | `List[Type]` | `[]` | `no` | `optional` | 179 | `accepted` |
| `--select-type-to-include` | `-select-type` | `List[Type]` | `[]` | `no` | `optional` | 177 | `accepted` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--set-filtered-gt-to-nocall` | `` | `boolean` | `false` | `no` | `optional` | 234 | `accepted` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | 253 | `accepted` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--variant` | `-V` | `GATKPath` | `null` | `yes` | `required` | 152 | `accepted` |
| `--variant-output-filtering` | `` | `Mode` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | 262 | `accepted` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--interval` | 159 |
| `--region` | 159 |
| `--output-manifest` | 174 |
| `--manifest` | 174 |
| `--sample-name-to-include` | 196 |
| `--sample-name-to-exclude` | 199 |
| `--select-expression` | 202 |
| `--select-genotype` | 206 |
| `--concordance-genotypes` | 215 |
| `--discordance-genotypes` | 215 |
| `--sample-level-concordance` | 216 |
| `--exclude-filtered-variants` | 222 |
| `--exclude-non-variant-sites` | 230 |
| `--java-options` | 262 |

## `fastgatk-gather-vcfs` ↔ `GatherVcfs (Picard)`

- Source: `fastgatk-native/src/gather_vcfs_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/picard_vcf_GatherVcfs.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--COMMENT` | `-CO` | `List[String]` | `[]` | `no` | `optional` | 143 | `accepted` |
| `--COMPRESSION_LEVEL` | `` | `int` | `5` | `no` | `common` | 170 | `accepted` |
| `--CREATE_INDEX` | `` | `Boolean` | `true` | `no` | `common` | 165 | `accepted` |
| `--CREATE_MD5_FILE` | `` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 102 | `accepted` |
| `--INPUT` | `-I` | `List[HtsPath]` | `[]` | `yes` | `required` | 117 | `accepted` |
| `--MAX_RECORDS_IN_RAM` | `` | `Integer` | `500000` | `no` | `common` | — | `unsupported` |
| `--OUTPUT` | `-O` | `File` | `null` | `yes` | `required` | 145 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 179 | `accepted` |
| `--REFERENCE_SEQUENCE` | `-R` | `PicardHtsPath` | `null` | `no` | `common` | 121 | `accepted` |
| `--REORDER_INPUT_BY_FIRST_VARIANT` | `-RI` | `boolean` | `false` | `no` | `optional` | 154 | `accepted` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--TMP_DIR` | `` | `List[File]` | `[]` | `no` | `common` | — | `unsupported` |
| `--USE_JDK_DEFLATER` | `-use_jdk_deflater` | `Boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--USE_JDK_INFLATER` | `-use_jdk_inflater` | `Boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--VALIDATION_STRINGENCY` | `` | `ValidationStringency` | `STRICT` | `no` | `common` | — | `unsupported` |
| `--VERBOSITY` | `` | `LogLevel` | `INFO` | `no` | `common` | 186 | `accepted` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--reference` | 121 |
| `-L` | 127 |
| `--intervals` | 127 |
| `--interval` | 128 |
| `--region` | 128 |
| `-isr` | 134 |
| `--interval-set-rule` | 134 |
| `--output-manifest` | 149 |
| `--manifest` | 149 |
| `--allow-overlaps` | 152 |
| `--reorder-input-by-first-variant` | 154 |
| `--create-output-variant-index` | 163 |
| `--disable-sequence-dictionary-validation` | 183 |
| `--java-options` | 185 |

## `fastgatk-left-align-trim` ↔ `LeftAlignAndTrimVariants`

- Source: `fastgatk-native/src/left_align_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_variantutils_LeftAlignAndTrimVariants.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | 207 | `accepted` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | 224 | `accepted` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--dont-trim-alleles` | `-no-trim` | `boolean` | `false` | `no` | `optional` | 188 | `accepted` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | 147 | `accepted` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 118 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | 160 | `accepted` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | 175 | `accepted` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | 154 | `accepted` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | 166 | `accepted` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 140 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--keep-original-ac` | `` | `boolean` | `false` | `no` | `optional` | 216 | `accepted` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-indel-length` | `` | `int` | `200` | `no` | `optional` | 194 | `accepted` |
| `--max-leading-bases` | `` | `int` | `1000` | `no` | `optional` | 199 | `accepted` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--output` | `-O` | `GATKPath` | `null` | `yes` | `required` | 183 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 224 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 138 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | 220 | `accepted` |
| `--split-multi-allelics` | `` | `boolean` | `false` | `no` | `optional` | 212 | `accepted` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--variant` | `-V` | `GATKPath` | `null` | `yes` | `required` | 136 | `accepted` |
| `--variant-output-filtering` | `` | `Mode` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | 226 | `accepted` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--interval` | 141 |
| `--region` | 141 |
| `--exclude-interval` | 148 |
| `--exclude-region` | 148 |
| `--output-manifest` | 185 |
| `--manifest` | 185 |
| `--no-trim` | 188 |
| `--java-options` | 226 |

## `fastgatk-funcotator` ↔ `Funcotator`

- Source: `fastgatk-native/src/funcotator_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_funcotator_Funcotator.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `unsupported` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `unsupported` |
| `--annotation-default` | `` | `List[String]` | `[]` | `no` | `optional` | 106 | `accepted` |
| `--annotation-override` | `` | `List[String]` | `[]` | `no` | `optional` | 109 | `accepted` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `unsupported` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `unsupported` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--custom-variant-classification-order` | `` | `GATKPath` | `null` | `no` | `optional` | — | `unsupported` |
| `--data-sources-path` | `` | `List[String]` | `[]` | `yes` | `required` | 97 | `accepted` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--exclude-field` | `` | `Set[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--five-prime-flank-size` | `` | `int` | `5000` | `no` | `optional` | — | `unsupported` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 74 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `unsupported` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `unsupported` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `unsupported` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | — | `unsupported` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--lookahead-cache-bp` | `` | `int` | `100000` | `no` | `optional` | — | `unsupported` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `unsupported` |
| `--min-num-bases-for-segment-funcotation` | `` | `int` | `150` | `no` | `advanced` | — | `unsupported` |
| `--output` | `-O` | `File` | `null` | `yes` | `required` | 89 | `accepted` |
| `--output-file-format` | `` | `OutputFormatType` | `null` | `yes` | `required` | 100 | `accepted` |
| `--prefer-mane-transcripts` | `` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `unsupported` |
| `--reannotate-vcf` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--ref-version` | `` | `String` | `null` | `yes` | `required` | 95 | `accepted` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 92 | `accepted` |
| `--remove-filtered-variants` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `unsupported` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `unsupported` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--splice-site-window-size` | `` | `int` | `2` | `no` | `optional` | — | `unsupported` |
| `--three-prime-flank-size` | `` | `int` | `0` | `no` | `optional` | — | `unsupported` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `unsupported` |
| `--transcript-list` | `` | `Set[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--transcript-selection-mode` | `` | `TranscriptSelectionMode` | `CANONICAL` | `no` | `optional` | 103 | `accepted` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--variant` | `-V` | `GATKPath` | `null` | `yes` | `required` | 86 | `accepted` |
| `--variant-output-filtering` | `` | `Mode` | `null` | `no` | `advanced` | — | `unsupported` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `unsupported` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

## `fastgatk-variant-annotator` ↔ `VariantAnnotator`

- Source: `fastgatk-native/src/variant_annotator_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_annotator_VariantAnnotator.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `unsupported` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `unsupported` |
| `--annotation` | `-A` | `List[String]` | `[]` | `no` | `optional` | 121 | `accepted` |
| `--annotation-group` | `-G` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--annotations-to-exclude` | `-AX` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--comparison` | `-comp` | `List[FeatureInput[VariantContext]]` | `[]` | `no` | `advanced` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `unsupported` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | 172 | `accepted` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--dbsnp` | `-D` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | 162 | `accepted` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--disable-tool-default-annotations` | `-disable-tool-default-annotations` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--enable-all-annotations` | `` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--expression` | `-E` | `Set[String]` | `[]` | `no` | `optional` | 124 | `accepted` |
| `--flow-order-for-annotations` | `` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--founder-id` | `-founder-id` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 92 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | 118 | `accepted` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `unsupported` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `unsupported` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `unsupported` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | — | `unsupported` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 164 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `unsupported` |
| `--min-base-quality-score` | `` | `byte` | `10` | `no` | `optional` | 169 | `accepted` |
| `--output` | `-O` | `File` | `null` | `yes` | `required` | 112 | `accepted` |
| `--pedigree` | `-ped` | `GATKPath` | `null` | `no` | `optional` | — | `unsupported` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `unsupported` |
| `--reference` | `-R` | `GATKPath` | `null` | `no` | `optional` | 115 | `accepted` |
| `--resource` | `` | `List[FeatureInput[VariantContext]]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--resource-allele-concordance` | `-rac` | `Boolean` | `false` | `no` | `optional` | 167 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `unsupported` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `unsupported` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `unsupported` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--variant` | `-V` | `GATKPath` | `null` | `yes` | `required` | 109 | `accepted` |
| `--variant-output-filtering` | `` | `Mode` | `null` | `no` | `advanced` | — | `unsupported` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `unsupported` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--interval` | 164 |

## `fastgatk-variant-filtration` ↔ `VariantFiltration`

- Source: `fastgatk-native/src/variant_filtration_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_filters_VariantFiltration.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--apply-allele-specific-filters` | `` | `boolean` | `false` | `no` | `optional` | 206 | `accepted` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--cluster-size` | `-cluster` | `Integer` | `3` | `no` | `optional` | 177 | `accepted` |
| `--cluster-window-size` | `-window` | `Integer` | `0` | `no` | `optional` | 183 | `accepted` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | 220 | `accepted` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | 224 | `accepted` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--filter-expression` | `-filter` | `List[String]` | `[]` | `no` | `optional` | 155 | `accepted` |
| `--filter-name` | `` | `List[String]` | `[]` | `no` | `optional` | 157 | `accepted` |
| `--filter-not-in-mask` | `` | `boolean` | `false` | `no` | `optional` | 201 | `accepted` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--genotype-filter-expression` | `-G-filter` | `List[String]` | `[]` | `no` | `optional` | 159 | `accepted` |
| `--genotype-filter-name` | `-G-filter-name` | `List[String]` | `[]` | `no` | `optional` | 162 | `accepted` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 110 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | 146 | `accepted` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 139 | `accepted` |
| `--invalidate-previous-filters` | `` | `boolean` | `false` | `no` | `optional` | 208 | `accepted` |
| `--invert-filter-expression` | `-invfilter` | `boolean` | `false` | `no` | `optional` | 202 | `accepted` |
| `--invert-genotype-filter-expression` | `-invG-filter` | `boolean` | `false` | `no` | `optional` | 204 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--mask` | `-mask` | `FeatureInput[Feature]` | `null` | `no` | `optional` | 165 | `accepted` |
| `--mask-description` | `` | `String` | `null` | `no` | `optional` | 169 | `accepted` |
| `--mask-extension` | `` | `Integer` | `0` | `no` | `optional` | 171 | `accepted` |
| `--mask-name` | `` | `String` | `Mask` | `no` | `optional` | 167 | `accepted` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--missing-values-evaluate-as-failing` | `` | `Boolean` | `false` | `no` | `optional` | 193 | `accepted` |
| `--output` | `-O` | `GATKPath` | `null` | `yes` | `required` | 135 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 224 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `no` | `optional` | 137 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--set-filtered-genotype-to-no-call` | `` | `boolean` | `false` | `no` | `optional` | 212 | `accepted` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--variant` | `-V` | `GATKPath` | `null` | `yes` | `required` | 133 | `accepted` |
| `--variant-output-filtering` | `` | `Mode` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | 226 | `accepted` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--interval` | 140 |
| `--region` | 140 |
| `--output-manifest` | 190 |
| `--manifest` | 190 |
| `--invfilter` | 202 |
| `--invG-filter` | 204 |
| `--java-options` | 226 |

## `fastgatk-apply-vqsr` ↔ `ApplyVQSR`

- Source: `fastgatk-native/src/apply_vqsr_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_vqsr_ApplyVQSR.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | 238 | `accepted` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | 258 | `accepted` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--exclude-filtered` | `` | `boolean` | `false` | `no` | `optional` | 254 | `accepted` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | 208 | `accepted` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 163 | `accepted` |
| `--ignore-all-filters` | `` | `boolean` | `false` | `no` | `optional` | 246 | `accepted` |
| `--ignore-filter` | `` | `List[String]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | 214 | `accepted` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | — | `ignored_via_catalog` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 202 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--lod-score-cutoff` | `` | `Double` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--mode` | `-mode` | `Mode` | `SNP` | `no` | `optional` | — | `ignored_via_catalog` |
| `--output` | `-O` | `GATKPath` | `null` | `yes` | `required` | 200 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 258 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--recal-file` | `` | `FeatureInput[VariantContext]` | `null` | `yes` | `required` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `no` | `optional` | — | `unsupported` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | 242 | `accepted` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--tranches-file` | `` | `GATKPath` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--truth-sensitivity-filter-level` | `-ts-filter-level` | `Double` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--use-allele-specific-annotations` | `-AS` | `boolean` | `false` | `no` | `optional` | 231 | `accepted` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--variant` | `-V` | `List[GATKPath]` | `[]` | `yes` | `required` | 183 | `accepted` |
| `--variant-output-filtering` | `` | `Mode` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--AS` | 233 |

## `fastgatk-variant-recalibrator` ↔ `VariantRecalibrator`

- Source: `fastgatk-native/src/variant_recalibrator_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_vqsr_VariantRecalibrator.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--aggregate` | `-aggregate` | `List[FeatureInput[VariantContext]]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--bad-lod-score-cutoff` | `-bad-lod-cutoff` | `double` | `-5.0` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | 471 | `accepted` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--debug-stdev-thresholding` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--dirichlet` | `` | `double` | `0.001` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | 494 | `accepted` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--dont-run-rscript` | `` | `boolean` | `false` | `no` | `advanced` | 486 | `accepted` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 302 | `accepted` |
| `--ignore-all-filters` | `` | `boolean` | `false` | `no` | `optional` | 482 | `accepted` |
| `--ignore-filter` | `` | `List[String]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--input-model` | `` | `GATKPath` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | — | `ignored_via_catalog` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--k-means-iterations` | `` | `int` | `100` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-attempts` | `` | `int` | `1` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--max-gaussians` | `` | `int` | `8` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--max-iterations` | `` | `int` | `150` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--max-negative-gaussians` | `` | `int` | `2` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--maximum-training-variants` | `` | `int` | `2500000` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--minimum-bad-variants` | `` | `int` | `1000` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--mode` | `-mode` | `Mode` | `SNP` | `no` | `optional` | — | `ignored_via_catalog` |
| `--mq-cap-for-logit-jitter-transform` | `-mq-cap` | `int` | `0` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--mq-jitter` | `` | `double` | `0.05` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--output` | `-O` | `GATKPath` | `null` | `yes` | `required` | 361 | `accepted` |
| `--output-model` | `` | `GATKPath` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--prior-counts` | `` | `double` | `20.0` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 479 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `no` | `optional` | 458 | `accepted` |
| `--resource` | `` | `List[FeatureInput[VariantContext]]` | `[]` | `yes` | `required` | 345 | `accepted` |
| `--rscript-file` | `` | `File` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--shrinkage` | `` | `double` | `1.0` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | 475 | `accepted` |
| `--standard-deviation-threshold` | `-std` | `double` | `10.0` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--target-titv` | `-titv` | `double` | `2.15` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--tranches-file` | `` | `File` | `null` | `yes` | `required` | — | `ignored_via_catalog` |
| `--trust-all-polymorphic` | `` | `boolean` | `false` | `no` | `advanced` | 490 | `accepted` |
| `--truth-sensitivity-tranche` | `-tranche` | `List[Double]` | `[100.0, 99.9, 99.0, 90.0]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--use-allele-specific-annotations` | `-AS` | `boolean` | `false` | `no` | `optional` | 461 | `accepted` |
| `--use-annotation` | `-an` | `List[String]` | `[]` | `yes` | `required` | 351 | `accepted` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--variant` | `-V` | `List[GATKPath]` | `[]` | `yes` | `required` | 336 | `accepted` |
| `--variant-output-filtering` | `` | `Mode` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--full-covariance` | 451 |
| `--AS` | 460 |
| `--output-tranches-for-scatter` | 467 |

## `fastgatk-sort-sam` ↔ `SortSam (Picard)`

- Source: `fastgatk-native/src/sort_sam_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/picard_sam_SortSam.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--COMPRESSION_LEVEL` | `` | `int` | `5` | `no` | `common` | 161 | `accepted` |
| `--CREATE_INDEX` | `` | `Boolean` | `false` | `no` | `common` | 191 | `accepted` |
| `--CREATE_MD5_FILE` | `` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 105 | `accepted` |
| `--INPUT` | `-I` | `File` | `null` | `yes` | `required` | 128 | `accepted` |
| `--MAX_RECORDS_IN_RAM` | `` | `Integer` | `500000` | `no` | `common` | 152 | `accepted` |
| `--OUTPUT` | `-O` | `File` | `null` | `yes` | `required` | 132 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 202 | `accepted` |
| `--REFERENCE_SEQUENCE` | `-R` | `PicardHtsPath` | `null` | `no` | `common` | 136 | `accepted` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--SORT_ORDER` | `-SO` | `SortOrder` | `null` | `yes` | `required` | 145 | `accepted` |
| `--TMP_DIR` | `` | `List[File]` | `[]` | `no` | `common` | 157 | `accepted` |
| `--USE_JDK_DEFLATER` | `-use_jdk_deflater` | `Boolean` | `false` | `no` | `common` | 208 | `accepted` |
| `--USE_JDK_INFLATER` | `-use_jdk_inflater` | `Boolean` | `false` | `no` | `common` | 208 | `accepted` |
| `--VALIDATION_STRINGENCY` | `` | `ValidationStringency` | `STRICT` | `no` | `common` | 168 | `accepted` |
| `--VERBOSITY` | `` | `LogLevel` | `INFO` | `no` | `common` | 217 | `accepted` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--reference` | 136 |
| `--output-manifest` | 142 |
| `--manifest` | 142 |
| `--sort-order` | 145 |
| `--max-records-in-memory` | 152 |
| `--tmp-dir` | 157 |
| `--resume-spill` | 177 |
| `--add-pg-tag` | 179 |
| `--program-record-id` | 183 |
| `--program-group-name` | 185 |
| `--program-group-version` | 187 |
| `--program-group-command-line` | 189 |
| `--create-output-bam-index` | 192 |
| `--create-output-variant-index` | 194 |
| `--disable-sequence-dictionary-validation` | 206 |
| `--use-jdk-deflater` | 207 |
| `--use-jdk-inflater` | 207 |
| `--java-options` | 216 |

## `fastgatk-mark-duplicates` ↔ `MarkDuplicates (Picard)`

- Source: `fastgatk-native/src/mark_duplicates_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/picard_sam_markduplicates_MarkDuplicates.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--ADD_PG_TAG_TO_READS` | `` | `boolean` | `true` | `no` | `common` | — | `unsupported` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--ASSUME_SORT_ORDER` | `-ASO` | `SortOrder` | `null` | `no` | `optional` | — | `unsupported` |
| `--ASSUME_SORTED` | `-AS` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--BARCODE_TAG` | `` | `String` | `null` | `no` | `optional` | — | `unsupported` |
| `--CLEAR_DT` | `` | `boolean` | `true` | `no` | `optional` | — | `unsupported` |
| `--COMMENT` | `-CO` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--COMPRESSION_LEVEL` | `` | `int` | `5` | `no` | `common` | — | `unsupported` |
| `--CREATE_INDEX` | `` | `Boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--CREATE_MD5_FILE` | `` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--DUPLEX_UMI` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--DUPLICATE_SCORING_STRATEGY` | `-DS` | `ScoringStrategy` | `SUM_OF_BASE_QUALITIES` | `no` | `optional` | — | `unsupported` |
| `--FLOW_DUP_STRATEGY` | `` | `FLOW_DUPLICATE_SELECTION_STRATEGY` | `FLOW_QUALITY_SUM_STRATEGY` | `no` | `optional` | — | `unsupported` |
| `--FLOW_EFFECTIVE_QUALITY_THRESHOLD` | `` | `int` | `15` | `no` | `optional` | — | `unsupported` |
| `--FLOW_MODE` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--FLOW_Q_IS_KNOWN_END` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--FLOW_SKIP_FIRST_N_FLOWS` | `` | `int` | `0` | `no` | `optional` | — | `unsupported` |
| `--FLOW_UNPAIRED_END_UNCERTAINTY` | `` | `int` | `0` | `no` | `optional` | — | `unsupported` |
| `--FLOW_UNPAIRED_START_UNCERTAINTY` | `` | `int` | `0` | `no` | `optional` | — | `unsupported` |
| `--FLOW_USE_END_IN_UNPAIRED_READS` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--FLOW_USE_UNPAIRED_CLIPPED_END` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 113 | `accepted` |
| `--INPUT` | `-I` | `List[String]` | `[]` | `yes` | `required` | 134 | `accepted` |
| `--MAX_FILE_HANDLES_FOR_READ_ENDS_MAP` | `-MAX_FILE_HANDLES` | `int` | `8000` | `no` | `optional` | — | `unsupported` |
| `--MAX_OPTICAL_DUPLICATE_SET_SIZE` | `` | `long` | `300000` | `no` | `optional` | — | `unsupported` |
| `--MAX_RECORDS_IN_RAM` | `` | `Integer` | `500000` | `no` | `common` | 149 | `accepted` |
| `--MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP` | `-MAX_SEQS` | `int` | `50000` | `no` | `optional` | — | `unsupported` |
| `--METRICS_FILE` | `-M` | `File` | `null` | `yes` | `required` | — | `unsupported` |
| `--MOLECULAR_IDENTIFIER_TAG` | `` | `String` | `null` | `no` | `optional` | — | `unsupported` |
| `--OPTICAL_DUPLICATE_PIXEL_DISTANCE` | `` | `int` | `100` | `no` | `optional` | — | `unsupported` |
| `--OUTPUT` | `-O` | `File` | `null` | `yes` | `required` | 136 | `accepted` |
| `--PROGRAM_GROUP_COMMAND_LINE` | `-PG_COMMAND` | `String` | `null` | `no` | `optional` | — | `unsupported` |
| `--PROGRAM_GROUP_NAME` | `-PG_NAME` | `String` | `MarkDuplicates` | `no` | `optional` | — | `unsupported` |
| `--PROGRAM_GROUP_VERSION` | `-PG_VERSION` | `String` | `null` | `no` | `optional` | — | `unsupported` |
| `--PROGRAM_RECORD_ID` | `-PG` | `String` | `MarkDuplicates` | `no` | `optional` | — | `unsupported` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 198 | `accepted` |
| `--READ_NAME_REGEX` | `` | `String` | `<optimized capture of last three ':' separated fields as numeric values>` | `no` | `optional` | — | `unsupported` |
| `--READ_ONE_BARCODE_TAG` | `` | `String` | `null` | `no` | `optional` | — | `unsupported` |
| `--READ_TWO_BARCODE_TAG` | `` | `String` | `null` | `no` | `optional` | — | `unsupported` |
| `--REFERENCE_SEQUENCE` | `-R` | `PicardHtsPath` | `null` | `no` | `common` | 138 | `accepted` |
| `--REMOVE_DUPLICATES` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--REMOVE_SEQUENCING_DUPLICATES` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--SORTING_COLLECTION_SIZE_RATIO` | `` | `double` | `0.25` | `no` | `optional` | — | `unsupported` |
| `--TAG_DUPLICATE_SET_MEMBERS` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--TAGGING_POLICY` | `` | `DuplicateTaggingPolicy` | `DontTag` | `no` | `optional` | — | `unsupported` |
| `--TMP_DIR` | `` | `List[File]` | `[]` | `no` | `common` | 146 | `accepted` |
| `--USE_JDK_DEFLATER` | `-use_jdk_deflater` | `Boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--USE_JDK_INFLATER` | `-use_jdk_inflater` | `Boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--VALIDATION_STRINGENCY` | `` | `ValidationStringency` | `STRICT` | `no` | `common` | — | `unsupported` |
| `--VERBOSITY` | `` | `LogLevel` | `INFO` | `no` | `common` | 202 | `accepted` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--reference` | 138 |
| `--metrics-file` | 140 |
| `--output-metrics` | 140 |
| `--output-manifest` | 143 |
| `--manifest` | 143 |
| `--tmp-dir` | 146 |
| `--max-records-in-memory` | 149 |
| `--read-name-regex` | 160 |
| `--tagging-policy` | 162 |
| `--add-pg-tag` | 164 |
| `--program-record-id` | 168 |
| `--program-group-name` | 170 |
| `--program-group-version` | 172 |
| `--program-group-command-line` | 174 |
| `--optical-duplicate-pixel-distance` | 176 |
| `--remove-duplicates` | 186 |
| `--remove-sequencing-duplicates` | 188 |
| `--resume-spill` | 190 |
| `--assume-sorted` | 192 |
| `--assume-sort-order` | 192 |
| `--create-output-bam-index` | 194 |
| `--disable-sequence-dictionary-validation` | 198 |
| `--use-jdk-deflater` | 199 |
| `--use-jdk-inflater` | 199 |
| `--clear-duplicates` | 200 |
| `--java-options` | 202 |

## `fastgatk-mutect2` ↔ `Mutect2`

- Source: `fastgatk-native/src/mutect2_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--active-probability-threshold` | `` | `double` | `0.002` | `no` | `advanced` | 1022 | `accepted` |
| `--adaptive-pruning-initial-error-rate` | `` | `double` | `0.001` | `no` | `advanced` | 1260 | `accepted` |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | 920 | `accepted` |
| `--af-of-alleles-not-in-resource` | `-default-af` | `double` | `-1.0` | `no` | `optional` | 971 | `accepted` |
| `--allele-informative-reads-overlap-margin` | `` | `int` | `2` | `no` | `advanced` | 1205 | `accepted` |
| `--alleles` | `` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--allow-non-unique-kmers-in-ref` | `` | `boolean` | `false` | `no` | `advanced` | 1293 | `accepted` |
| `--annotation` | `-A` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--annotation-group` | `-G` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--annotations-to-exclude` | `-AX` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--assembly-region-out` | `` | `String` | `null` | `no` | `optional` | 866 | `accepted` |
| `--assembly-region-padding` | `` | `int` | `100` | `no` | `optional` | 1031 | `accepted` |
| `--bam-output` | `-bamout` | `String` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--bam-writer-type` | `` | `WriterType` | `CALLED_HAPLOTYPES` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--base-qual-correction-factor` | `` | `int` | `5` | `no` | `advanced` | 1088 | `accepted` |
| `--base-quality-score-threshold` | `` | `byte` | `18` | `no` | `optional` | 1065 | `accepted` |
| `--callable-depth` | `` | `int` | `10` | `no` | `optional` | 926 | `accepted` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | 916 | `accepted` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--debug-assembly` | `-debug` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--disable-adaptive-pruning` | `` | `boolean` | `false` | `no` | `advanced` | 1250 | `accepted` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-cap-base-qualities-to-map-quality` | `` | `boolean` | `false` | `no` | `advanced` | 1106 | `accepted` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | 1353 | `accepted` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--disable-symmetric-hmm-normalizing` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--disable-tool-default-annotations` | `-disable-tool-default-annotations` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | 1358 | `accepted` |
| `--dont-increase-kmer-sizes-for-cycles` | `` | `boolean` | `false` | `no` | `advanced` | 1225 | `accepted` |
| `--dont-use-dragstr-pair-hmm-scores` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--dont-use-soft-clipped-bases` | `` | `boolean` | `false` | `no` | `optional` | 1154 | `accepted` |
| `--downsampling-stride` | `-stride` | `int` | `1` | `no` | `optional` | 1185 | `accepted` |
| `--dragstr-het-hom-ratio` | `` | `int` | `2` | `no` | `optional` | — | `ignored_via_catalog` |
| `--dragstr-params-path` | `` | `GATKPath` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--emit-ref-confidence` | `-ERC` | `ReferenceConfidenceMode` | `NONE` | `no` | `advanced` | 884 | `accepted` |
| `--enable-all-annotations` | `` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--enable-dynamic-read-disqualification-for-genotyping` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--expected-mismatch-rate-for-read-disqualification` | `` | `double` | `0.02` | `no` | `advanced` | 1143 | `accepted` |
| `--f1r2-max-depth` | `` | `int` | `200` | `no` | `optional` | — | `ignored_via_catalog` |
| `--f1r2-median-mq` | `` | `int` | `50` | `no` | `optional` | — | `ignored_via_catalog` |
| `--f1r2-min-bq` | `` | `int` | `20` | `no` | `optional` | — | `ignored_via_catalog` |
| `--f1r2-tar-gz` | `` | `File` | `null` | `no` | `optional` | 933 | `accepted` |
| `--flow-assembly-collapse-partial-mode` | `` | `boolean` | `false` | `no` | `advanced` | 1141 | `accepted` |
| `--flow-disallow-probs-larger-than-call` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--flow-fill-empty-bins-value` | `` | `double` | `0.001` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--flow-filter-alleles` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--flow-filter-alleles-qual-threshold` | `` | `float` | `30.0` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--flow-filter-alleles-sor-threshold` | `` | `float` | `3.0` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--flow-filter-lone-alleles` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--flow-lump-probs` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--flow-matrix-mods` | `` | `String` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--flow-mode` | `` | `FlowMode` | `NONE` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--flow-order-for-annotations` | `` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--flow-probability-scaling-factor` | `` | `int` | `10` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--flow-quantization-bins` | `` | `int` | `121` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--flow-remove-non-single-base-pair-indels` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--flow-remove-one-zero-probs` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--flow-report-insertion-or-deletion` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--flow-retain-max-n-probs-base-format` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--flow-symmetric-indel-probs` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--flow-use-t0-tag` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--force-active` | `` | `boolean` | `false` | `no` | `advanced` | 869 | `accepted` |
| `--force-call-filtered-alleles` | `-genotype-filtered-alleles` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--founder-id` | `-founder-id` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--genotype-germline-sites` | `` | `boolean` | `false` | `no` | `optional` | 964 | `accepted` |
| `--genotype-germline-sites-fraction` | `` | `double` | `1.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--genotype-pon-sites` | `` | `boolean` | `false` | `no` | `optional` | 960 | `accepted` |
| `--germline-resource` | `` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | 954 | `accepted` |
| `--graph-output` | `-graph` | `String` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--gvcf-lod-band` | `-LODB` | `List[Double]` | `[-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0]` | `no` | `advanced` | 905 | `accepted` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 725 | `accepted` |
| `--ignore-itr-artifacts` | `` | `boolean` | `false` | `no` | `optional` | 1158 | `accepted` |
| `--independent-mates` | `` | `boolean` | `false` | `no` | `advanced` | 1162 | `accepted` |
| `--initial-tumor-lod` | `-init-lod` | `double` | `2.0` | `no` | `optional` | 984 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `yes` | `required` | 820 | `accepted` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | 853 | `accepted` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 849 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--keep-boundary-flows` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--kmer-size` | `` | `List[Integer]` | `[10, 25]` | `no` | `advanced` | 1214 | `accepted` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--likelihood-calculation-engine` | `` | `Implementation` | `PairHMM` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--linked-de-bruijn-graph` | `` | `boolean` | `false` | `no` | `advanced` | 1254 | `accepted` |
| `--max-assembly-region-size` | `` | `int` | `300` | `no` | `optional` | 1037 | `accepted` |
| `--max-mnp-distance` | `-mnp-dist` | `int` | `1` | `no` | `advanced` | 1341 | `accepted` |
| `--max-num-haplotypes-in-population` | `` | `int` | `128` | `no` | `advanced` | 1322 | `accepted` |
| `--max-population-af` | `-max-af` | `double` | `0.01` | `no` | `optional` | 968 | `accepted` |
| `--max-prob-propagation-distance` | `` | `int` | `50` | `no` | `advanced` | 1040 | `accepted` |
| `--max-reads-per-alignment-start` | `` | `int` | `50` | `no` | `optional` | 1179 | `accepted` |
| `--max-suspicious-reads-per-alignment-start` | `` | `int` | `0` | `no` | `advanced` | 1192 | `accepted` |
| `--max-unpruned-variants` | `` | `int` | `100` | `no` | `advanced` | 1273 | `accepted` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--min-assembly-region-size` | `` | `int` | `50` | `no` | `optional` | 1034 | `accepted` |
| `--min-base-quality-score` | `-mbq` | `byte` | `10` | `no` | `optional` | 1096 | `accepted` |
| `--min-dangling-branch-length` | `` | `int` | `4` | `no` | `advanced` | 1279 | `accepted` |
| `--min-pruning` | `` | `int` | `2` | `no` | `advanced` | 1236 | `accepted` |
| `--minimum-allele-fraction` | `-min-AF` | `double` | `0.0` | `no` | `advanced` | 941 | `accepted` |
| `--mitochondria-mode` | `` | `Boolean` | `false` | `no` | `optional` | 1363 | `accepted` |
| `--native-pair-hmm-threads` | `` | `int` | `4` | `no` | `optional` | 1008 | `accepted` |
| `--native-pair-hmm-use-double-precision` | `` | `boolean` | `false` | `no` | `optional` | 1012 | `accepted` |
| `--normal-lod` | `` | `double` | `2.2` | `no` | `optional` | 994 | `accepted` |
| `--normal-sample` | `-normal` | `List[String]` | `[]` | `no` | `optional` | 835 | `accepted` |
| `--num-pruning-samples` | `` | `int` | `1` | `no` | `advanced` | 1242 | `accepted` |
| `--output` | `-O` | `GATKPath` | `null` | `yes` | `required` | 864 | `accepted` |
| `--pair-hmm-gap-continuation-penalty` | `` | `int` | `10` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--pair-hmm-implementation` | `-pairHMM` | `Implementation` | `FASTEST_AVAILABLE` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--pair-hmm-results-file` | `` | `GATKPath` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--panel-of-normals` | `-pon` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | 957 | `accepted` |
| `--pcr-indel-model` | `` | `PCRErrorModel` | `CONSERVATIVE` | `no` | `advanced` | 1116 | `accepted` |
| `--pcr-indel-qual` | `` | `int` | `40` | `no` | `optional` | 1080 | `accepted` |
| `--pcr-snv-qual` | `` | `int` | `40` | `no` | `optional` | 1073 | `accepted` |
| `--pedigree` | `-ped` | `GATKPath` | `null` | `no` | `optional` | — | `unsupported` |
| `--permutect-alt-downsample` | `` | `int` | `20` | `no` | `optional` | — | `ignored_via_catalog` |
| `--permutect-dataset-mode` | `` | `PermutectDatasetMode` | `ILLUMINA` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--permutect-non-artifact-ratio` | `` | `int` | `1` | `no` | `optional` | — | `ignored_via_catalog` |
| `--permutect-ref-downsample` | `` | `int` | `10` | `no` | `optional` | — | `ignored_via_catalog` |
| `--permutect-test-dataset` | `` | `File` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--permutect-test-truth` | `` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--permutect-training-dataset` | `` | `File` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--permutect-training-truth` | `` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--phred-scaled-global-read-mismapping-rate` | `` | `int` | `45` | `no` | `advanced` | 1108 | `accepted` |
| `--pileup-detection` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--pruning-lod-threshold` | `` | `double` | `2.302585092994046` | `no` | `advanced` | 1263 | `accepted` |
| `--pruning-seeding-lod-threshold` | `` | `double` | `9.210340371976184` | `no` | `advanced` | 1268 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | 1349 | `accepted` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--recover-all-dangling-branches` | `` | `boolean` | `false` | `no` | `advanced` | 1295 | `accepted` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 847 | `accepted` |
| `--reference-model-deletion-quality` | `` | `byte` | `30` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | 873 | `accepted` |
| `--smith-waterman` | `` | `Implementation` | `FASTEST_AVAILABLE` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-dangling-end-gap-extend-penalty` | `` | `int` | `-6` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-dangling-end-gap-open-penalty` | `` | `int` | `-110` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-dangling-end-match-value` | `` | `int` | `25` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-dangling-end-mismatch-penalty` | `` | `int` | `-50` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-haplotype-to-reference-gap-extend-penalty` | `` | `int` | `-11` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-haplotype-to-reference-gap-open-penalty` | `` | `int` | `-260` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-haplotype-to-reference-match-value` | `` | `int` | `200` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-haplotype-to-reference-mismatch-penalty` | `` | `int` | `-150` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-read-to-haplotype-gap-extend-penalty` | `` | `int` | `-5` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-read-to-haplotype-gap-open-penalty` | `` | `int` | `-30` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-read-to-haplotype-match-value` | `` | `int` | `10` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--smith-waterman-read-to-haplotype-mismatch-penalty` | `` | `int` | `-15` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--soft-clip-low-quality-ends` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--tumor-lod-to-emit` | `-emit-lod` | `double` | `3.0` | `no` | `optional` | 989 | `accepted` |
| `--tumor-sample` | `-tumor` | `String` | `null` | `no` | `optional` | 831 | `accepted` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-pdhmm` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--use-pdhmm-overlap-optimization` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--annotate-allele-specific` | 825 |
| `--normal-input` | 833 |
| `--normal` | 835 |
| `--region` | 849 |
| `--stats` | 924 |
| `--contamination` | 935 |
| `--contamination-fraction` | 935 |
| `--somatic-prior` | 948 |
| `--germline-prior` | 950 |
| `--artifact-prior` | 952 |
| `--population-af` | 972 |
| `--default-af` | 973 |
| `--output-manifest` | 997 |
| `--manifest` | 997 |
| `--batch-records` | 1000 |
| `--stream-by-region` | 1002 |
| `--threads` | 1008 |
| `--min-depth` | 1018 |
| `--min-alt-support` | 1020 |
| `--minimum-mapping-quality` | 1043 |
| `--min-read-length` | 1049 |
| `--max-read-length` | 1057 |
| `--min-base-quality` | 1097 |
| `--flow-assembly-collapse-hmer-size` | 1132 |
| `--include-duplicates` | 1152 |
| `--do-not-correct-overlapping-quality` | 1168 |
| `--do-not-correct-overlapping-base-qualities` | 1170 |
| `--max-reads-per-locus` | 1178 |
| `--downsampling-seed` | 1202 |
| `--min-kmer-count` | 1230 |
| `--adaptive-pruning` | 1248 |
| `--disable-artificial-haplotype-recovery` | 1256 |
| `--enable-legacy-graph-cycle-detection` | 1258 |
| `--min-dangling-matching-bases` | 1286 |
| `--error-correct-reads` | 1300 |
| `--error-correction-log-odds` | 1304 |
| `--kmer-length-for-read-error-correction` | 1307 |
| `--min-observations-for-kmer-to-be-solid` | 1314 |
| `--max-haplotype-paths` | 1321 |
| `--max-haplotype-depth` | 1329 |
| `--max-haplotype-combination-alleles` | 1334 |
| `--mnp-dist` | 1341 |

## `fastgatk-combine-gvcfs` ↔ `CombineGVCFs`

- Source: `fastgatk-native/src/combine_gvcf_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_CombineGVCFs.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--annotation` | `-A` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--annotation-group` | `-G` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--annotations-to-exclude` | `-AX` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--break-bands-at-multiples-of` | `` | `int` | `0` | `no` | `optional` | 185 | `accepted` |
| `--call-genotypes` | `` | `boolean` | `false` | `no` | `optional` | 173 | `accepted` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--combine-variants-distance` | `` | `int` | `0` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--convert-to-base-pair-resolution` | `` | `boolean` | `false` | `no` | `optional` | 179 | `accepted` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | 163 | `accepted` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--dbsnp` | `-D` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | 194 | `accepted` |
| `--disable-tool-default-annotations` | `-disable-tool-default-annotations` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--drop-somatic-filtering-annotations` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--enable-all-annotations` | `` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--flow-order-for-annotations` | `` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--founder-id` | `-founder-id` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 119 | `accepted` |
| `--ignore-variants-starting-outside-interval` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--input-is-somatic` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | 149 | `accepted` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 142 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-distance` | `` | `int` | `2147483647` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--output` | `-O` | `GATKPath` | `null` | `yes` | `required` | 158 | `accepted` |
| `--pedigree` | `-ped` | `GATKPath` | `null` | `no` | `optional` | — | `unsupported` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 194 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--ref-padding` | `` | `int` | `1` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 140 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | 167 | `accepted` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--variant` | `-V` | `List[GATKPath]` | `[]` | `yes` | `required` | 138 | `accepted` |
| `--variant-output-filtering` | `` | `Mode` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | 196 | `accepted` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--interval` | 142 |
| `--region` | 143 |
| `--output-manifest` | 160 |
| `--manifest` | 160 |
| `--fastgatk-materialize-genotypes` | 171 |
| `--stream-merge` | 183 |
| `--stream-by-locus` | 183 |
| `--java-options` | 196 |

## `fastgatk-filter-mutect-calls` ↔ `FilterMutectCalls`

- Source: `fastgatk-native/src/filter_mutect_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_mutect_filtering_FilterMutectCalls.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--contamination-estimate` | `` | `double` | `0.0` | `no` | `optional` | 532 | `accepted` |
| `--contamination-table` | `` | `List[File]` | `[]` | `no` | `optional` | 540 | `accepted` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | 682 | `accepted` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--distance-on-haplotype` | `` | `int` | `100` | `no` | `optional` | 664 | `accepted` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--f-score-beta` | `` | `double` | `1.0` | `no` | `optional` | 655 | `accepted` |
| `--false-discovery-rate` | `` | `double` | `0.05` | `no` | `optional` | 653 | `accepted` |
| `--filtering-stats` | `` | `String` | `null` | `no` | `optional` | 516 | `accepted` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 431 | `accepted` |
| `--initial-threshold` | `` | `double` | `0.1` | `no` | `optional` | 651 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | 494 | `accepted` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | 505 | `accepted` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 498 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--log-artifact-prior` | `` | `double` | `-2.302585092994046` | `no` | `optional` | 612 | `accepted` |
| `--log-indel-prior` | `` | `double` | `-16.11809565095832` | `no` | `optional` | 608 | `accepted` |
| `--log-snv-prior` | `` | `double` | `-13.815510557964275` | `no` | `optional` | 604 | `accepted` |
| `--long-indel-length` | `` | `int` | `5` | `no` | `optional` | 601 | `accepted` |
| `--max-alt-allele-count` | `` | `int` | `1` | `no` | `optional` | 561 | `accepted` |
| `--max-events-in-haplotype` | `` | `int` | `2` | `no` | `optional` | 660 | `accepted` |
| `--max-events-in-region` | `` | `int` | `3` | `no` | `optional` | 657 | `accepted` |
| `--max-median-fragment-length-difference` | `` | `int` | `10000` | `no` | `optional` | 575 | `accepted` |
| `--max-n-ratio` | `` | `double` | `Infinity` | `no` | `optional` | 578 | `accepted` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--microbial-mode` | `` | `boolean` | `false` | `no` | `optional` | 673 | `accepted` |
| `--min-allele-fraction` | `` | `double` | `0.0` | `no` | `optional` | 548 | `accepted` |
| `--min-median-base-quality` | `` | `int` | `20` | `no` | `optional` | 564 | `accepted` |
| `--min-median-mapping-quality` | `` | `int` | `-1` | `no` | `optional` | 567 | `accepted` |
| `--min-median-read-position` | `` | `int` | `1` | `no` | `optional` | 572 | `accepted` |
| `--min-reads-per-strand` | `` | `int` | `0` | `no` | `optional` | 552 | `accepted` |
| `--min-slippage-length` | `` | `int` | `8` | `no` | `optional` | 589 | `accepted` |
| `--mitochondria-mode` | `` | `boolean` | `false` | `no` | `optional` | 669 | `accepted` |
| `--normal-p-value-threshold` | `` | `double` | `0.001` | `no` | `optional` | 581 | `accepted` |
| `--orientation-bias-artifact-priors` | `-ob-priors` | `List[File]` | `[]` | `no` | `optional` | 630 | `accepted` |
| `--output` | `-O` | `String` | `null` | `yes` | `required` | 496 | `accepted` |
| `--pcr-slippage-rate` | `` | `double` | `0.1` | `no` | `optional` | 596 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 688 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 490 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | 518 | `accepted` |
| `--stats` | `` | `String` | `null` | `no` | `optional` | 514 | `accepted` |
| `--threshold-strategy` | `` | `Strategy` | `OPTIMAL_F_SCORE` | `no` | `optional` | 645 | `accepted` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--tumor-segmentation` | `` | `List[File]` | `[]` | `no` | `optional` | 542 | `accepted` |
| `--unique-alt-read-count` | `-unique` | `int` | `0` | `no` | `optional` | 558 | `accepted` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--variant` | `-V` | `GATKPath` | `null` | `yes` | `required` | 492 | `accepted` |
| `--variant-output-filtering` | `` | `Mode` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | 689 | `accepted` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--interval` | 499 |
| `--region` | 499 |
| `--min-tlod` | 529 |
| `--contamination-fraction` | 531 |
| `--max-contamination-probability` | 537 |
| `--tumor-sample` | 544 |
| `--min-orientation-balance` | 546 |
| `--min-af` | 548 |
| `--min-reads-on-each-strand` | 553 |
| `--unique` | 558 |
| `--normal-pileup-p-value-threshold` | 580 |
| `--filter-error-probability-threshold` | 586 |
| `--min-pcr-slippage-size` | 590 |
| `--slippage-rate` | 595 |
| `--log-somatic-prior` | 614 |
| `--min-somatic-probability` | 618 |
| `--max-germline-probability` | 621 |
| `--max-germline-posterior` | 622 |
| `--max-artifact-probability` | 627 |
| `--ob-priors` | 631 |
| `--max-strand-artifact-probability` | 636 |
| `--max-orientation-artifact-probability` | 637 |
| `--normal-artifact-lod` | 642 |
| `--max-intra-haplotype-distance` | 663 |
| `--output-manifest` | 677 |
| `--manifest` | 677 |
| `--java-options` | 688 |

## `fastgatk-genomicsdb-import` ↔ `GenomicsDBImport`

- Source: `fastgatk-native/src/genomicsdb_import_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_genomicsdb_GenomicsDBImport.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--avoid-nio` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--batch-size` | `` | `int` | `0` | `no` | `optional` | 235 | `accepted` |
| `--bypass-feature-reader` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `0` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `0` | `no` | `optional` | — | `unsupported` |
| `--consolidate` | `` | `Boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--genomicsdb-segment-size` | `` | `long` | `1048576` | `no` | `optional` | — | `ignored_via_catalog` |
| `--genomicsdb-shared-posixfs-optimizations` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--genomicsdb-update-workspace-path` | `` | `String` | `null` | `yes` | `required` | 168 | `accepted` |
| `--genomicsdb-use-gcs-hdfs-connector` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--genomicsdb-vcf-buffer-size` | `` | `long` | `16384` | `no` | `optional` | — | `ignored_via_catalog` |
| `--genomicsdb-workspace-path` | `` | `String` | `null` | `yes` | `required` | 156 | `accepted` |
| `--header` | `` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 141 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | 214 | `accepted` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 200 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-num-intervals-to-import-in-parallel` | `` | `int` | `1` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--merge-contigs-into-num-partitions` | `-merge-contigs-into-num-partitions` | `int` | `0` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--merge-input-intervals` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--output-interval-list-to-file` | `` | `String` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--overwrite-existing-genomicsdb-workspace` | `` | `Boolean` | `false` | `no` | `optional` | 267 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reader-threads` | `` | `int` | `1` | `no` | `advanced` | 243 | `accepted` |
| `--reference` | `-R` | `GATKPath` | `null` | `no` | `optional` | — | `unsupported` |
| `--sample-name-map` | `` | `String` | `null` | `no` | `advanced` | 192 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | 251 | `accepted` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--validate-sample-name-map` | `` | `Boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--variant` | `-V` | `List[String]` | `[]` | `no` | `optional` | 181 | `accepted` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--interval` | 201 |
| `--region` | 201 |
| `--output-manifest` | 230 |
| `--manifest` | 230 |
| `--fastgatk-native-workspace` | 258 |
| `--resume-native-workspace` | 262 |

## `fastgatk-gather-bqsr-reports` ↔ `GatherBQSRReports`

- Source: `fastgatk-native/src/gather_bqsr_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_bqsr_GatherBQSRReports.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 62 | `accepted` |
| `--input` | `-I` | `List[File]` | `[]` | `yes` | `required` | 71 | `accepted` |
| `--output` | `-O` | `File` | `null` | `yes` | `required` | 75 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 81 | `accepted` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | 85 | `accepted` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | 94 | `accepted` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | 95 | `accepted` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | 87 | `accepted` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--output-manifest` | 79 |
| `--jdk-deflater` | 94 |
| `--jdk-inflater` | 95 |

## `fastgatk-analyze-covariates` ↔ `AnalyzeCovariates`

- Source: `fastgatk-native/src/analyze_covariates_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_bqsr_AnalyzeCovariates.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--after-report-file` | `-after` | `File` | `null` | `no` | `optional` | 181 | `accepted` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--before-report-file` | `-before` | `File` | `null` | `no` | `optional` | 178 | `accepted` |
| `--bqsr-recal-file` | `-bqsr` | `File` | `null` | `no` | `optional` | 190 | `accepted` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 168 | `accepted` |
| `--ignore-last-modification-times` | `` | `boolean` | `false` | `no` | `optional` | 203 | `accepted` |
| `--intermediate-csv-file` | `-csv` | `File` | `null` | `no` | `optional` | 194 | `accepted` |
| `--plots-report-file` | `-plots` | `File` | `null` | `no` | `optional` | 197 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 205 | `accepted` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--disable-sequence-dictionary-validation` | 205 |

## `fastgatk-variants-to-table` ↔ `VariantsToTable`

- Source: `fastgatk-native/src/variants_to_table_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_variantutils_VariantsToTable.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--asFieldsToTake` | `-ASF` | `List[String]` | `[]` | `no` | `optional` | 184 | `accepted` |
| `--asGenotypeFieldsToTake` | `-ASGF` | `List[String]` | `[]` | `no` | `optional` | 190 | `accepted` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | 203 | `accepted` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--error-if-missing-data` | `-EMD` | `boolean` | `false` | `no` | `advanced` | 199 | `accepted` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | 171 | `accepted` |
| `--fields` | `-F` | `List[String]` | `[]` | `no` | `optional` | 180 | `accepted` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--genotype-fields` | `-GF` | `List[String]` | `[]` | `no` | `optional` | 182 | `accepted` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 107 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | 177 | `accepted` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | 151 | `accepted` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | 174 | `accepted` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | 142 | `accepted` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 135 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--moltenize` | `-moltenize` | `boolean` | `false` | `no` | `advanced` | 198 | `accepted` |
| `--output` | `-O` | `String` | `null` | `yes` | `required` | 131 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 203 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `no` | `optional` | 133 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--show-filtered` | `-raw` | `boolean` | `false` | `no` | `advanced` | 197 | `accepted` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--split-multi-allelic` | `-SMA` | `boolean` | `false` | `no` | `optional` | 196 | `accepted` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--variant` | `-V` | `GATKPath` | `null` | `yes` | `required` | 129 | `accepted` |
| `--variant-output-filtering` | `` | `Mode` | `null` | `no` | `advanced` | 159 | `accepted` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | 206 | `accepted` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--interval` | 136 |
| `--region` | 136 |
| `--allele-specific-fields` | 184 |
| `--allele-specific-genotype-fields` | 190 |
| `--raw` | 197 |
| `--output-manifest` | 200 |
| `--manifest` | 200 |
| `--java-options` | 206 |

## `fastgatk-variant-eval` ↔ `VariantEval`

- Source: `fastgatk-native/src/variant_eval_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_varianteval_VariantEval.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--ancestral-alignments` | `-aa` | `File` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--combine-variants-distance` | `` | `int` | `0` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--comparison` | `-comp` | `List[FeatureInput[VariantContext]]` | `[]` | `no` | `optional` | 107 | `accepted` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--dbsnp` | `-D` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | 164 | `accepted` |
| `--do-not-use-all-standard-modules` | `-no-ev` | `Boolean` | `false` | `no` | `optional` | 150 | `accepted` |
| `--do-not-use-all-standard-stratifications` | `-no-st` | `Boolean` | `false` | `no` | `optional` | 152 | `accepted` |
| `--eval` | `-eval` | `List[FeatureInput[VariantContext]]` | `[]` | `yes` | `required` | 105 | `accepted` |
| `--eval-module` | `-EV` | `List[String]` | `[]` | `no` | `optional` | 141 | `accepted` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--gold-standard` | `-gold` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | 176 | `accepted` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 82 | `accepted` |
| `--ignore-variants-starting-outside-interval` | `` | `boolean` | `false` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | — | `ignored_via_catalog` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 125 | `accepted` |
| `--keep-ac0` | `-keep-ac0` | `boolean` | `false` | `no` | `optional` | 156 | `accepted` |
| `--known-cnvs` | `-known-cnvs` | `FeatureInput[Feature]` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--knownNames` | `-known-name` | `Set[String]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--list` | `-ls` | `Boolean` | `false` | `no` | `optional` | 102 | `accepted` |
| `--max-distance` | `` | `int` | `2147483647` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--mendelian-violation-qual-threshold` | `-mvq` | `double` | `50.0` | `no` | `optional` | 116 | `accepted` |
| `--merge-evals` | `-merge-evals` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--min-phase-quality` | `-mpq` | `double` | `10.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--output` | `-O` | `File` | `null` | `yes` | `required` | 109 | `accepted` |
| `--pedigree` | `-ped` | `GATKPath` | `null` | `no` | `optional` | 113 | `accepted` |
| `--pedigreeValidationType` | `-pedValidationType` | `PedigreeValidationType` | `STRICT` | `no` | `optional` | — | `ignored_via_catalog` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 164 | `accepted` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--ref-padding` | `` | `int` | `1` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 111 | `accepted` |
| `--require-strict-allele-match` | `-strict` | `boolean` | `false` | `no` | `optional` | 154 | `accepted` |
| `--sample` | `-sn` | `Set[String]` | `[]` | `no` | `optional` | 177 | `accepted` |
| `--sample-ploidy` | `-ploidy` | `int` | `2` | `no` | `optional` | 179 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--selectExps` | `-select` | `ArrayList[String]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--selectNames` | `-select-name` | `ArrayList[String]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--strat-intervals` | `-strat-intervals` | `FeatureInput[Feature]` | `null` | `no` | `optional` | 178 | `accepted` |
| `--stratification-module` | `-ST` | `List[String]` | `[]` | `no` | `optional` | 144 | `accepted` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--variant-output-filtering` | `` | `Mode` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--gatk-report` | 103 |
| `--gatk-compatible-report` | 103 |
| `--ped` | 113 |
| `-S` | 143 |
| `--known-name` | 176 |
| `--num-samples` | 179 |

## `fastgatk-validate-variants` ↔ `ValidateVariants`

- Source: `fastgatk-native/src/validate_variants_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_variantutils_ValidateVariants.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--dbsnp` | `-D` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | 97 | `accepted` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | 190 | `accepted` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--do-not-validate-filtered-records` | `-do-not-validate-filtered-records` | `Boolean` | `false` | `no` | `optional` | 146 | `accepted` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | 119 | `accepted` |
| `--fail-gvcf-on-overlap` | `-no-overlaps` | `Boolean` | `false` | `no` | `optional` | 177 | `accepted` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 74 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | 131 | `accepted` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | 122 | `accepted` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | 109 | `accepted` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 102 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 187 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `no` | `optional` | 95 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--validate-GVCF` | `-gvcf` | `Boolean` | `false` | `no` | `optional` | 167 | `accepted` |
| `--validation-type-to-exclude` | `-Xtype` | `List[ValidationType]` | `[]` | `no` | `optional` | 141 | `accepted` |
| `--variant` | `-V` | `GATKPath` | `null` | `yes` | `required` | 93 | `accepted` |
| `--variant-output-filtering` | `` | `Mode` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--warn-on-errors` | `-warn-on-errors` | `Boolean` | `false` | `no` | `optional` | 159 | `accepted` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `-O` | 100 |
| `--do-not-validate-filtered` | 147 |
| `--gvcf` | 167 |
| `--no-overlaps` | 177 |

## `fastgatk-get-pileup-summaries` ↔ `GetPileupSummaries`

- Source: `fastgatk-native/src/get_pileup_summaries_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_contamination_GetPileupSummaries.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | 303 | `accepted` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | 293 | `accepted` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | 296 | `accepted` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | 239 | `accepted` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 205 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `yes` | `required` | 224 | `accepted` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | 254 | `accepted` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | 251 | `accepted` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | 242 | `accepted` |
| `--intervals` | `-L` | `List[String]` | `[]` | `yes` | `required` | 232 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | 306 | `accepted` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-depth-per-sample` | `-max-depth-per-sample` | `int` | `0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--maximum-population-allele-frequency` | `-max-af` | `double` | `0.2` | `no` | `optional` | 261 | `accepted` |
| `--minimum-population-allele-frequency` | `-min-af` | `double` | `0.01` | `no` | `optional` | 257 | `accepted` |
| `--output` | `-O` | `File` | `null` | `yes` | `required` | 228 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 291 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | 301 | `accepted` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `no` | `optional` | 230 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | 310 | `accepted` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--variant` | `-V` | `FeatureInput[VariantContext]` | `null` | `yes` | `required` | 226 | `accepted` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | 309 | `accepted` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--interval` | 233 |
| `--region` | 233 |
| `--min-af` | 257 |
| `--max-af` | 261 |
| `--minimum-mapping-quality` | 265 |
| `--min-read-length` | 268 |
| `--max-read-length` | 275 |
| `--batch-records` | 283 |
| `--threads` | 286 |
| `--output-manifest` | 288 |
| `--manifest` | 288 |
| `--include-duplicates` | 299 |
| `--java-options` | 309 |

## `fastgatk-calculate-contamination` ↔ `CalculateContamination`

- Source: `fastgatk-native/src/calculate_contamination_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_contamination_CalculateContamination.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 83 | `accepted` |
| `--high-coverage-ratio-threshold` | `` | `double` | `3.0` | `no` | `optional` | 113 | `accepted` |
| `--input` | `-I` | `File` | `null` | `yes` | `required` | 93 | `accepted` |
| `--low-coverage-ratio-threshold` | `` | `double` | `0.5` | `no` | `optional` | 110 | `accepted` |
| `--matched-normal` | `-matched` | `File` | `null` | `no` | `optional` | 95 | `accepted` |
| `--output` | `-O` | `File` | `null` | `yes` | `required` | 101 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 121 | `accepted` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--tumor-segmentation` | `-segments` | `File` | `null` | `no` | `optional` | 103 | `accepted` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | 123 | `accepted` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--matched` | 95 |
| `--segments` | 103 |
| `--threads` | 116 |
| `--output-manifest` | 118 |
| `--manifest` | 118 |
| `--disable-sequence-dictionary-validation` | 121 |
| `--java-options` | 123 |
| `--seconds-between-progress-updates` | 124 |

## `fastgatk-gather-pileup-summaries` ↔ `GatherPileupSummaries`

- Source: `fastgatk-native/src/gather_pileup_summaries_tool.cpp`
- GATK reference: `<none>` (native-only bridge binary)
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

No gatkdoc JSON is shipped with this GATK version; the tool is native-only.

## `fastgatk-learn-read-orientation-model` ↔ `LearnReadOrientationModel`

- Source: `fastgatk-native/src/learn_read_orientation_model_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_readorientation_LearnReadOrientationModel.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--convergence-threshold` | `` | `double` | `1.0E-4` | `no` | `optional` | 136 | `accepted` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 116 | `accepted` |
| `--input` | `-I` | `List[File]` | `[]` | `yes` | `required` | 130 | `accepted` |
| `--max-depth` | `` | `int` | `200` | `no` | `optional` | 144 | `accepted` |
| `--num-em-iterations` | `` | `int` | `20` | `no` | `optional` | 139 | `accepted` |
| `--output` | `-O` | `File` | `null` | `yes` | `required` | 132 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 151 | `accepted` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | 155 | `accepted` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | 157 | `accepted` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | 158 | `accepted` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | 167 | `accepted` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--sample` | 134 |
| `--max-em-iterations` | 139 |
| `--threads` | 146 |
| `--output-manifest` | 148 |
| `--manifest` | 148 |
| `--disable-sequence-dictionary-validation` | 164 |
| `--java-options` | 166 |

## `fastgatk-denoise-read-counts` ↔ `DenoiseReadCounts`

- Source: `fastgatk-native/src/denoise_read_counts_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_copynumber_DenoiseReadCounts.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--annotated-intervals` | `` | `File` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--count-panel-of-normals` | `` | `File` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--denoised-copy-ratios` | `` | `File` | `null` | `yes` | `required` | — | `ignored_via_catalog` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 92 | `accepted` |
| `--input` | `-I` | `File` | `null` | `yes` | `required` | 105 | `accepted` |
| `--number-of-eigensamples` | `` | `Integer` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 134 | `accepted` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--standardized-copy-ratios` | `` | `File` | `null` | `yes` | `required` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `-O` | 107 |
| `--disable-sequence-dictionary-validation` | 134 |
| `--disable-tool-default-read-filters` | 135 |
| `--dont-trim-intervals` | 136 |

## `fastgatk-create-read-count-panel-of-normals` ↔ `CreateReadCountPanelOfNormals`

- Source: `fastgatk-native/src/create_read_count_panel_of_normals_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_copynumber_CreateReadCountPanelOfNormals.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--annotated-intervals` | `` | `File` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--conf` | `` | `List[String]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--do-impute-zeros` | `` | `boolean` | `true` | `no` | `optional` | — | `ignored_via_catalog` |
| `--extreme-outlier-truncation-percentile` | `` | `double` | `0.1` | `no` | `optional` | — | `ignored_via_catalog` |
| `--extreme-sample-median-percentile` | `` | `double` | `2.5` | `no` | `optional` | — | `ignored_via_catalog` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 81 | `accepted` |
| `--input` | `-I` | `List[File]` | `[]` | `yes` | `required` | 97 | `accepted` |
| `--maximum-chunk-size` | `` | `int` | `16777215` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--maximum-zeros-in-interval-percentage` | `` | `double` | `5.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--maximum-zeros-in-sample-percentage` | `` | `double` | `5.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--minimum-interval-median-percentile` | `` | `double` | `10.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--number-of-eigensamples` | `` | `int` | `20` | `no` | `optional` | — | `ignored_via_catalog` |
| `--output` | `-O` | `File` | `null` | `yes` | `required` | 99 | `accepted` |
| `--program-name` | `` | `String` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 130 | `accepted` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--spark-master` | `` | `String` | `local[*]` | `no` | `optional` | — | `unsupported` |
| `--spark-verbosity` | `` | `String` | `null` | `no` | `optional` | — | `unsupported` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

## `fastgatk-call-copy-ratio-segments` ↔ `CallCopyRatioSegments`

- Source: `fastgatk-native/src/call_copy_ratio_segments_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_copynumber_CallCopyRatioSegments.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--calling-copy-ratio-z-score-threshold` | `` | `double` | `2.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 93 | `accepted` |
| `--input` | `-I` | `File` | `null` | `yes` | `required` | 104 | `accepted` |
| `--neutral-segment-copy-ratio-lower-bound` | `` | `double` | `0.9` | `no` | `optional` | — | `ignored_via_catalog` |
| `--neutral-segment-copy-ratio-upper-bound` | `` | `double` | `1.1` | `no` | `optional` | — | `ignored_via_catalog` |
| `--outlier-neutral-segment-copy-ratio-z-score-threshold` | `` | `double` | `2.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--output` | `-O` | `File` | `null` | `yes` | `required` | 106 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 133 | `accepted` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--disable-sequence-dictionary-validation` | 133 |
| `--disable-tool-default-read-filters` | 134 |

## `fastgatk-model-segments` ↔ `ModelSegments`

- Source: `fastgatk-native/src/model_segments_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_copynumber_ModelSegments.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--allelic-counts` | `` | `List[File]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--denoised-copy-ratios` | `` | `List[File]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--genotyping-base-error-rate` | `` | `double` | `0.05` | `no` | `optional` | — | `ignored_via_catalog` |
| `--genotyping-homozygous-log-ratio-threshold` | `` | `double` | `-10.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 243 | `accepted` |
| `--kernel-approximation-dimension` | `` | `int` | `100` | `no` | `optional` | — | `ignored_via_catalog` |
| `--kernel-scaling-allele-fraction` | `` | `double` | `1.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--kernel-variance-allele-fraction` | `` | `double` | `0.025` | `no` | `optional` | — | `ignored_via_catalog` |
| `--kernel-variance-copy-ratio` | `` | `double` | `0.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--maximum-number-of-segments-per-chromosome` | `` | `int` | `1000` | `no` | `optional` | — | `ignored_via_catalog` |
| `--maximum-number-of-smoothing-iterations` | `` | `int` | `25` | `no` | `optional` | — | `ignored_via_catalog` |
| `--minimum-total-allele-count-case` | `` | `int` | `0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--minimum-total-allele-count-normal` | `` | `int` | `30` | `no` | `optional` | — | `ignored_via_catalog` |
| `--minor-allele-fraction-prior-alpha` | `` | `double` | `25.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--normal-allelic-counts` | `` | `File` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--number-of-burn-in-samples-allele-fraction` | `` | `int` | `50` | `no` | `optional` | — | `ignored_via_catalog` |
| `--number-of-burn-in-samples-copy-ratio` | `` | `int` | `50` | `no` | `optional` | — | `ignored_via_catalog` |
| `--number-of-changepoints-penalty-factor` | `` | `double` | `1.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--number-of-samples-allele-fraction` | `` | `int` | `100` | `no` | `optional` | — | `ignored_via_catalog` |
| `--number-of-samples-copy-ratio` | `` | `int` | `100` | `no` | `optional` | — | `ignored_via_catalog` |
| `--number-of-smoothing-iterations-per-fit` | `` | `int` | `0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--output` | `-O` | `File` | `null` | `yes` | `required` | 306 | `accepted` |
| `--output-prefix` | `` | `String` | `null` | `yes` | `required` | — | `ignored_via_catalog` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 419 | `accepted` |
| `--segments` | `` | `File` | `null` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--smoothing-credible-interval-threshold-allele-fraction` | `` | `double` | `2.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--smoothing-credible-interval-threshold-copy-ratio` | `` | `double` | `2.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--window-size` | `` | `List[Integer]` | `[8, 16, 32, 64, 128, 256]` | `no` | `optional` | — | `ignored_via_catalog` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `-I` | 280 |
| `--disable-sequence-dictionary-validation` | 419 |
| `--disable-tool-default-read-filters` | 420 |

## `fastgatk-gather-tranches` ↔ `GatherTranches`

- Source: `fastgatk-native/src/gather_tranches_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_vqsr_GatherTranches.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 135 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `yes` | `required` | 144 | `accepted` |
| `--mode` | `-mode` | `Mode` | `null` | `yes` | `required` | — | `ignored_via_catalog` |
| `--output` | `-O` | `GATKPath` | `null` | `yes` | `required` | 160 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 165 | `accepted` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--truth-sensitivity-tranche` | `-tranche` | `List[Double]` | `[100.0, 99.9, 99.0, 90.0]` | `no` | `optional` | 146 | `accepted` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--disable-sequence-dictionary-validation` | 165 |

## `fastgatk-annotate-intervals` ↔ `AnnotateIntervals`

- Source: `fastgatk-native/src/annotate_intervals_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_copynumber_AnnotateIntervals.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | 114 | `accepted` |
| `--feature-query-lookahead` | `` | `int` | `1000000` | `no` | `optional` | — | `ignored_via_catalog` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 90 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | — | `ignored_via_catalog` |
| `--intervals` | `-L` | `List[String]` | `[]` | `yes` | `required` | 108 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--mappability-track` | `` | `FeatureInput[BEDFeature]` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--output` | `-O` | `File` | `null` | `yes` | `required` | 117 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 164 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 106 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--segmental-duplication-track` | `` | `FeatureInput[BEDFeature]` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

## `fastgatk-count-bases-in-reference` ↔ `CountBasesInReference`

- Source: `fastgatk-native/src/count_bases_in_reference_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_fasta_CountBasesInReference.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 58 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | — | `ignored_via_catalog` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 69 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--output` | `-O` | `GATKPath` | `null` | `no` | `optional` | 79 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 87 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 66 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--interval` | 69 |
| `--region` | 70 |
| `--threads` | 77 |
| `--output-manifest` | 82 |
| `--manifest` | 82 |

## `fastgatk-compare-references` ↔ `CompareReferences`

- Source: `fastgatk-native/src/compare_references_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_reference_CompareReferences.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--base-comparison` | `` | `BaseComparisonMode` | `NO_BASE_COMPARISON` | `no` | `optional` | 84 | `accepted` |
| `--base-comparison-output` | `` | `GATKPath` | `null` | `no` | `optional` | 86 | `accepted` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--display-only-differing-sequences` | `` | `boolean` | `false` | `no` | `optional` | 82 | `accepted` |
| `--display-sequences-by-name` | `` | `boolean` | `false` | `no` | `optional` | 80 | `accepted` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 55 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | — | `ignored_via_catalog` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--md5-calculation-mode` | `-md5-calculation-mode` | `MD5CalculationMode` | `RECALCULATE_IF_MISSING` | `no` | `optional` | 78 | `accepted` |
| `--output` | `-O` | `GATKPath` | `null` | `no` | `optional` | 75 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 94 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 68 | `accepted` |
| `--references-to-compare` | `-refcomp` | `List[GATKPath]` | `[]` | `yes` | `required` | 70 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--refcomp` | 70 |
| `--threads` | 88 |
| `--output-manifest` | 90 |
| `--manifest` | 90 |

## `fastgatk-check-reference-compatibility` ↔ `CheckReferenceCompatibility`

- Source: `fastgatk-native/src/check_reference_compatibility_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_reference_CheckReferenceCompatibility.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 50 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | 58 | `accepted` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | — | `ignored_via_catalog` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--output` | `-O` | `GATKPath` | `null` | `no` | `optional` | 66 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 72 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `no` | `optional` | — | `unsupported` |
| `--references-to-compare` | `-refcomp` | `List[GATKPath]` | `[]` | `yes` | `required` | 62 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--variant` | `-V` | `GATKPath` | `null` | `no` | `optional` | 60 | `accepted` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--refcomp` | 62 |
| `--output-manifest` | 68 |
| `--manifest` | 68 |

## `fastgatk-fasta-reference-maker` ↔ `FastaReferenceMaker`

- Source: `fastgatk-native/src/fasta_reference_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_fasta_FastaReferenceMaker.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 73 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | — | `ignored_via_catalog` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 97 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--line-width` | `` | `int` | `60` | `no` | `optional` | 107 | `accepted` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--output` | `-O` | `String` | `null` | `yes` | `required` | 105 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 131 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 94 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--interval` | 97 |
| `--region` | 98 |
| `--threads` | 109 |
| `--output-manifest` | 111 |
| `--manifest` | 111 |
| `-V` | 117 |
| `--variant` | 117 |
| `--variants` | 117 |
| `--snp-mask` | 124 |
| `--snp-mask-priority` | 126 |
| `--use-iupac-sample` | 128 |

## `fastgatk-fasta-alternate-reference-maker` ↔ `FastaAlternateReferenceMaker`

- Source: `fastgatk-native/src/fasta_reference_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_fasta_FastaAlternateReferenceMaker.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 73 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | — | `ignored_via_catalog` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 97 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--line-width` | `` | `int` | `60` | `no` | `optional` | 107 | `accepted` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--output` | `-O` | `String` | `null` | `yes` | `required` | 105 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 131 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 94 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--snp-mask` | `` | `FeatureInput[VariantContext]` | `null` | `no` | `optional` | 124 | `accepted` |
| `--snp-mask-priority` | `` | `boolean` | `false` | `no` | `optional` | 126 | `accepted` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-iupac-sample` | `` | `String` | `null` | `no` | `optional` | 128 | `accepted` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--variant` | `-V` | `FeatureInput[VariantContext]` | `null` | `yes` | `required` | 117 | `accepted` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--interval` | 97 |
| `--region` | 98 |
| `--threads` | 109 |
| `--output-manifest` | 111 |
| `--manifest` | 111 |
| `--variants` | 117 |

## `fastgatk-shift-fasta` ↔ `ShiftFasta`

- Source: `fastgatk-native/src/shift_fasta_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_fasta_ShiftFasta.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 62 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-file-name` | `` | `String` | `null` | `no` | `optional` | 80 | `accepted` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | — | `ignored_via_catalog` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--line-width` | `` | `int` | `60` | `no` | `optional` | 82 | `accepted` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--output` | `-O` | `String` | `null` | `yes` | `required` | 74 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 91 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 72 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--shift-back-output` | `` | `String` | `null` | `yes` | `required` | 76 | `accepted` |
| `--shift-offset-list` | `` | `List[Integer]` | `[]` | `no` | `optional` | 78 | `accepted` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--threads` | 84 |
| `--output-manifest` | 86 |
| `--manifest` | 86 |

## `fastgatk-index-feature-file` ↔ `IndexFeatureFile`

- Source: `fastgatk-native/src/index_feature_file_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_IndexFeatureFile.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 55 | `accepted` |
| `--input` | `-I` | `GATKPath` | `null` | `yes` | `required` | 65 | `accepted` |
| `--output` | `-O` | `GATKPath` | `null` | `no` | `optional` | 68 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 78 | `accepted` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--threads` | 71 |
| `--output-manifest` | 73 |
| `--manifest` | 73 |

## `fastgatk-count-reads` ↔ `CountReads`

- Source: `fastgatk-native/src/read_metrics_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_CountReads.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | 297 | `accepted` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | 349 | `accepted` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | 281 | `accepted` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 193 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `yes` | `required` | 226 | `accepted` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | 268 | `accepted` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | 249 | `accepted` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | 260 | `accepted` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | 239 | `accepted` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 232 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | 291 | `accepted` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--output` | `-O` | `GATKPath` | `null` | `no` | `optional` | 400 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 415 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | 285 | `accepted` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `no` | `optional` | 229 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--interval` | 232 |
| `--region` | 233 |
| `--read-name` | 303 |
| `--keep-read-group` | 305 |
| `--read-group-black-list` | 309 |
| `--read-filter-tag` | 319 |
| `--read-filter-tag-comp` | 322 |
| `--read-filter-tag-op` | 331 |
| `--min-read-length` | 354 |
| `--max-read-length` | 361 |
| `--minimum-mapping-quality` | 368 |
| `--maximum-mapping-quality` | 376 |
| `--min-fragment-length` | 384 |
| `--max-fragment-length` | 392 |
| `--batch-records` | 403 |
| `--threads` | 408 |
| `--output-manifest` | 410 |
| `--manifest` | 410 |

## `fastgatk-flag-stat` ↔ `FlagStat`

- Source: `fastgatk-native/src/read_metrics_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_FlagStat.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | 297 | `accepted` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | 349 | `accepted` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | 281 | `accepted` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 193 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `yes` | `required` | 226 | `accepted` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | 268 | `accepted` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | 249 | `accepted` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | 260 | `accepted` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | 239 | `accepted` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 232 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | 291 | `accepted` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--output` | `-O` | `GATKPath` | `null` | `no` | `optional` | 400 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 415 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | 285 | `accepted` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `no` | `optional` | 229 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--interval` | 232 |
| `--region` | 233 |
| `--read-name` | 303 |
| `--keep-read-group` | 305 |
| `--read-group-black-list` | 309 |
| `--read-filter-tag` | 319 |
| `--read-filter-tag-comp` | 322 |
| `--read-filter-tag-op` | 331 |
| `--min-read-length` | 354 |
| `--max-read-length` | 361 |
| `--minimum-mapping-quality` | 368 |
| `--maximum-mapping-quality` | 376 |
| `--min-fragment-length` | 384 |
| `--max-fragment-length` | 392 |
| `--batch-records` | 403 |
| `--threads` | 408 |
| `--output-manifest` | 410 |
| `--manifest` | 410 |

## `fastgatk-split-intervals` ↔ `SplitIntervals`

- Source: `fastgatk-native/src/split_intervals_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_SplitIntervals.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | 208 | `accepted` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--dont-mix-contigs` | `` | `boolean` | `false` | `no` | `optional` | 188 | `accepted` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | 169 | `accepted` |
| `--extension` | `` | `String` | `-scattered.interval_list` | `no` | `optional` | 190 | `accepted` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 146 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-file-num-digits` | `` | `int` | `4` | `no` | `optional` | 195 | `accepted` |
| `--interval-file-prefix` | `` | `String` | `""` | `no` | `optional` | 192 | `accepted` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | — | `ignored_via_catalog` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 166 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | 209 | `accepted` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--min-contig-size` | `` | `int` | `0` | `no` | `optional` | 199 | `accepted` |
| `--output` | `-O` | `File` | `null` | `yes` | `required` | 173 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 208 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `no` | `optional` | 163 | `accepted` |
| `--scatter-count` | `-scatter` | `int` | `1` | `no` | `optional` | 176 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--subdivision-mode` | `-mode` | `IntervalListScatterMode` | `INTERVAL_SUBDIVISION` | `no` | `optional` | 181 | `accepted` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | 212 | `accepted` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--scatter` | 176 |
| `--mode` | 181 |
| `--output-manifest` | 203 |
| `--manifest` | 203 |
| `--java-options` | 212 |

## `fastgatk-filter-intervals` ↔ `FilterIntervals`

- Source: `fastgatk-native/src/filter_intervals_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_copynumber_FilterIntervals.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--annotated-intervals` | `` | `File` | `null` | `no` | `optional` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | 140 | `accepted` |
| `--extreme-count-filter-maximum-percentile` | `` | `double` | `99.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--extreme-count-filter-minimum-percentile` | `` | `double` | `1.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--extreme-count-filter-percentage-of-samples` | `` | `double` | `90.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 113 | `accepted` |
| `--input` | `-I` | `List[File]` | `[]` | `no` | `optional` | 145 | `accepted` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | 179 | `accepted` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | 153 | `accepted` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | 169 | `accepted` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | 158 | `accepted` |
| `--intervals` | `-L` | `List[String]` | `[]` | `yes` | `required` | 134 | `accepted` |
| `--low-count-filter-count-threshold` | `` | `int` | `10` | `no` | `optional` | — | `ignored_via_catalog` |
| `--low-count-filter-percentage-of-samples` | `` | `double` | `50.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--maximum-gc-content` | `` | `double` | `0.9` | `no` | `optional` | — | `ignored_via_catalog` |
| `--maximum-mappability` | `` | `double` | `1.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--maximum-segmental-duplication-content` | `` | `double` | `0.5` | `no` | `optional` | — | `ignored_via_catalog` |
| `--minimum-gc-content` | `` | `double` | `0.1` | `no` | `optional` | — | `ignored_via_catalog` |
| `--minimum-mappability` | `` | `double` | `0.9` | `no` | `optional` | — | `ignored_via_catalog` |
| `--minimum-segmental-duplication-content` | `` | `double` | `0.0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--output` | `-O` | `File` | `null` | `yes` | `required` | 148 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 217 | `accepted` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

## `fastgatk-preprocess-intervals` ↔ `PreprocessIntervals`

- Source: `fastgatk-native/src/preprocess_intervals_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_copynumber_PreprocessIntervals.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--bin-length` | `` | `int` | `1000` | `no` | `optional` | — | `ignored_via_catalog` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | 99 | `accepted` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 79 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `no` | `common` | — | `unsupported` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | — | `ignored_via_catalog` |
| `--intervals` | `-L` | `List[String]` | `[]` | `no` | `optional` | 93 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--output` | `-O` | `File` | `null` | `yes` | `required` | 102 | `accepted` |
| `--padding` | `` | `int` | `250` | `no` | `optional` | — | `ignored_via_catalog` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 135 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 91 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

## `fastgatk-collect-read-counts` ↔ `CollectReadCounts`

- Source: `fastgatk-native/src/collect_read_counts_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_copynumber_CollectReadCounts.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | 179 | `accepted` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | 192 | `accepted` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | 193 | `accepted` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | 119 | `accepted` |
| `--format` | `` | `Format` | `HDF5` | `no` | `optional` | 124 | `accepted` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 92 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `yes` | `required` | 108 | `accepted` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | 151 | `accepted` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | 141 | `accepted` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | 150 | `accepted` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | — | `ignored_via_catalog` |
| `--intervals` | `-L` | `List[String]` | `[]` | `yes` | `required` | 112 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--output` | `-O` | `File` | `null` | `yes` | `required` | 122 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 164 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | 168 | `accepted` |
| `--reference` | `-R` | `GATKPath` | `null` | `no` | `optional` | 110 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | 198 | `accepted` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | 183 | `accepted` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | 185 | `accepted` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | 186 | `accepted` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | 197 | `accepted` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--interval` | 113 |
| `--region` | 113 |
| `--sample` | 126 |
| `--minimum-mapping-quality` | 128 |
| `--batch-records` | 131 |
| `--threads` | 134 |
| `--output-manifest` | 136 |
| `--manifest` | 136 |
| `--include-duplicates` | 139 |
| `--java-options` | 196 |

## `fastgatk-collect-f1r2-counts` ↔ `CollectF1R2Counts`

- Source: `fastgatk-native/src/collect_f1r2_counts_tool.cpp`
- GATK reference: `<none>` (native-only bridge binary)
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

No gatkdoc JSON is shipped with this GATK version; the tool is native-only.

## `fastgatk-collect-allelic-counts` ↔ `CollectAllelicCounts`

- Source: `fastgatk-native/src/collect_allelic_counts_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_copynumber_CollectAllelicCounts.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | 260 | `accepted` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | 275 | `accepted` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | 254 | `accepted` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | 208 | `accepted` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 177 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `yes` | `required` | 197 | `accepted` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | 234 | `accepted` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | 229 | `accepted` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | 233 | `accepted` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | 243 | `accepted` |
| `--intervals` | `-L` | `List[String]` | `[]` | `yes` | `required` | 201 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-depth-per-sample` | `-max-depth-per-sample` | `int` | `0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--minimum-base-quality` | `` | `int` | `20` | `no` | `optional` | — | `ignored_via_catalog` |
| `--output` | `-O` | `File` | `null` | `yes` | `required` | 211 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 275 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | 270 | `accepted` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 199 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--include-duplicates` | 251 |

## `fastgatk-depth-of-coverage` ↔ `DepthOfCoverage`

- Source: `fastgatk-native/src/depth_of_coverage_tool.cpp`
- GATK arguments live in the JSON at the canonical `gatkdoc/org_broadinstitute_hellbender_tools_walkers_coverage_DepthOfCoverage.json` location; this generator reads `arguments[]` directly.
- Status column meanings:
    - `accepted` — native argv parser handles the option literal.
    - `ignored_via_catalog` — option appears in `fastgatk::cli::kCatalog` and is eaten by `fastgatk::cli::consume(...)`.
    - `unsupported` — GATK exposes the option but native does not handle it.

| GATK option | Synonyms | Type | Default | Required | Kind | Native parsing line(s) | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `--add-output-sam-program-record` | `-add-output-sam-program-record` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--add-output-vcf-command-line` | `-add-output-vcf-command-line` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--arguments_file` | `` | `List[File]` | `[]` | `no` | `optional` | — | `unsupported` |
| `--calculate-coverage-over-genes` | `-gene-list` | `List[String]` | `[]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--cloud-index-prefetch-buffer` | `-CIPB` | `int` | `-1` | `no` | `optional` | — | `unsupported` |
| `--cloud-prefetch-buffer` | `-CPB` | `int` | `40` | `no` | `optional` | — | `unsupported` |
| `--count-type` | `` | `CountPileupType` | `COUNT_READS` | `no` | `optional` | — | `ignored_via_catalog` |
| `--create-output-bam-index` | `-OBI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-bam-md5` | `-OBM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-index` | `-OVI` | `boolean` | `true` | `no` | `common` | — | `ignored_via_catalog` |
| `--create-output-variant-md5` | `-OVM` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--disable-bam-index-caching` | `-DBIC` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--disable-read-filter` | `-DF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--disable-sequence-dictionary-validation` | `-disable-sequence-dictionary-validation` | `boolean` | `false` | `no` | `optional` | 264 | `accepted` |
| `--disable-tool-default-read-filters` | `-disable-tool-default-read-filters` | `boolean` | `false` | `no` | `common` | — | `unsupported` |
| `--exclude-intervals` | `-XL` | `List[String]` | `[]` | `no` | `common` | 173 | `accepted` |
| `--gatk-config-file` | `` | `String` | `null` | `no` | `common` | — | `unsupported` |
| `--gcs-max-retries` | `-gcs-retries` | `int` | `20` | `no` | `optional` | — | `unsupported` |
| `--gcs-project-for-requester-pays` | `` | `String` | `""` | `no` | `optional` | — | `unsupported` |
| `--help` | `-h` | `boolean` | `false` | `no` | `optional` | 134 | `accepted` |
| `--ignore-deletion-sites` | `` | `boolean` | `false` | `no` | `advanced` | 227 | `accepted` |
| `--include-deletions` | `` | `boolean` | `false` | `no` | `advanced` | 225 | `accepted` |
| `--include-ref-n-sites` | `` | `boolean` | `false` | `no` | `advanced` | 247 | `accepted` |
| `--input` | `-I` | `List[GATKPath]` | `[]` | `yes` | `required` | 164 | `accepted` |
| `--interval-exclusion-padding` | `-ixp` | `int` | `0` | `no` | `common` | 190 | `accepted` |
| `--interval-merging-rule` | `-imr` | `IntervalMergingRule` | `ALL` | `no` | `optional` | — | `ignored_via_catalog` |
| `--interval-padding` | `-ip` | `int` | `0` | `no` | `common` | 185 | `accepted` |
| `--interval-set-rule` | `-isr` | `IntervalSetRule` | `UNION` | `no` | `common` | 178 | `accepted` |
| `--intervals` | `-L` | `List[String]` | `[]` | `yes` | `required` | 168 | `accepted` |
| `--inverted-read-filter` | `-XRF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--lenient` | `-LE` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--max-base-quality` | `` | `byte` | `127` | `no` | `optional` | — | `ignored_via_catalog` |
| `--max-depth-per-sample` | `-max-depth-per-sample` | `int` | `0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--max-variants-per-shard` | `` | `int` | `0` | `no` | `common` | — | `ignored_via_catalog` |
| `--min-base-quality` | `` | `byte` | `0` | `no` | `optional` | — | `ignored_via_catalog` |
| `--nBins` | `` | `int` | `499` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--omit-depth-output-at-each-base` | `` | `boolean` | `false` | `no` | `optional` | 235 | `accepted` |
| `--omit-genes-not-entirely-covered-by-traversal` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--omit-interval-statistics` | `` | `boolean` | `false` | `no` | `optional` | 239 | `accepted` |
| `--omit-locus-table` | `` | `boolean` | `false` | `no` | `optional` | 231 | `accepted` |
| `--omit-per-sample-statistics` | `` | `boolean` | `false` | `no` | `optional` | 243 | `accepted` |
| `--output` | `-O` | `String` | `null` | `yes` | `required` | 197 | `accepted` |
| `--output-format` | `` | `DEPTH_OF_COVERAGE_OUTPUT_FORMAT` | `CSV` | `no` | `optional` | — | `ignored_via_catalog` |
| `--partition-type` | `-pt` | `EnumSet[Partition]` | `[sample]` | `no` | `optional` | — | `ignored_via_catalog` |
| `--print-base-counts` | `` | `boolean` | `false` | `no` | `optional` | 223 | `accepted` |
| `--QUIET` | `` | `Boolean` | `false` | `no` | `common` | 264 | `accepted` |
| `--read-filter` | `-RF` | `List[String]` | `[]` | `no` | `common` | — | `unsupported` |
| `--read-index` | `-read-index` | `List[GATKPath]` | `[]` | `no` | `common` | — | `ignored_via_catalog` |
| `--read-validation-stringency` | `-VS` | `ValidationStringency` | `SILENT` | `no` | `common` | — | `ignored_via_catalog` |
| `--reference` | `-R` | `GATKPath` | `null` | `yes` | `required` | 166 | `accepted` |
| `--seconds-between-progress-updates` | `-seconds-between-progress-updates` | `double` | `10.0` | `no` | `common` | — | `ignored_via_catalog` |
| `--sequence-dictionary` | `-sequence-dictionary` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--showHidden` | `-showHidden` | `boolean` | `false` | `no` | `advanced` | — | `unsupported` |
| `--sites-only-vcf-output` | `` | `boolean` | `false` | `no` | `optional` | — | `ignored_via_catalog` |
| `--start` | `` | `int` | `1` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--stop` | `` | `int` | `500` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--summary-coverage-threshold` | `` | `List[Integer]` | `[15]` | `no` | `advanced` | — | `ignored_via_catalog` |
| `--tmp-dir` | `` | `GATKPath` | `null` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-deflater` | `-jdk-deflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--use-jdk-inflater` | `-jdk-inflater` | `boolean` | `false` | `no` | `common` | — | `ignored_via_catalog` |
| `--verbosity` | `-verbosity` | `LogLevel` | `INFO` | `no` | `common` | — | `ignored_via_catalog` |
| `--version` | `` | `boolean` | `false` | `no` | `optional` | — | `unsupported` |

### Native-only options (no GATK counterpart in 4.6.2.0)

| Native literal | Native parsing line(s) |
| --- | --- |
| `--omit-genes` | 249 |
| `--omit-intervals` | 253 |
| `--omit-sample-summary` | 259 |
