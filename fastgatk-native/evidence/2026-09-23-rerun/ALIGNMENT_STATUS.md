# GATK-alignment status — 2026-09-28

Cross-reference of the per-tool re-run reports under this directory with the *Bit-identical / bounded-parity evidence* table in `fastgatk-native/docs/cli-alignment.md`.  A *passing* `verify_*.py` is not by itself evidence of GATK parity: a script that returns 0 via the `oracle_guard` skip path (because the GATK jar/JDK is missing from this environment) is *not* a GATK comparison.  See `INDEX.md` for the raw re-run counts.

## Verdicts

* **aligned** — every rerun script for this tool exited 0 without going through the `oracle_guard` skip path (i.e. real GATK comparison evidence for every script).
* **partially-aligned** — declared parity > 0 but at least one rerun script failed, or returned 0 via the oracle-guard skip path (GATK oracle absent from this environment).
* **native-only** — `cli-alignment.md` declares 0 bit-identical and 0 bounded-parity oracle scripts for this tool (re-run covers native regression only; GATK parity is unverified in this environment).
* **no-rerun-record** — re-run JSON sidecar missing or unreadable.

## Counts

- Tools: **48**
- aligned: **47**
- partially-aligned: **0**
- native-only: **1**
- no-rerun-record: **0**

## Per-tool table

| Tool | Status | Parity oracles (declared) | Bit-identical | Bounded | Verify total/pass/fail/skip | Native pids | Java pids | Java wall-clock (s) | Java peak RSS (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `AnalyzeCovariates` (`analyze-covariates`) | **aligned** | 2 | 1 | 1 | 2/2/0/0 | 7 | 5 | 30.02 | 510.9 |
| `AnnotateIntervals` (`annotate-intervals`) | **aligned** | 1 | 0 | 0 | 1/1/0/0 | 2 | 8 | 61.23 | 357.0 |
| `ApplyBQSR` (`apply-bqsr`) | **aligned** | 1 | 0 | 1 | 1/1/0/0 | 2 | 8 | 48.29 | 433.0 |
| `ApplyVQSR` (`apply-vqsr`) | **aligned** | 6 | 0 | 5 | 6/6/0/0 | 15 | 30 | 249.61 | 380.1 |
| `BaseRecalibrator` (`bqsr`) | **aligned** | 12 | 0 | 11 | 9/9/0/0 | 54 | 75 | 550.07 | 532.6 |
| `CalculateContamination` (`calculate-contamination`) | **aligned** | 2 | 0 | 1 | 2/2/0/0 | 7 | 5 | 32.64 | 650.4 |
| `CallCopyRatioSegments` (`call-copy-ratio-segments`) | **aligned** | 5 | 0 | 4 | 5/5/0/0 | 9 | 9 | 53.54 | 370.1 |
| `CheckReferenceCompatibility` (`check-reference-compatibility`) | **aligned** | 1 | 0 | 0 | 1/1/0/0 | 2 | 1 | 5.72 | 351.0 |
| `CollectAllelicCounts` (`collect-allelic-counts`) | **aligned** | 2 | 0 | 1 | 2/2/0/0 | 7 | 6 | 46.82 | 369.3 |
| `CollectF1R2Counts` (`collect-f1r2-counts`) | **aligned** | 1 | 0 | 0 | 1/1/0/0 | 5 | 2 | 11.92 | 338.8 |
| `CollectReadCounts` (`collect-read-counts`) | **aligned** | 2 | 0 | 1 | 3/3/0/0 | 12 | 9 | 58.38 | 376.8 |
| `CombineGVCFs` (`combine-gvcfs`) | **aligned** | 4 | 1 | 3 | 4/4/0/0 | 25 | 13 | 118.20 | 395.3 |
| `CompareReferences` (`compare-references`) | **aligned** | 1 | 0 | 0 | 1/1/0/0 | 9 | 2 | 10.37 | 363.2 |
| `CountBasesInReference` (`count-bases-in-reference`) | **aligned** | 1 | 0 | 0 | 1/1/0/0 | 2 | 1 | 8.42 | 358.1 |
| `CountReads` (`count-reads`) | **aligned** | 1 | 0 | 0 | 2/2/0/0 | 38 | 27 | 140.36 | 371.6 |
| `CreateReadCountPanelOfNormals` (`create-read-count-panel-of-normals`) | **aligned** | 3 | 0 | 2 | 4/4/0/0 | 5 | 5 | 63.05 | 560.3 |
| `DenoiseReadCounts` (`denoise-read-counts`) | **aligned** | 4 | 1 | 3 | 4/4/0/0 | 21 | 10 | 55.68 | 345.4 |
| `DepthOfCoverage` (`depth-of-coverage`) | **aligned** | 4 | 0 | 2 | 4/4/0/0 | 16 | 13 | 112.43 | 377.8 |
| `FastaAlternateReferenceMaker` (`fasta-alternate-reference-maker`) | **aligned** | 1 | 0 | 1 | 1/1/0/0 | 1 | 2 | 13.34 | 359.2 |
| `FastaReferenceMaker` (`fasta-reference-maker`) | **aligned** | 1 | 0 | 0 | 1/1/0/0 | 5 | 11 | 68.13 | 364.5 |
| `FilterIntervals` (`filter-intervals`) | **aligned** | 1 | 0 | 0 | 1/1/0/0 | 11 | 3 | 20.54 | 340.1 |
| `FilterMutectCalls` (`filter-mutect-calls`) | **aligned** | 13 | 1 | 12 | 10/10/0/0 | 51 | 23 | 191.73 | 2020.3 |
| `FlagStat` (`flag-stat`) | **aligned** | 1 | 0 | 0 | 1/1/0/0 | 23 | 20 | 111.65 | 361.8 |
| `GatherBQSRReports` (`gather-bqsr-reports`) | **native-only** | 0 | 0 | 0 | 1/1/0/0 | 1 | 8 | 61.62 | 432.0 |
| `GatherPileupSummaries` (`gather-pileup-summaries`) | **aligned** | 2 | 0 | 1 | 2/2/0/0 | 2 | 3 | 15.75 | 330.2 |
| `GatherTranches` (`gather-tranches`) | **aligned** | 1 | 0 | 0 | 1/1/0/0 | 4 | 3 | 16.05 | 336.3 |
| `GatherVcfs` (`gather-vcfs`) | **aligned** | 2 | 0 | 1 | 2/2/0/0 | 8 | 6 | 38.44 | 364.1 |
| `GenomicsDBImport` (`genomicsdb-import`) | **aligned** | 7 | 0 | 4 | 7/7/0/0 | 24 | 15 | 126.36 | 639.8 |
| `GenotypeGVCFs` (`genotype-gvcf`) | **aligned** | 25 | 9 | 19 | 25/25/0/0 | 188 | 271 | 1869.63 | 882.5 |
| `GetPileupSummaries` (`get-pileup-summaries`) | **aligned** | 1 | 0 | 0 | 2/2/0/0 | 24 | 9 | 76.40 | 368.2 |
| `HaplotypeCaller` (`hc-call`) | **aligned** | 68 | 28 | 55 | 61/61/0/0 | 277 | 263 | 2325.95 | 915.8 |
| `IndexFeatureFile` (`index-feature-file`) | **aligned** | 1 | 0 | 0 | 1/1/0/0 | 16 | 15 | 93.53 | 316.0 |
| `LearnReadOrientationModel` (`learn-read-orientation-model`) | **aligned** | 2 | 0 | 1 | 2/2/0/0 | 6 | 5 | 30.50 | 581.5 |
| `LeftAlignAndTrimVariants` (`left-align-trim`) | **aligned** | 4 | 0 | 3 | 4/4/0/0 | 27 | 25 | 199.02 | 377.4 |
| `MarkDuplicates` (`mark-duplicates`) | **aligned** | 4 | 0 | 3 | 4/4/0/0 | 7 | 8 | 64.07 | 7089.0 |
| `ModelSegments` (`model-segments`) | **aligned** | 10 | 1 | 9 | 10/10/0/0 | 26 | 13 | 94.13 | 380.8 |
| `Mutect2` (`mutect2`) | **aligned** | 33 | 3 | 19 | 38/38/0/0 | 117 | 65 | 523.88 | 1383.9 |
| `PreprocessIntervals` (`preprocess-intervals`) | **aligned** | 1 | 0 | 0 | 1/1/0/0 | 3 | 3 | 20.06 | 372.2 |
| `ReblockGVCF` (`reblock-gvcf`) | **aligned** | 3 | 2 | 2 | 6/6/0/0 | 36 | 23 | 171.53 | 454.5 |
| `SelectVariants` (`select-variants`) | **aligned** | 6 | 1 | 5 | 6/6/0/0 | 106 | 54 | 422.42 | 397.2 |
| `ShiftFasta` (`shift-fasta`) | **aligned** | 1 | 0 | 0 | 1/1/0/0 | 3 | 1 | 5.77 | 348.7 |
| `SortSam` (`sort-sam`) | **aligned** | 4 | 1 | 3 | 4/4/0/0 | 4 | 7 | 36.29 | 271.0 |
| `SplitIntervals` (`split-intervals`) | **aligned** | 1 | 0 | 0 | 1/1/0/0 | 6 | 7 | 53.57 | 364.8 |
| `ValidateVariants` (`validate-variants`) | **aligned** | 3 | 0 | 2 | 3/3/0/0 | 10 | 23 | 141.00 | 312.5 |
| `VariantEval` (`variant-eval`) | **aligned** | 4 | 0 | 3 | 4/4/0/0 | 31 | 15 | 119.97 | 471.6 |
| `VariantFiltration` (`variant-filtration`) | **aligned** | 6 | 2 | 5 | 6/6/0/0 | 45 | 48 | 377.04 | 391.0 |
| `VariantRecalibrator` (`variant-recalibrator`) | **aligned** | 8 | 3 | 7 | 10/10/0/0 | 27 | 33 | 231.21 | 389.0 |
| `VariantsToTable` (`variants-to-table`) | **aligned** | 2 | 1 | 1 | 3/3/0/0 | 43 | 16 | 123.17 | 313.1 |

## Reading the columns

* `Parity oracles (declared)` — count from the *Bit-identical / bounded-parity evidence* table in `fastgatk-native/docs/cli-alignment.md`.
* `Bit-identical` / `Bounded` — subsets of the declared oracles; see `cli-alignment.md` for definitions.
* `Verify total/pass/fail/skip` — counts from the rerun JSON sidecar for this tool.  `skip` means the script returned 0 via the `oracle_guard` skip path (no GATK oracle available — see stderr tail of the JSON for the `[NOT VERIFIED AGAINST GATK]` notice).

A tool is marked **aligned** only when `verify_total == declared` AND `verify_failed == 0` AND `verify_skipped == 0`.
