# GATK-alignment status — 2026-09-28

Cross-reference of the per-tool re-run reports under this directory with the *Bit-identical / bounded-parity evidence* table in `fastgatk-native/docs/cli-alignment.md`.  A *passing* `verify_*.py` is not by itself evidence of GATK parity: a script that returns 0 via the `oracle_guard` skip path (because the GATK jar/JDK is missing from this environment) is *not* a GATK comparison.  See `INDEX.md` for the raw re-run counts.

## Verdicts

* **aligned** — every rerun script for this tool exited 0 without going through the `oracle_guard` skip path (i.e. real GATK comparison evidence for every script).
* **partially-aligned** — declared parity > 0 but at least one rerun script failed, or returned 0 via the oracle-guard skip path (GATK oracle absent from this environment).
* **native-only** — `cli-alignment.md` declares 0 bit-identical and 0 bounded-parity oracle scripts for this tool (re-run covers native regression only; GATK parity is unverified in this environment).
* **no-rerun-record** — re-run JSON sidecar missing or unreadable.

## Counts

- Tools: **48**
- aligned: **0**
- partially-aligned: **0**
- native-only: **0**
- no-rerun-record: **48**

## Per-tool table

| Tool | Status | Parity oracles (declared) | Bit-identical | Bounded | Verify total/pass/fail/skip | Native pids | Java pids | Java wall-clock (s) | Java peak RSS (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `AnalyzeCovariates` (`analyze-covariates`) | no-rerun-record | 2 | 1 | 1 | — | — | — | — | — |
| `AnnotateIntervals` (`annotate-intervals`) | no-rerun-record | 1 | 0 | 0 | — | — | — | — | — |
| `ApplyBQSR` (`apply-bqsr`) | no-rerun-record | 1 | 0 | 1 | — | — | — | — | — |
| `ApplyVQSR` (`apply-vqsr`) | no-rerun-record | 6 | 0 | 5 | — | — | — | — | — |
| `BaseRecalibrator` (`bqsr`) | no-rerun-record | 12 | 0 | 11 | — | — | — | — | — |
| `CalculateContamination` (`calculate-contamination`) | no-rerun-record | 2 | 0 | 1 | — | — | — | — | — |
| `CallCopyRatioSegments` (`call-copy-ratio-segments`) | no-rerun-record | 5 | 0 | 4 | — | — | — | — | — |
| `CheckReferenceCompatibility` (`check-reference-compatibility`) | no-rerun-record | 1 | 0 | 0 | — | — | — | — | — |
| `CollectAllelicCounts` (`collect-allelic-counts`) | no-rerun-record | 2 | 0 | 1 | — | — | — | — | — |
| `CollectF1R2Counts` (`collect-f1r2-counts`) | no-rerun-record | 1 | 0 | 0 | — | — | — | — | — |
| `CollectReadCounts` (`collect-read-counts`) | no-rerun-record | 2 | 0 | 1 | — | — | — | — | — |
| `CombineGVCFs` (`combine-gvcfs`) | no-rerun-record | 4 | 1 | 3 | — | — | — | — | — |
| `CompareReferences` (`compare-references`) | no-rerun-record | 1 | 0 | 0 | — | — | — | — | — |
| `CountBasesInReference` (`count-bases-in-reference`) | no-rerun-record | 1 | 0 | 0 | — | — | — | — | — |
| `CountReads` (`count-reads`) | no-rerun-record | 1 | 0 | 0 | — | — | — | — | — |
| `CreateReadCountPanelOfNormals` (`create-read-count-panel-of-normals`) | no-rerun-record | 3 | 0 | 2 | — | — | — | — | — |
| `DenoiseReadCounts` (`denoise-read-counts`) | no-rerun-record | 4 | 1 | 3 | — | — | — | — | — |
| `DepthOfCoverage` (`depth-of-coverage`) | no-rerun-record | 4 | 0 | 2 | — | — | — | — | — |
| `FastaAlternateReferenceMaker` (`fasta-alternate-reference-maker`) | no-rerun-record | 1 | 0 | 1 | — | — | — | — | — |
| `FastaReferenceMaker` (`fasta-reference-maker`) | no-rerun-record | 1 | 0 | 0 | — | — | — | — | — |
| `FilterIntervals` (`filter-intervals`) | no-rerun-record | 1 | 0 | 0 | — | — | — | — | — |
| `FilterMutectCalls` (`filter-mutect-calls`) | no-rerun-record | 13 | 1 | 12 | — | — | — | — | — |
| `FlagStat` (`flag-stat`) | no-rerun-record | 1 | 0 | 0 | — | — | — | — | — |
| `GatherBQSRReports` (`gather-bqsr-reports`) | no-rerun-record | 0 | 0 | 0 | — | — | — | — | — |
| `GatherPileupSummaries` (`gather-pileup-summaries`) | no-rerun-record | 2 | 0 | 1 | — | — | — | — | — |
| `GatherTranches` (`gather-tranches`) | no-rerun-record | 1 | 0 | 0 | — | — | — | — | — |
| `GatherVcfs` (`gather-vcfs`) | no-rerun-record | 2 | 0 | 1 | — | — | — | — | — |
| `GenomicsDBImport` (`genomicsdb-import`) | no-rerun-record | 7 | 0 | 4 | — | — | — | — | — |
| `GenotypeGVCFs` (`genotype-gvcf`) | no-rerun-record | 25 | 9 | 19 | — | — | — | — | — |
| `GetPileupSummaries` (`get-pileup-summaries`) | no-rerun-record | 1 | 0 | 0 | — | — | — | — | — |
| `HaplotypeCaller` (`hc-call`) | no-rerun-record | 68 | 28 | 55 | — | — | — | — | — |
| `IndexFeatureFile` (`index-feature-file`) | no-rerun-record | 1 | 0 | 0 | — | — | — | — | — |
| `LearnReadOrientationModel` (`learn-read-orientation-model`) | no-rerun-record | 2 | 0 | 1 | — | — | — | — | — |
| `LeftAlignAndTrimVariants` (`left-align-trim`) | no-rerun-record | 4 | 0 | 3 | — | — | — | — | — |
| `MarkDuplicates` (`mark-duplicates`) | no-rerun-record | 4 | 0 | 3 | — | — | — | — | — |
| `ModelSegments` (`model-segments`) | no-rerun-record | 10 | 1 | 9 | — | — | — | — | — |
| `Mutect2` (`mutect2`) | no-rerun-record | 33 | 3 | 19 | — | — | — | — | — |
| `PreprocessIntervals` (`preprocess-intervals`) | no-rerun-record | 1 | 0 | 0 | — | — | — | — | — |
| `ReblockGVCF` (`reblock-gvcf`) | no-rerun-record | 3 | 2 | 2 | — | — | — | — | — |
| `SelectVariants` (`select-variants`) | no-rerun-record | 6 | 1 | 5 | — | — | — | — | — |
| `ShiftFasta` (`shift-fasta`) | no-rerun-record | 1 | 0 | 0 | — | — | — | — | — |
| `SortSam` (`sort-sam`) | no-rerun-record | 4 | 1 | 3 | — | — | — | — | — |
| `SplitIntervals` (`split-intervals`) | no-rerun-record | 1 | 0 | 0 | — | — | — | — | — |
| `ValidateVariants` (`validate-variants`) | no-rerun-record | 3 | 0 | 2 | — | — | — | — | — |
| `VariantEval` (`variant-eval`) | no-rerun-record | 4 | 0 | 3 | — | — | — | — | — |
| `VariantFiltration` (`variant-filtration`) | no-rerun-record | 6 | 2 | 5 | — | — | — | — | — |
| `VariantRecalibrator` (`variant-recalibrator`) | no-rerun-record | 8 | 3 | 7 | — | — | — | — | — |
| `VariantsToTable` (`variants-to-table`) | no-rerun-record | 2 | 1 | 1 | — | — | — | — | — |

## Reading the columns

* `Parity oracles (declared)` — count from the *Bit-identical / bounded-parity evidence* table in `fastgatk-native/docs/cli-alignment.md`.
* `Bit-identical` / `Bounded` — subsets of the declared oracles; see `cli-alignment.md` for definitions.
* `Verify total/pass/fail/skip` — counts from the rerun JSON sidecar for this tool.  `skip` means the script returned 0 via the `oracle_guard` skip path (no GATK oracle available — see stderr tail of the JSON for the `[NOT VERIFIED AGAINST GATK]` notice).

A tool is marked **aligned** only when `verify_total == declared` AND `verify_failed == 0` AND `verify_skipped == 0`.
