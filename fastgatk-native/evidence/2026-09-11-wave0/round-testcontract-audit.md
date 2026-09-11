# Round audit: native-internal registered tests that could block a GATK-parity fix

Scope: every test registered by `fastgatk-native/CMakeLists.txt` (`add_test`, 273 total) that is
**not** a GATK oracle — i.e. tests that assert **native-internal** expectations and would go red if
native were changed to match GATK, even though the assertion was never validated against GATK.

Read-only audit. No file in the repo was modified except this report. All scratch work under
`tempfile.TemporaryDirectory()`.

Toolchain used for confirmations:
`third_party/jdk17/bin/java -Xmx1g -jar third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar`
against `fastgatk-native/build/*`. The pinned GATK **source tree** (`gatk-source/`) was used to
classify semantics without running anything.

---

## 0. Executive summary

* 273 `add_test` entries; **94** do not carry `oracle` in the test name. Of those 94, **~30 are
  oracle-style by construction** (they byte-compare against the pinned jar or against checked-in
  GATK fixture files) and are therefore out of scope; the remaining ~64 are genuine
  native-internal contracts.
* **Six NEW class-(c) divergences are pinned as contracts, all measured against the pinned jar**
  (three in the gVCF subsystem that already produced the `verify_indel.py:173` incident, three in
  variant manipulation):
  1. `verify_reblock_gvcf.py:64,82-83` pins that native *emits* a reference block for a
     `--drop-low-quals` site; pinned GATK 4.6.2.0 **drops the record entirely** (§3.2).
  2. `verify_genotype_gvcf.py:391` pins that native emits a `*,G` spanning-deletion record;
     pinned GATK 4.6.2.0 emits **0 records** by default and, with `-all-sites`, emits a record
     with `ALT=.` / `GT=./.` (§3.3).
  3. `verify_reblock_gvcf.py:283-287` pins native's triploid "keep the call + compact PLs";
     pinned GATK 4.6.2.0 **reblocks the site to a hom-ref block** (§3.2).
  4. `verify_select_variants.py:820` `assert ref_only_records == []` — GATK **keeps** the
     all-hom-ref record as `ALT=.` (§4.4.1).
  5. `verify_variant_recalibrator.py:84,422,437,456` `culprit=full-covariance-gmm` — the string
     does not exist anywhere in the GATK jar; GATK writes `culprit=MQ` (§4.4.2).
  6. `verify_variant_filtration.py:504-505` `AS_FilterStatus=LowASQD,PASS` — GATK emits
     `AS_FilterStatus=SITE|SITE` and applies no allele filter (§4.4.3).
  A further three measured class-(c) items were found in somatic/coverage tools (MBQ/MMQ arity,
  `--min-tlod`/`low_tlod`/`FAIL`, DepthOfCoverage per-sample columns) plus four input-envelope
  divergences (§4.1).
* A second, subtler failure mode is pervasive: many `-contract` tests invoke the pinned jar
  **only as a reader** (`SelectVariants`) or behind `if java.is_file() and jar.is_file():` with no
  assertion that the jar exists. They *look* like oracles but their expected values are
  native's own recorded output. `verify_reblock_gvcf.py`, `verify_genotype_gvcf.py`,
  `verify_combine_gvcfs.py`, `verify_select_variants.py`, `verify_gather_vcfs.py`,
  `verify_variant_recalibrator.py` are all in this class; 13 of the 18 variant-manipulation
  scripts are dual-mode, and the `if java_oracle:` pattern is used in at least 11 scripts (§4.3).
* Highest-value cheap fix: every item above was settled with a single short GATK run on the test's
  own fixture (§3.1–§3.6, §4.4). They should be converted into real oracles or their expectations
  corrected.

---

## 1. Classification key

* **(a)** genuine GATK contract — the asserted value is produced by GATK itself (confirmed by
  running the pinned jar, or read out of `gatk-source/`).
* **(b)** native-internal implementation detail that happens to match GATK today, or that has no
  GATK counterpart (native-only flags/manifest/telemetry). Fixing a GATK divergence *can* break it.
* **(c)** KNOWN/likely to actively lock in a divergence from GATK — a parity fix would turn the
  test red, and the test is currently the only thing asserting the divergent behaviour.

An assertion is listed as "record-level" if it constrains VCF/gVCF record content: columns
2–10 (POS/ID/REF/ALT/QUAL/FILTER/INFO/FORMAT/sample), record count, or record order.

---

## 2. Table — every audited non-oracle-named registered test

CMake line = line of `add_test(NAME ...)` in `fastgatk-native/CMakeLists.txt`.

### 2.1 Scaffolding / build / infra (no VCF records) — class (b) or n/a

| line | test name | script / target | record-level? | class | notes |
|---|---|---|---|---|---|
| 1030 | `fastgatk-launcher-contract` | `verify_launcher_contract.py` | no (argv/exit codes) | (b) | pins native `expect=2` + `UNSUPPORTED_PARAMETER`; GATK uses exit 1 + its own text (§4.1) |
| 1033 | `fastgatk-progress-score-contract` | `compute_progress_score.py` | no | (b) | score arithmetic |
| 1036 | `fastgatk-fixture-digest-contract` | `verify_fixture_digests.py` | no (input hashes) | n/a | hashes inputs only — cannot detect an unverified expectation |
| 1039 | `fastgatk-kokkos-api-boundary` | `verify_kokkos_api_boundary.py` | no (source regex) | n/a | |
| 1047 | `fastgatk-kokkos-backend-config` | `verify_kokkos_backend_config.py` | no (CMake cache) | n/a | |
| 1101 | `fastgatk-remote-staging-contract` | `verify_remote_staging.py` | TSV row text | (b) | `:127,:157` `["chr1\t1\t8\t0.500000"]` — genuine GATK output but **no GATK call in this file** |
| 1104 | `fastgatk-resource-limits-contract` | `verify_resource_limits.py` | no (JSON summary) | (b) | `:42 effective_threads==2`, `:44 host_hard_bytes==1048576` |
| 1109 | `fastgatk-slurm-resource-wrapper-contract` | `verify_slurm_wrapper.py` | no (env/argv) | (b) | |
| 1112 | `fastgatk-scheduler-retry-contract` | `verify_scheduler_retry_contract.py` | no (file names) | (b) | fake shards write header-only VCFs |
| 1119 | `fastgatk-workflow-local-contract` | `verify_pipeline_local.sh` | header line only | (b) | `:70-71` `##fileformat=VCFv4.2`, `#CHROM...` |
| 2496 | `fastgatk-genomicsdb-import-benchmark` | `benchmark_genomicsdb_import.py` | no | n/a | benchmark |
| 1785 | `fastgatk-reblock-gvcf-java-benchmark-command-contract` | `benchmark_reblock_gvcf_java.py --dry-run` | no | n/a | |
| 286 | `fastgatk-hts-reader-indel-quality` | `fastgatk-hts-reader-test` (C++) | no | n/a | SAM reader unit test |
| 295 | `fastgatk-hts-reader-unbounded-index-capability` | `fastgatk-hts-reader-test` (C++) | no | n/a | |
| 350 | `fastgatk-calling-reference-boundary-smoke` | C++ | no | n/a | asserts activity_loci/region spans, not VCF |
| 362 | `fastgatk-haplotype-cigar-trim-smoke` | C++ | no | n/a | asserts `sequence_begin/end` |
| 394 | `fastgatk-calling-indel-smoke` | C++ | no | n/a | asserts prior vectors |
| 750 | `fastgatk-contamination-kernel-segmenter-contract` | C++ | no | n/a | asserts changepoint step vectors |

### 2.2 Interval / FASTA / read-stat tools

| line | test name | script | record-level? | class | key assertions |
|---|---|---|---|---|---|
| 1054 | `fastgatk-annotate-intervals-contract` | `verify_annotate_intervals.py` | TSV rows | (a) | `:66` rows incl. `chr2 1 4 NaN`; java-compared `:98-100,:151-153` |
| 1057 | `fastgatk-count-bases-in-reference-contract` | `verify_count_bases_in_reference.py` | text lines | (a) + **(c)** | `:57,:92` java-verified; **`:84-85` 3-field `.interval_list` accepted by native, GATK exits 3**; **`:100,:108` plain-gzip inputs, GATK exits 2** |
| 1060 | `fastgatk-compare-references-contract` | `verify_compare_references.py` | table text | (a) | compares to checked-in GATK fixture tables |
| 1063 | `fastgatk-check-reference-compatibility-contract` | `verify_check_reference_compatibility.py` | table text | (a) | GATK fixtures + live java; `:73-78` native-only `.bgz` extension |
| 1066 | `fastgatk-fasta-reference-maker-contract` | `verify_fasta_reference_tools.py` | byte compare | (a) | **guarded by `if java_oracle:` — no jar presence assertion** |
| 1072 | `fastgatk-shift-fasta-contract` | `verify_shift_fasta.py` | byte compare | (a) | GATK fixtures + live java |
| 1075 | `fastgatk-index-feature-file-contract` | `verify_index_feature_file.py` | `.idx` bytes | (a) | `:87,:125,:171,:206,:227` `read_bytes() == java_index.read_bytes()` |
| 1080 | `fastgatk-count-reads-contract` | `verify_count_reads.py` | counts + native text | (a) + **(c)** | counts java-verified `:421-585`; **`:143-150` rejects `-ip -1` that GATK accepts**; **`:193-212` accepts unindexed gzip GATK rejects** |
| 1085 | `fastgatk-flag-stat-contract` | `verify_flag_stat.py` | 12-line text block | (a) + **(c)** | `:458` `java_lines == EXPECTED`; **`:136-155` gz inputs**, **`:259-266` negative padding** |
| 1090 | `fastgatk-split-intervals-contract` | `verify_split_intervals.py` | interval rows | (a) + **(c)** | `:159` matches GATK exactly (confirmed); **`:140-145` gz BED+`-XL`: GATK exits 2** |
| 1095 | `fastgatk-filter-intervals-contract` | `verify_filter_intervals.py` | interval rows | (a) | `:63,:168-169,:235-237` java-verified (oracle behind `if java.exists()`) |
| 1098 | `fastgatk-preprocess-intervals-contract` | `verify_preprocess_intervals.py` | interval rows | (a) | `:71,:108` java-verified (same guard) |

### 2.3 HaplotypeCaller contracts

| line | test name | script | record-level? | class | key assertions |
|---|---|---|---|---|---|
| 1219 | `fastgatk-indel-contract` | `verify_indel.py` | **yes** | **(a)** | the corrected test. All record expectations confirmed byte-identical to GATK (§3.1): `:56,:87-88,:102,:116-120` (`17:11 … DP=0 END=32`), `:172` POS=10, `:178` "no END on a concrete gVCF variant record" |
| 1222 | `fastgatk-hc-ploidy-contract` | `verify_hc_ploidy.py` | yes | (b) | `:52` PL len == ploidy+1, `:54` AN==ploidy, `:77` triploid candidate PL len 10 — GATK-consistent but never compared; `:57` native-only manifest flag |
| 1225 | `fastgatk-hc-likelihood-filter-contract` | `verify_hc_likelihood_filter.py` | yes | (a)/(b) | `:117-119` `native_line == gatk_line`, `values == gatk_values`, `info == gatk_info` → **oracle**; but `:64-68` (native *default* threshold path) are native-only. Confirmed: GATK at default `-stand-call-conf` also emits `69067 T G 33.48` (§3.4) |
| 1230 | `fastgatk-hc-variant-annotations-contract` | `verify_hc_variant_annotations.py` | **yes** | **(a)** | `:64-73` requires `MQ/MLEAC/MLEAF/QD/FS/SOR/MQRankSum/ReadPosRankSum/BaseQRankSum/ExcessHet` on the sentinel and `QD==round(QUAL/DP,2)`. Native and GATK produce identical INFO on all three records of this fixture (§3.4) |
| 1235 | `fastgatk-hc-genotype-priors-contract` | `verify_hc_genotype_priors.py` | **yes** | **(b) high risk** | `:128` `PL == "174,87,75,87,0,75"`, `:130` `GQ=="49"`, `:163` `"174,87,…,174"`, `:181` `GT=="0/1/2"`, `:183-184` poly PL/GQ — all native-observed literals on a synthetic SAM. GATK is used only for `GenotypePriorCalculator` vectors (`:24-58`), **never to produce a record** |
| 1340 | `fastgatk-hc-chr17-independent-regions` | `verify_hc_chr17_independent_regions.py` | yes | (a) | `:91` `("17",69067,"T","G")`, `:105` `("17",69298,"A","T")` with `--dont-use-soft-clipped-bases` — both **confirmed identical to GATK** (§3.4); `:97,:111` byte-determinism |
| 1369 | `fastgatk-hc-rcm-realignment-contract` | `verify_hc_rcm_realignment.py` | yes | (b) | `:51` `baseline == candidate` byte equality between two native modes; `:54-58` telemetry |
| 1374 | `fastgatk-flow-hmer-option-contract` | `verify_flow_hmer_option.py` | yes | (b) | synthetic fixture; PL/AD/GT literals |
| 1389 | `fastgatk-hc-streaming-contract` | `verify_hc_streaming.py` | yes | (b) | `:65` `normal == streamed` byte equality; `:67` exactly 2 records, contig order `["chr1","chr2"]`; `--stream-by-contig` is a native-only flag |
| 1394 | `fastgatk-hc-region-streaming-contract` | `verify_hc_region_streaming.py` | yes | (b), documented | `:50,:162` byte equality; `:86` key-set equality (already annotated in-file as fixture-scoped); GATK-compared only at `:184` |
| 1409 | `fastgatk-hc-gvcf-stream-overlap-diagnostic` | `verify_gvcf_stream_overlapping_indels_gatk_oracle.py` | yes | (a) | **oracle by construction (name lacks `contract`) — out of target** |
| 1524 | `fastgatk-hc-multi-input-contract` | `verify_hc_multi_input.py` | yes | (b) + oracle | `:226` bounded-key check; GATK `HaplotypeCaller`/`SelectVariants` used |
| 1541 | `fastgatk-activity-region-controls-contract` | `verify_activity_region_controls.py` | no | (b) | `:39` `activity_filter_size==50`, `:40` `…==100` — GATK-derived defaults echoed as native telemetry |
| 1546 | `fastgatk-overlapping-quality-correction-contract` | `verify_overlapping_quality_correction.py` | no | (b) | engine-level |
| 1551 | `fastgatk-real-assembly-graph-contract` | `verify_real_assembly_graph.py` | AD only | (b) | |
| 1556 | `fastgatk-hc-chr20-real-contract` | `verify_hc_chr20_real_contract.py` | **yes** | (a) | actually byte/oracle compares against pinned GATK rows+header `:237-293`; out of target |
| 1615 | `fastgatk-vqsr-scatter-joint-contract` | `verify_vqsr_scatter_joint.py` | yes | (b) | QUAL/FILTER/PL/AD/GT/QD/culprit literals |
| 1124 | `fastgatk-bqsr-contract` | `verify_bqsr.py` | table rows | (a) | GATK `BaseRecalibrator` used; `:126-127,:257` table literals |
| 1204 | `fastgatk-analyze-covariates-contract` | `verify_analyze_covariates.py` | plot files | (a) | GATK `AnalyzeCovariates` used |

### 2.4 gVCF materialization / reblocking (highest-risk family)

| line | test name | script | record-level? | class | key assertions |
|---|---|---|---|---|---|
| 1653 | `fastgatk-genotype-gvcf-contract` | `verify_genotype_gvcf.py` | **yes (1011 lines)** | **(c)** + (b) | **`:391` `len(star_records)==1 and "*" in star_records[0][4]` — GATK emits 0 records (§3.3)**; `:353-367` multi-sample GT/GQ/AC/AN/AF from PL (GATK-consistent, but the fixture header omits `GQ` so GATK cannot even run it — §3.5); `:87,:995` `RCQ=`/`RCP=` presence; `:232-234` tag/FORMAT exclusions; `:762` GQ==60; `:907,:951` record counts; GATK is used **only** as a reader (`SelectVariants` `:104`) |
| 1750 | `fastgatk-reblock-gvcf-contract` | `verify_reblock_gvcf.py` | **yes (609 lines)** | **(c)** + (b) | **`:64,:82-83` third record with `END=20` — GATK drops it (§3.2)**; **`:283-287` triploid `0/1/1` + compacted PL — GATK reblocks to `0/0/0` (§3.2)**; **`:77`,`:80` positional FORMAT assumption (GATK order is `GT:AD:DP:GQ:PL`, native `GT:DP:AD:PL:GQ`) (§3.2)**; `:556-560` RAW_MQandDP/MQ_DP = **(a)**, source-confirmed `ReblockGVCF.java:1000-1018`; `:323-332` multi-sample block merge (GATK refuses multi-sample ReblockGVCF entirely, §3.2); `:374-386` deletion→block END = (a), source-confirmed `:625`; `:425-435` trim/gap block = (a), source-confirmed `addRefBlockIfNecessary` `:703-746`; `:485-497` annotation cleanup; `:583` FILTER tolerant `in {"PASS","."}` |
| 2287 | `fastgatk-combine-gvcfs-contract` | `verify_combine_gvcfs.py` | **yes** | (b), fixture-scoped | `:70-71` `len(records)==2`, `"END=7" in records[0] and "DP=10" in records[0]` — native coalesces blocks 1-3 + 4-7; GATK on the identical **sample-less** fixture emits 3 records and no INFO (§3.6). `:325-327` `--break-bands-at-multiples-of 5` → `[1,5,10]`/`END 4,9,12` = **(a)** (matches GATK exactly); `:339` BP resolution 12 records, no END = **(a)** |
| 2294 | `fastgatk-combine-gvcfs-gatk-oracle` | — | yes | (a) | oracle, out of scope |

### 2.5 Mutect2 / somatic

| line | test name | script | record-level? | class | key assertions |
|---|---|---|---|---|---|
| 2042 | `fastgatk-mutect2-contract` | `verify_mutect2.py` | **yes** | **(c)** samples + (b) | **`:650` `synthetic_header[9:] == ["NORMAL","OTHER","TUMOR"]` — GATK emits `NORMAL OTHER OTHER_NORMAL TUMOR` (§4.1)**; `:192/:210` call positions (b, value matches GATK); `:76-77` QUAL `.` + NALOD/NLOD (a); `:647/:678/:684` ALT order / `0/1/2` / `DP=0` — **not GATK-reproducible** (`--min-depth`, `--min-alt-support`, `--min-kmer-count`, `--max-haplotype-depth` are not GATK Mutect2 options) |
| 2317 | `fastgatk-filter-mutect-calls-contract` | `verify_filter_mutect_calls.py` | **yes** | **(c)** | **`:156-157` `base_qual`/`map_qual` from `Number=1` MBQ/MMQ — GATK fires no filter (§4.1)**; **`:64` `FAIL` + `:664` `low_tlod` — `--min-tlod` is not a GATK option and `low_tlod` is not a GATK constant (§4.1)**; filter *names* elsewhere = (a) |
| 2051 | `fastgatk-mutect2-itr-artifact-gatk-contract` | `verify_mutect2_itr_artifact_gatk_contract.py` | no (sidecar) | (b) | counters |
| 2060 | `fastgatk-mutect2-pcr-overlap-gatk-contract` | `verify_mutect2_pcr_overlap_gatk_contract.py` | no | (b) | `:71-85` counters |
| 2068 | `fastgatk-mutect2-independent-mates-gatk-contract` | `verify_mutect2_independent_mates_gatk_contract.py` | no | (b) | counters |
| 2076 | `fastgatk-mutect2-gvcf-reference-blocks-gatk-contract` | `verify_mutect2_gvcf_reference_blocks_gatk_contract.py` | yes | (a) | `:71 native_rows == gatk_rows`; `:77-80` literal is redundant with that |
| 2083 | `fastgatk-mutect2-bp-resolution-gatk-contract` | `verify_mutect2_bp_resolution_gatk_contract.py` | yes | (a) | `:60 observed == expected` |
| 2091 | `fastgatk-mutect2-gvcf-eventmap-gatk-contract` | `verify_mutect2_gvcf_eventmap_gatk_contract.py` | yes | (a) | `:108-116` GATK equality |
| 2099 | `fastgatk-mutect2-reference-confidence-eventmap-full-gatk-contract` | `…full_gatk_contract.py` | yes | (a) | `:86 observed == expected` |
| 2106/2111 | `…gvcf-matched-normal…` / `…gvcf-multitumor…` | both | yes (negative) | (a) | both tools must fail with the same semantics |
| 2158 | `fastgatk-mutect2-mismapping-rate-boundary` | `verify_mutect2_mismapping_rate_boundary.py` | no | (a) | GATK error text + jshell probe |
| 2171 | `fastgatk-mutect2-tlod-formula` | `verify_mutect2_tlod_formula.py` | no | (b) | hand-transcribed formula |

### 2.6 Variant manipulation

| line | test name | script | record-level? | class | key assertions |
|---|---|---|---|---|---|
| 1791 | `fastgatk-select-variants-contract` | `verify_select_variants.py` | **yes** | **(c)** + (b) | **`:820` `assert ref_only_records == []` — GATK keeps the record with `ALT=.` (§4.4.1)**; **`:420` `vc.getType() == VariantContext.Type.SNP` → `["1","2"]`, GATK selects `[]`**; `:412,:413,:416,:423` JEXL expressions GATK **rejects** (rc 2); **`:116` `…split(":")[-1]=="80"` and `:790-795` positional FORMAT reads** (native `GT:AD:PL:GQ` vs GATK `GT:AD:GQ:PL`); `:112-115,:122-124` values GATK-identical; `:117-118` substring `AF=0.5` is format-tolerant and does **not** pin `0.5` vs `0.500` |
| 1826 | `fastgatk-gather-vcfs-contract` | `verify_gather_vcfs.py` | yes (header lines) | **(b)** | **`:125-126,:130` two `##GatherVcfs.comment=` lines (`comments==2`) — GATK writes exactly one for two `-CO`/`--COMMENT` values (§4.4.7)**; `:70`,`:188` positions `[1,2]` (oracle-checked `:150-199`); `:403` block kept by `-L` = (a) confirmed |
| 1837 | `fastgatk-left-align-trim-contract` | `verify_left_align.py` | yes | (a) + **(b) positional** | `:159-160` symbolic records compared to GATK incl. sample column; `:344-360` triploid compared by field name — values identical (`GQ 30/60`); **but `:111-113` `split(":")[-1]` and `:326-327` positional FORMAT lists pin native's `GT:AD:PL:GQ` order** (§4.4.4) |
| 1865 | `fastgatk-variant-filtration-contract` | `verify_variant_filtration.py` | **yes** | **(c)** + (b) | **`:504-505` `AS_FilterStatus=LowASQD,PASS` / `PASS,LowASQD` — GATK emits `AS_FilterStatus=SITE\|SITE` and does not allele-filter (§4.4.3)**; **`:77` `"Het;LowDP;LowQD;LowQual"` — GATK exits 3 on that command line (`NumberFormatException: "1.0"`); its 3-filter subset gives the same ordering, so the `LowQD` term is native-only (§4.4)**; the comparator `normalized_filter_records:33-34` also sorts filter names and drops INFO, so ordering is not oracle-verified; `:114-115,:619-630` are oracle-checked |
| 1894 | `fastgatk-apply-vqsr-contract` | `verify_apply_vqsr.py` | yes | (a) + (b) | Real GATK `ApplyVQSR` oracle at `:423-437` and `:486-501`; `:95` tranche `+` suffix and `:348-349` `NA` padding confirmed as **GATK** behaviour (`ApplyVQSR.java:197,458`); `:93-96,:256-258,:294-349` literals sit on fixtures the script documents as **not valid Java oracle inputs** (`:444-449`), i.e. (b); `:410` `BACKEND_UNAVAILABLE` is a native fail-closed string |
| 1936 | `fastgatk-variant-recalibrator-contract` | `verify_variant_recalibrator.py` | **yes (INFO/model)** | **(c)** + (b) | **`:84,:422,:437,:456` `culprit=full-covariance-gmm` / `AS_culprit=full-covariance-gmm` — the string does not occur in the GATK jar (0 hits); GATK writes `culprit=MQ` (§4.4.2)**; **`:276-277` lowercase `e+01` where GATK's `%.16E` writes `E+01`**; `:274` `#:GATKReport.v1.1:8`, `:85`, `:86-88`, `:233-235`, `:372` `gatk_random_seed` confirmed (a) |
| 1989 | `fastgatk-sort-sam-contract` | `verify_sort_sam.py` | yes | (a) + (b) | `:55` `@HD VN:1.6 SO:coordinate` and `:56-57,:84` oracle-checked; `:98-100` `@PG` `PN/VN/CL` is native-rendered from CLI flags |
| 2014 | `fastgatk-mark-duplicates-contract` | `verify_mark_duplicates.py` | yes (SAM/metrics) | (a) + **(b)** | **`:191,:193` histogram rows `["1","0","0","1"]` / `["2","1","1","0"]` — GATK writes `1.0`/`2.0` (§4.4.6)**; `:190` `## HISTOGRAM java.lang.Double` and `:118` `DT:Z:LB` confirmed (a); `:144-146` native `@PG` |
| 2499 | `fastgatk-variants-to-table-contract` | `verify_variants_to_table.py` | yes (TSV) | (a) | 11 byte-identical GATK TSV comparisons; `:305-306`,`:429-431`,`:264` confirmed identical |
| 2510 | `fastgatk-variant-eval-contract` | `verify_variant_eval.py` | yes (report) | (a) + (b) | GATK `VariantEval` oracle `:216-234`; `:100-112,:134,:469-477` parse native's own `##table=` schema (b) |
| 2535 | `fastgatk-validate-variants-contract` | `verify_validate_variants.py` | no (native TSV) | **(b)** | `:52,:71` assert a native report TSV (`"\tPASS\t4\t3\t0\t0"`); GATK `ValidateVariants` declares **no OUTPUT argument**, so this cannot be a GATK contract; `:188` native text vs GATK's `"not observed at all"` |
| 2322 | `fastgatk-optional-boolean-contract` | `verify_optional_boolean_contract.py` | no | (b) | `:77` 20 tools × 3 spellings must print `"invalid boolean"`; no jar invocation |
| 2435 | `fastgatk-genomicsdb-import-contract` | `verify_genomicsdb_import.py` | yes | (b) | `:208` record-index row; `:372-374` `len(row_fields)==11`, GQ literal `"7"` |
| 2490 | `fastgatk-genomicsdb-native-storage-boundary` | `verify_genomicsdb_native_storage_boundary.py` | metadata | (b)/(c deliberate) | native-only sparse workspace GATK must reject |

### 2.7 Coverage / pileup / CNV

| line | test name | script | record-level? | class | key assertions |
|---|---|---|---|---|---|
| 2551 | `fastgatk-get-pileup-summaries-contract` | `verify_get_pileup_summaries.py` | table | (a) | `:55-58` GATK table format (`PileupSummary.java`, `TableWriter.java`) |
| 2567 | `fastgatk-gather-pileup-summaries-contract` | `verify_gather_pileup_summaries.py` | row order | (a) | `:53-55` — `GatherPileupSummaries` sorts by first record |
| 2561 | `fastgatk-calculate-contamination-contract` | `verify_calculate_contamination.py` | table | (a) | `:54,:153-157` GATK byte-identical |
| 2801 | `fastgatk-collect-allelic-counts-contract` | `verify_collect_allelic_counts.py` | TSV rows | (b) | `:79-81` counts/rows |
| 2809 | `fastgatk-collect-f1r2-counts-contract` | `verify_collect_f1r2_counts.py` | dual | (a)/(b) | oracle at `:242-271`; extras `:140-152` |
| 2583 | `fastgatk-collect-read-counts-contract` | `verify_collect_read_counts.py` | TSV rows | (b) | `:64-65,:122,:150,:161` |
| 2604 | `fastgatk-denoise-read-counts-contract` | `verify_denoise_read_counts.py` | values | (b) | `:113-116,:164` |
| 2633 | `fastgatk-create-read-count-panel-of-normals-contract` | `verify_create_read_count_panel_of_normals.py` | HDF5 | (a)/(b) | `:68` native format string; Java roundtrip `:160` |
| 2667 | `fastgatk-call-copy-ratio-segments-contract` | `verify_call_copy_ratio_segments.py` | call column | (b) | `:55` `["0","0","0","+","+","-"]` |
| 2712 | `fastgatk-model-segments-contract` | `verify_model_segments.py` | dual | (a)/(b) | GATK-verified het sidecars `:229-232`; segment counts `:70-73` are a native extension |
| 2814 | `fastgatk-depth-of-coverage-contract` | `verify_depth_of_coverage.py` | table | (a) | byte-compared to GATK; `:352,:361` are a native-only `COUNT_FRAGMENTS` extension |
| 2817 | **`fastgatk-depth-of-coverage-multisample-contract`** | `verify_depth_of_coverage_multisample.py` | table | **(c)** | `:52` `…Average_Depth_S1,Depth_for_S1,Average_Depth_S2,Depth_for_S2` and `:53-59` `17:69000,3,2.00,2,1.00,1` — **GATK prints one `Average_Depth_sample` column and `1.50` (§4.1)** |
| 2599 | `fastgatk-hdf5-simple-count-collection-metadata` | `verify_hdf5_simple_count_collection.py` | HDF5 | (a)/(b) | |
| 2575 | `fastgatk-learn-read-orientation-model-contract` | `verify_learn_read_orientation_model.py` | dual | (a)/(b) | GATK row-order equality `:185-197` |
| 2798 | `fastgatk-gather-tranches-contract` | `verify_gather_tranches.py` | dual | (a) | GATK byte compare `:109,:136` |

---

## 3. Confirmations performed (literal output)

All comparisons were made with the **same fixtures the tests use**, run through both the native binary
the test invokes and pinned GATK 4.6.2.0.

### 3.1 `verify_indel.py` — class (a) confirmed (this test is *correct*)

Region `17:1-32`, the test's own 32 bp non-repetitive reference and 4×`10M1I10M` reads
(native run with `-L 17:1-32`; GATK on the same data as BAM, no `-L` needed for a 32 bp contig):

```
NATIVE: 17 1 A <NON_REF> . . END=9 | GT:DP:GQ:MIN_DP:PL 0/0:4:12:4:0,12,179
NATIVE: 17 10 A AC,<NON_REF> 147.09 . DP=4;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=14400,4 | GT:AD:DP:GQ:PL:SB 1/1:0,4,0:4:12:161,12,0,161,12,161:0,0,4,0
NATIVE: 17 11 G <NON_REF> . . END=32 | GT:DP:GQ:MIN_DP:PL 0/0:0:0:0:0,0,0
GATK  : 17 1 A <NON_REF> . . END=9 | GT:DP:GQ:MIN_DP:PL 0/0:4:12:4:0,12,179
GATK  : 17 10 A AC,<NON_REF> 147.09 . DP=4;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=14400,4 | GT:AD:DP:GQ:PL:SB 1/1:0,4,0:4:12:161,12,0,161,12,161:0,0,4,0
GATK  : 17 11 G <NON_REF> . . END=32 | GT:DP:GQ:MIN_DP:PL 0/0:0:0:0:0,0,0
```

The three riskiest assertions in that file — `:116-120` (a record at POS 11 with `DP == "0"` and
`END=32`), `:110` (`block_pl[0] == 0`), and `:178` (`"END" not in info` for the concrete
`17:10 A>AC,<NON_REF>` record) — are all genuine GATK behaviour. Note the concrete variant record
carries `RAW_MQandDP`, *not* `END`, on both sides.

### 3.2 `verify_reblock_gvcf.py` — class (c) ×2 and class (b) ×2 confirmed

Fixture A (script lines 41-51), exactly the script's flags
`-GQB 20 -GQB 100 --drop-low-quals --rgq-threshold 10 --floor-blocks`, both tools given the same
100 bp reference (adding `-R` to native does not change its output):

```
--- A NATIVE(-R)
   POS 1  REF A ALT <NON_REF> FILTER . INFO END=10                   FMT GT:DP:AD:GQ       VAL 0/0:9:10,0:0
   POS 12 REF A ALT G,<NON_REF> FILTER . INFO DP=20;RAW_MQandDP=72000,20;RAW_GT_COUNT=0,1,0
                                                                    FMT GT:DP:AD:PL:GQ    VAL 0/1:20:12,8,0:50,0,80,99,99,20:20
   POS 20 REF C ALT <NON_REF> FILTER . INFO DP=5;END=20              FMT GT:DP:GQ          VAL 0/0:5:0
--- A GATK
   POS 1  REF A ALT <NON_REF> FILTER PASS INFO END=10                FMT GT:DP:GQ          VAL 0/0:9:0
   POS 12 REF A ALT G,<NON_REF> FILTER . INFO DP=20;RAW_GT_COUNT=0,1,0;RAW_MQandDP=72000,20
                                                                    FMT GT:AD:DP:GQ:PL    VAL 0/1:12,8,0:20:20:50,0,80,99,99,20
```

* (c) **`assert len(records) == 3` (`:64`) + `records[2][4] == "<NON_REF>"` / `"END=20"` (`:82-83`)**
  pins the third record. GATK emits only 2 records. Cause: `ReblockGVCF.java:542-547`
  (`lowQualVariantToGQ0HomRef` returns `null` when `dropLowQuals` is set and
  `isMonomorphicCallWithAlts` is false — here `PL=[2,0,40]`, so `minElementIndex != 0`), so the
  site is **dropped**, not converted to a block.
* (b) **`records[1][9].split(":")[3]` (`:77`) and `…split(":")[-1] == "20"` (`:80`)** depend on
  native's FORMAT order `GT:DP:AD:PL:GQ`. GATK's order is `GT:AD:DP:GQ:PL`; index 3 is `PL` in
  native and `GQ` in GATK, and the last field is `GQ` in native but `PL` in GATK. Adopting GATK's
  FORMAT ordering silently breaks both.
* Also visible: native writes `FILTER=.` on the merged reference block where GATK writes `PASS`.
  Not asserted at `:65-66`, and `:583` is deliberately tolerant (`in {"PASS","."}`), so this
  divergence is currently *not* pinned.

Triploid fixture (script lines 269-273), no extra flags:

```
--- T NATIVE(-R)
   POS 30 REF A ALT C,<NON_REF> INFO DP=17;RAW_MQandDP=61200,17;RAW_GT_COUNT=0,1,0
                                 FMT GT:DP:AD:PL:GQ     VAL 0/1/1:17:12,5,0:0,10,20,30,100,110,120,160,170,190:10
--- T GATK
   POS 30 REF A ALT <NON_REF>    INFO END=30
                                 FMT GT:DP:GQ:MIN_DP:PL VAL 0/0/0:17:10:17:0,10,20,30
```

* (c) **`:283-287`** (`ALT == "C,<NON_REF>"`, `GT == "0/1/1"`, `AD == "12,5,0"`,
  `PL == "0,10,20,30,100,110,120,160,170,190"`) pins native keeping the call. GATK's
  `shouldBeReblocked` (`ReblockGVCF.java:515-531`) reblocks the site because the minimum-PL
  genotype is hom-ref (`genotypeHasConcreteAlt` false), producing a reference block with
  `GT 0/0/0`, `END=30` and 4 PLs. A parity fix for this is exactly what the test forbids.
* (c/b) **`:318-332` multi-sample** — pinned GATK 4.6.2.0 **rejects multi-sample input outright**:
  `A USER ERROR has occurred: Bad input: ReblockGVCF can take multiple input GVCFs, but they must
  be non-overlapping shards from the same sample. Found samples [S1, S2]`. The whole multi-sample
  block (block merge, per-sample MIN_DP/DP/GQ) is therefore a native-only extension with no GATK
  reference behaviour.

Confirmed **(a)** in the same file, from GATK source read directly:

```
ReblockGVCF.java:1005-1017
  final int rawMqValue = sourceVC.hasAttribute(RAW_RMS_MAPPING_QUALITY_DEPRECATED) ?
          (int)Math.round(...RAW_RMS_MAPPING_QUALITY_DEPRECATED, 0.0) :
          (int)Math.round(...MQ... * ...MQ... * ...DP..., 0.0);
  ...put(RAW_MAPPING_QUALITY_WITH_DEPTH_KEY, rawMqValue + "," + DP);
```
→ `:556` `RAW_MQandDP=25000,10` (MQ=50,DP=10 → 50²·10) and `:557-558`
`RAW_MQandDP=123,10` + `MQ_DP=10` are genuine GATK output. Likewise `:374-386` (END on a
*block* produced from a dropped deletion) matches `changeCallToHomRefVersusNonRef` `:618-625`,
and `:425-435` (trim + gap block) matches `addRefBlockIfNecessary` `:703-746`.

### 3.3 `verify_genotype_gvcf.py:391` — class (c) confirmed

The script's own fixture (`:378-381`), `chr1 2 . A *,G,<NON_REF> … GT:DP:AD:PL 0/1:20:12,8,0,0:0,0,100,100,100,100,100,100,100,100`;
native is run exactly as the script runs it (no `-R`), GATK with `-R` (required):

```
=== exact  native(no -R) rc=0 | native(-R) rc=0 | GATK(-R) rc=0
    NATIVE no-R   : 1 record
      chr1 2 . A *,G 0 PASS DP=20;MLEAC=0,0;MLEAF=0,0;AC=0,0;AN=2;AF=0,0;RCQ=0.00434077;RCP=0.999001;ExcessHet=-0
      GT:DP:AD:PL:GQ  0/0:20:12,8,0:0,0,100,100,100,100:0
    GATK          : 0 records
=== with a variant MLE (same ALT set, PL min at 0/1)
    NATIVE        : 1 record   chr1 3 . A *,G 0 PASS … RCQ=70;RCP=1e-07;QD=0 … GT 0/1 … PL 100,0,100,100,100,100 GQ 99
    GATK          : 0 records
=== GATK -all-sites on the exact fixture
    GATK          : 1 record   chr1 2 . A . 127.78 . DP=20;MLEAC=.;MLEAF=.   GT ./. 
```

So for this fixture GATK has no `*,G` record in any mode: default = no record, `-all-sites` = a
record with `ALT=.` / `GT=./.`. Native emits `ALT=*,G` with `GT=0/0`. The test's
`assert len(star_records) == 1 and "*" in star_records[0][4]` (and
`compatibility["spanning_deletion_nonvariant_set"]`) therefore pins native-only behaviour. GATK's
gate is `GenotypeGVCFs.java:325-329` → `GenotypeGVCFsEngine.regenotypeVC` (`:154-196`), which
returns `null` unless the regenotyped site is polymorphic in samples or `-all-sites`/`-include-non-variants`
is set.

### 3.4 HC annotations / threshold / soft-clip — class (a) confirmed

`verify_hc_variant_annotations.py` fixture (`NA12878.chr17_69k_70k` + 1 Mb chr17 ref, `-L 17:69000-70000`,
`-stand-call-conf 0`, native `--min-depth 1 --min-alt-support 1`):

```
NATIVE: 69067 T G 33.48 . AC=2;AF=1.00;AN=2;DP=1;ExcessHet=0.0000;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;QD=33.48;SOR=1.609
NATIVE: 69368 G C 670.64 . AC=1;AF=0.500;AN=2;BaseQRankSum=-3.891;DP=42;ExcessHet=0.0000;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=59.56;MQRankSum=1.096;QD=15.97;ReadPosRankSum=-2.008;SOR=0.242
NATIVE: 69631 C T 677.64 . AC=1;AF=0.500;AN=2;BaseQRankSum=4.461;DP=37;ExcessHet=0.0000;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=57.83;MQRankSum=0.518;QD=19.36;ReadPosRankSum=0.877;SOR=0.399
GATK  : (identical, all three records, INFO and QUAL byte-for-byte)
```

`verify_hc_chr17_independent_regions.py` / `verify_hc_likelihood_filter.py`:

```
NATIVE default (17:69000-69100): 1 record -> [('69067','T','G','33.48')]
GATK   default (-stand-call-conf 30): 1 record -> [('69067','T','G','33.48')]
NATIVE --dont-use-soft-clipped-bases: [('69067','T','G','33.48'),('69298','A','T','81.64'),('69368','G','C','701.64'),('69631','C','T','677.64')]
GATK   --dont-use-soft-clipped-bases: [('69067','T','G','33.48'),('69298','A','T','81.64'),('69368','G','C','701.64'),('69631','C','T','677.64')]
```

### 3.5 `verify_genotype_gvcf.py` multi-sample fixture is not GATK-runnable

`GenotypeGVCFs` on the script's `multi.g.vcf` fixture fails before producing output:

```
GATK-multi rc=3 (0 data records)
Caused by: java.lang.IllegalStateException: Key GQ found in VariantContext field FORMAT at chr1:1
but this key isn't defined in the VCFHeader. We require all VCFs to have complete VCF headers by default.
```

The fixture header (`:328-338`) declares `GT/DP/AD/PL` but not `GQ`, while GATK's output adds GQ.
The test's `assert "GQ" in format_names` (`:355`) is therefore a native-only claim as written: the
fixture cannot be validated against GATK until `##FORMAT=<ID=GQ,...>` is added.

### 3.6 `verify_combine_gvcfs.py` block coalescing

Schedule fixture, `--break-bands-at-multiples-of 5` and `--convert-to-base-pair-resolution`
(**class (a)** — match exactly):

```
=== split  NATIVE rc=0 | GATK rc=0
    NATIVE: 3 records  chr1 1 … DP=10;END=4 | chr1 5 … DP=10;END=9 | chr1 10 … DP=10;END=12
    GATK  : 3 records  chr1 1 … END=4        | chr1 5 … END=9        | chr1 10 … END=12
=== bp     NATIVE rc=0 | GATK rc=0   (both 12 records at POS 1..12, no END)
```

Shard-merge fixture (script `:42-47`, sample-less header):

```
=== sample-less NATIVE rc=0 : 2 records
      chr1 1 . A <NON_REF> . . DP=10;END=7
      chr1 8 . A G         . . DP=10
=== sample-less GATK   rc=0 : 3 records
      chr1 1 . A <NON_REF> . . END=3
      chr1 4 . A <NON_REF> . . END=7
      chr1 8 . A <NON_REF> . . .
```

`:70-71` (`len(records)==2`, `"END=7" in records[0] and "DP=10" in records[0]`) pins native's
coalescing of 1-3 + 4-7 and its INFO/DP retention. GATK does neither on this input (it has no
genotypes to merge on, and it rewrites the concrete ALT at POS 8 to `<NON_REF>`). Classification:
**(b), fixture-scoped** — the assertion is not a GATK contract, but a GATK run cannot produce a
comparable expectation either, so it cannot be "fixed" by pointing at GATK.

---

## 4. Confirmed findings that came from delegated parallel audits

The audit was fanned out across three delegated sub-audits (non-VCF tools; somatic/CNV/pileup;
variant manipulation), each of which ran the pinned jar on the scripts' own fixtures. All three
completed and reported; their executed evidence is reproduced here because the failures are of the
same shape — a native-internal expectation asserted as if it were a contract. §4.4 is the
variant-manipulation batch. The task instructions said not to run the full `ctest` suite; none was
run — every confirmation below is a targeted invocation.

### 4.1 Confirmed class-(c) divergences pinned as contracts

1. `verify_filter_mutect_calls.py:156-157` — `assert "base_qual" in defaults_records[0][6]`,
   `assert "map_qual" in defaults_records[1][6]`, on a fixture whose header declares
   `MBQ`/`MMQ` as `Number=1` (`:29-30`, `:139-140`). GATK on the line-for-line same VCF:
   `10 PASS`, `11 PASS`, `12 position`, `13 fragment`, `14 PASS`. Cause
   (`BaseQualityFilter.java`, `MappingQualityFilter.java`): GATK reads MBQ/MMQ as `Number=R` and
   `skip(1)`/`remove(0)` the ref entry, so a `Number=1` value leaves no ALT value and no filter
   fires (`MBQ=30,19` → `base_qual`; `MMQ=30,29` → `map_qual`).
2. `verify_filter_mutect_calls.py:63-65,664` — `"\tFAIL\tTLOD=-1"` and
   `AS_FilterStatus=low_tlod`. `min-tlod is not a recognized option` in GATK 4.6.2.0;
   `M2FiltersArgumentCollection.java` has no TLOD filter (GATK's is `weak_evidence`), `low_tlod`
   is not a constant, and `FAIL` is reserved for "all alleles filtered for different reasons"
   (`Mutect2FilteringEngine.java:227`).
3. `verify_depth_of_coverage_multisample.py:52-59` — native
   `Locus,Total_Depth,Average_Depth_S1,Depth_for_S1,Average_Depth_S2,Depth_for_S2` with
   `17:69000,3,2.00,2,1.00,1`; GATK prints `Locus,Total_Depth,Average_Depth_sample,Depth_for_S1,Depth_for_S2`
   with `17:69000,3,1.50,2,1` (`CoverageOutputWriter.java:222-231,294`; `DoCOutputType.java:8`).
   The single-sample sibling `verify_depth_of_coverage.py:65` expects `Average_Depth_sample`, so
   native is also internally inconsistent between the two paths.
4. `verify_mutect2.py:650` — native `#CHROM … NORMAL OTHER TUMOR`; GATK for the same two inputs
   `#CHROM … NORMAL OTHER OTHER_NORMAL TUMOR` (GATK lists every non-normal header sample).
5. `verify_count_reads.py:143-150` / `verify_flag_stat.py:259-266` — assert native **rejects**
   `-ip -1`. GATK accepts and honours it: `CountReads -L 17:69990-69990 -ip -1` → rc 0, count `11`.
6. `verify_count_reads.py:193-212`, `verify_flag_stat.py:136-155`,
   `verify_count_bases_in_reference.py:100,108`, `verify_split_intervals.py:115-145` — assert
   native consumes **unindexed, non-BGZF** `.gz` interval/BED files; GATK exits 2
   (`An index is required but was not found …`) and cannot index them either.
7. `verify_count_bases_in_reference.py:84-85` — asserts a **3-field** `.interval_list` works;
   GATK exits 3 (`htsjdk.tribble.TribbleException: Invalid interval record contains 3 fields`).

### 4.2 Confirmed native error-text/exit-code pinning (class (b))

`verify_count_reads.py:155,182,254,416`, `verify_flag_stat.py:134,271,283,303,314,353,390`,
`verify_launcher_contract.py:143,152,158` pin native's private diagnostic vocabulary and exit
code `2`, where GATK uses exit `1` and its own wording (`Unrecognized read filter name: …`,
`Argument "/minimum-mapping-quality" is only valid when the argument "MappingQualityReadFilter" is specified`,
`Argument interval-merging-rule has a bad value: TOUCHING …`). Same shape as the already-catalogued
`--alleles` exit-code divergence.

### 4.3 Structural risk: tests that *look* like oracles but are not

* At least 11 scripts wrap the GATK comparison in `if java.exists() and jar.exists():` /
  `if java_oracle:` **without asserting the jar exists** (`verify_filter_intervals.py:177`,
  `verify_preprocess_intervals.py:137`, `verify_annotate_intervals.py:89,134,182,208`,
  `verify_fasta_reference_tools.py:105-236`, `verify_count_bases_in_reference.py:121`,
  `verify_compare_references.py:70`, `verify_check_reference_compatibility.py:84`,
  `verify_left_align.py:147,271`, `verify_validate_variants.py:115,191`,
  `verify_variant_recalibrator.py` etc.). Without the pinned jar or JDK the test still passes, and
  its (a)-classified literals silently degrade into unverified native expectations — exactly the
  `verify_indel.py:173` failure mode, but undetectable because the test stays green.
  `verify_fixture_digests.py` does not help: it hashes only *inputs*.
* Many `-contract` tests invoke the jar only as a **reader** (`SelectVariants`), never as a
  comparator: `verify_genotype_gvcf.py:104`, `verify_reblock_gvcf.py:104`,
  `verify_combine_gvcfs.py:120`, `verify_select_variants.py`, `verify_gather_vcfs.py`,
  `verify_variant_recalibrator.py`, `verify_apply_vqsr.py`, `verify_filter_mutect_calls.py`,
  `verify_mutect2.py`, `verify_validate_variants.py`. A `SelectVariants` read-back proves the file
  is *parseable*, not that any value in it is right.
* `verify_remote_staging.py:127,157` pins `["chr1\t1\t8\t0.500000"]` with **no GATK invocation in
  the file at all** (the value is currently correct — confirmed equal to GATK's AnnotateIntervals
  output — but nothing in the test enforces that).

### 4.4 Confirmed class-(c)/(b) divergences in the variant-manipulation family

These were measured by running the pinned jar on the scripts' own fixtures and diffing against the
native binary (all commands capped; scratch in a temp dir). Literal output as reported:

**(1) `verify_select_variants.py:820` — `assert ref_only_records == []` (class (c), uncatalogued).**
Fixture `:801-812`, `--remove-unused-alternates` leaving an all-hom-ref record:

```
NATIVE: []
GATK  : chr1 20 . A . 50 PASS AN=2  GT:AD:GQ:PL 0/0:30:30:0
```

GATK **keeps** the record (ALT `.`); native drops it. The in-file comment at `:798-800` asserts the
opposite of GATK's behaviour. A parity fix that emits the ref-only record fails `:820`.

**(2) `verify_variant_recalibrator.py:84,422,437,456` — `culprit=full-covariance-gmm` (class (c)).**

```
GATK   model first line: '#:GATKReport.v1.1:8'          (script :274 — identical)
GATK   'MQ  5.7666666666666664E+01'                     (uppercase E, VariantRecalibrator.java:839 "%.16E")
NATIVE 'MQ  5.7666666666666664e+01'                     (script :276-277 — lowercase e)
GATK   recal row: chr1 1 . N <VQSR> . . END=1;NEGATIVE_TRAIN_SITE;POSITIVE_TRAIN_SITE;VQSLOD=1.4860;culprit=MQ
NATIVE recal row: chr1 1 . N <VQSR> . . VQSLOD=1.4860;culprit=full-covariance-gmm;POSITIVE_TRAIN_SITE;NEGATIVE_TRAIN_SITE;END=1
```

`grep -c full-covariance-gmm` against the pinned jar → **0 occurrences**. GATK writes the *name of
the worst annotation*. Four assertions pin the native-only provenance string; the sibling oracle
test pins a third value (`serialized-gmm`).

**(3) `verify_variant_filtration.py:504-505` — `AS_FilterStatus=LowASQD,PASS` /
`PASS,LowASQD` (class (c)).**

```
GATK   (expr 'vc.getAttribute("AS_QD") < 2'): chr1 1 . A C,G 50 PASS AS_FilterStatus=SITE|SITE;AS_QD=1.0,3.0
GATK   (expr 'AS_QD < 2'):                    chr1 1 . A C,G 50 PASS AS_FilterStatus=SITE|SITE;AS_QD=1.0,3.0
NATIVE (expr 'AS_QD < 2'):                    chr1 1 . A C,G 50 PASS AS_QD=1,3;AS_FilterStatus=LowASQD,PASS
```

GATK emits the `SITE|SITE` placeholder and does not allele-filter at all for either spelling.

**(4) Positional FORMAT-order assertions (class (b)→(c)).**
`verify_select_variants.py:116` (`…split(":")[-1] == "80"`) and `:790-795`;
`verify_left_align.py:111-113` (`split_gq["C"]=="30"`, `split_gq["G"]=="60"`) and `:326-327`.
The *values* match GATK; only the positions differ, because native preserves input FORMAT order
while GATK always rewrites:

```
GATK  : chr1 10 . A C 70 PASS AC=0;AF=0.00;AN=2  GT:AD:GQ:PL 0/0:8,2:30:30,0,40
NATIVE: chr1 10 . A C 70 PASS AC=0;AN=2;AF=0      GT:AD:PL:GQ 0/0:8,2:30,0,40:30
```

Any fix adopting GATK's `GT:AD:GQ:PL` order silently breaks all four assertions.
(Incidentally native emits `AF=0` where GATK emits `AF=0.00` on this path — the catalogued
AF-zero rule is apparently not applied here; not asserted by the test.)

**(5) `verify_select_variants.py:420` and `:412,413,416,423` (class (b)/(c)).**

```
DIFF 'vc.getType() == VariantContext.Type.SNP'  gatk=[]  script=['1','2']   (script :420, GATK rc=0)
ERR  'vc.isMultiallelic()'       rc=2  A USER ERROR: Invalid JEXL expression detected for select-0  (:412)
ERR  'vc.isTransition()'         rc=2   (:413)
ERR  'vc.isPass()'               rc=2   (:416)
ERR  'vc.hasAlternateAllele(1)'  rc=2   (:423)
```

Native's JEXL is a superset of GATK's; restricting it to GATK's grammar (the parity direction)
breaks `:420` and four others.

**(6) `verify_mark_duplicates.py:191,193` (class (b)).** Native histogram rows `["1","0","0","1"]` /
`["2","1","1","0"]`; GATK writes `1.0\t0\t0\t1` / `2.0\t1\t1\t0` (Picard `set_size` is a Double).
The `## HISTOGRAM\tjava.lang.Double` header (`:190`) is identical, so only the value rendering differs.

**(7) `verify_gather_vcfs.py:125-126,130` (class (b)).** GATK writes exactly **one**
`##GatherVcfs.comment=hello` line for `-CO hello --COMMENT world`, `-CO hello -CO world`, and
`--COMMENT hello --COMMENT world` alike; native writes two and asserts `comments == 2`.

**(8) `verify_variant_filtration.py:77` (class (b)).**

```
GATK 4-filter rc=3  java.lang.NumberFormatException: For input string: "1.0"
GATK 3-filter FILTERs: ['Het;LowDP;LowQual', 'PASS', 'Het']
NATIVE 4-filter FILTERs: ['Het;LowDP;LowQD;LowQual', 'PASS', 'Het;LowQD']
```

The ordering is the same lexicographically, but GATK cannot execute the 4-filter command line at
all — so the `LowQD` term is native-only and the assertion is unoracleable as written.

**Confirmed as genuine (a) in the same batch — do not relax:** `verify_bqsr.py:122`
(`#:GATKReport.v1.1:5`) and `:126` (94 quantized rows); `verify_sort_sam.py:55`
(`@HD VN:1.6 SO:coordinate`); `verify_apply_vqsr.py:95` (tranche `+` suffix) and `:348-349`
(`NA` = GATK's `emptyStringValue`, `ApplyVQSR.java:197`); `verify_variants_to_table.py:305-306`,
`:429-431`, `:264`; `verify_variant_recalibrator.py:274`, `:85`, `:86-88`, `:372`
(`gatk_random_seed == 47382911`, `Utils.java:52`); `verify_vqsr_scatter_joint.py:99-111`
(tranche `Version number 6`/`5`); and the great majority of native's JEXL selections in
`verify_select_variants.py:403-458,584-625`.

---

## 5. Ranked list — tests most likely to BLOCK a future GATK-parity fix

| # | test / assertion | what it pins | class + evidence | cheapest confirmation |
|---|---|---|---|---|
| 1 | `verify_reblock_gvcf.py:64` `assert len(records) == 3` (+`:82-83` record 3 is a `<NON_REF>` block with `END=20`) | `--drop-low-quals` converts a low-PL[0] site to a reference block instead of dropping it | **(c)**, **measured** (§3.2): GATK emits 2 records; `ReblockGVCF.java:542` returns null | re-run §3.2 command pair (~4 s) |
| 2 | `verify_genotype_gvcf.py:391` `assert len(star_records) == 1 and "*" in star_records[0][4]` | native emits a `*` + concrete-ALT record GATK never emits | **(c)**, **measured** (§3.3): GATK default 0 records, `-all-sites` gives `ALT=.`/`GT=./.` | re-run §3.3 (~4 s) |
| 3 | `verify_reblock_gvcf.py:283-287` triploid `ALT=="C,<NON_REF>"`, `GT=="0/1/1"`, `PL=="0,10,20,30,100,110,120,160,170,190"` | native keeps a variant whose min-PL genotype is hom-ref; GATK reblocks it to `0/0/0` | **(c)**, **measured** (§3.2); mechanism `ReblockGVCF.java:515-531` | re-run §3.2 triploid pair (~4 s) |
| 4 | `verify_filter_mutect_calls.py:156-157` `base_qual`/`map_qual` (and `:396,:398,:506`, `:216,:237`) | MBQ/MMQ read as `Number=1` still fire allele-specific filters | **(c)**, measured by delegated audit | `FilterMutectCalls` with no options on `defaults_body` |
| 5 | `verify_filter_mutect_calls.py:64` `FAIL` + `:664` `low_tlod` | a native-only `--min-tlod` filter vocabulary | **(c)**, measured: `min-tlod is not a recognized option` | one `FilterMutectCalls --min-tlod` run |
| 6 | `verify_depth_of_coverage_multisample.py:52-59` per-sample `Average_Depth_S1/S2` columns and `1.00` denominators | native's per-sample averaging where GATK emits one `Average_Depth_sample` | **(c)**, measured (§4.1.3) | `SortSam` + `DepthOfCoverage` on the 4-read SAM |
| 7 | `verify_reblock_gvcf.py:77` `records[1][9].split(":")[3]` and `:80` `…[-1] == "20"` | positional dependence on native FORMAT order `GT:DP:AD:PL:GQ` | **(b), high** — GATK's order is `GT:AD:DP:GQ:PL` (measured, §3.2); adopting GATK's order breaks both lines | compare FORMAT columns in §3.2 output |
| 8 | `verify_hc_genotype_priors.py:128,130,163,181,183,184` exact `PL`/`GQ`/`GT` strings on a synthetic 2-ALT SAM | native PairHMM/prior numerics recorded as expectation; GATK is used only for `GenotypePriorCalculator` vectors, never for a record | **(b), high** — any parity fix touching assembly/PairHMM/prior normalization changes these values | run GATK `HaplotypeCaller` on the same SAM (needs a `.fai`/`.dict`; the fixture is synthetic) |
| 9 | `verify_mutect2.py:650` sample-column set (plus `:647,:678,:684`) | native drops the `OTHER_NORMAL` sample column | **(c)**, measured (§4.1.4) | reproduce `:587-640` fixture + GATK `Mutect2` |
| 10 | `verify_combine_gvcfs.py:70-71` `len(records)==2` + coalesced `END=7` + `DP=10` | native coalesces adjacent blocks on a sample-less gVCF | **(b), fixture-scoped** — measured (§3.6): GATK emits 3 records and no INFO; fixture has no genotypes so GATK cannot merge | re-run §3.6 |
| 11 | `verify_count_reads.py:143-150`, `verify_flag_stat.py:259-266` negative `-ip` rejection; `:193-212`/`:136-155`/`verify_count_bases_in_reference.py:100,108`/`verify_split_intervals.py:115-145` unindexed-gz acceptance; `verify_count_bases_in_reference.py:84-85` 3-field interval_list | native input-validation envelope GATK does not share | **(c)**, measured (§4.1.5-7) | one command per case; all documented in §4.1 |
| 12 | `verify_genotype_gvcf.py:353-367` multi-sample GT/GQ/AC/AN/AF | joint PL→GT materialization | **(b)** — cannot be checked against GATK until `GQ` is added to the fixture header (§3.5) | add the header line, re-run `GenotypeGVCFs` |
| 13 | `verify_reblock_gvcf.py:318-332` multi-sample block merge (`sample_count==2`, per-sample MIN_DP/DP/GQ) | a native extension GATK refuses outright | **(b)** — measured: GATK `A USER ERROR … Found samples [S1, S2]` | re-run §3.2 M fixture |
| 14 | `verify_hc_streaming.py:65,67`, `verify_hc_region_streaming.py:50,162`, `verify_hc_rcm_realignment.py:51` | byte-equality between native's own modes and native's own default-vs-flag output; record counts | **(b)** — native-only flags (`--stream-by-contig`) or documented divergence (`--stream-by-region`) | already documented in `:67-85` of the region-streaming script |
| 15 | `verify_genotype_gvcf.py:87,995` `RCQ=`/`RCP=` presence, `:232-234` tag removal, `:762` `GQ=="60"`, `:907,:951` record counts | ReferenceConfidenceModel annotations | **(b)** — no GATK comparison; RCQ/RCP are GATK annotations but the values/counts are native's | run GATK `GenotypeGVCFs` on the script's HC gVCF |
| 16 | **`verify_select_variants.py:820`** `assert ref_only_records == []` | native **drops** the all-hom-ref record that `--remove-unused-alternates` leaves; GATK keeps it as `ALT=.` — **uncatalogued divergence** | **(c)**, measured (§4.4.1) | 6-line fixture + one native and one GATK run |
| 17 | **`verify_variant_recalibrator.py:84,422,437,456`** `culprit=full-covariance-gmm` / `AS_culprit=full-covariance-gmm` | a provenance string that exists nowhere in GATK; GATK writes the worst annotation name (`culprit=MQ`) | **(c)**, measured (§4.4.2): `grep -c full-covariance-gmm <jar>` → 0 | one `VariantRecalibrator` run, or the grep |
| 18 | **`verify_variant_filtration.py:504-505`** `AS_FilterStatus=LowASQD,PASS` / `PASS,LowASQD` | native's allele-level filter encoding; GATK emits `SITE|SITE` and applies no allele filter | **(c)**, measured (§4.4.3) | one `VariantFiltration --apply-allele-specific-filters` run |
| 19 | Positional FORMAT reads: `verify_select_variants.py:116,790-795`; `verify_left_align.py:111-113,326-327` | native's preserved input FORMAT order; GATK rewrites to `GT:AD:GQ:PL` | **(b)→(c)**, measured (§4.4.4) — values identical, positions not | diff FORMAT columns from §4.4.4 |
| 20 | `verify_select_variants.py:420` `vc.getType()==VariantContext.Type.SNP`; `:412,:413,:416,:423` | native JEXL is a superset of GATK's; GATK selects `[]` or rejects (rc 2) | **(b)/(c)**, measured (§4.4.5) | re-run those `--select` expressions under GATK |
| 21 | `verify_variant_recalibrator.py:276-277` lowercase `e+01` vs GATK `%.16E`; `verify_mark_duplicates.py:191,193` `1`/`2` vs `1.0`/`2.0`; `verify_gather_vcfs.py:125-126,130` two comment lines vs one | exact numeric/line rendering | **(b)**, all measured (§4.4.2, §4.4.6, §4.4.7) | diff the corresponding GATK outputs |
| 22 | `verify_variant_filtration.py:77` `"Het;LowDP;LowQD;LowQual"` | a 4-filter command line GATK **cannot execute** (`NumberFormatException: "1.0"`, rc 3); the in-file oracle also sorts filter names and drops INFO, so ordering is unverified | **(b)** — measured (§4.4.8) + comparator reading | one GATK `VariantFiltration` run |
| 23 | remaining count/ratio literals: `verify_collect_read_counts.py:64-65`, `verify_collect_allelic_counts.py:79-81`, `verify_denoise_read_counts.py:113-116,164`, `verify_call_copy_ratio_segments.py:55`, `verify_gather_vcfs.py:70,188`, `verify_apply_vqsr.py:93-96,256-258`, `verify_validate_variants.py:52,71`, `verify_variant_eval.py:100-112` | native pileup/annotation/report arithmetic recorded as expectation, no comparator for those values | **(b)** | run the matching GATK tool on the same fixture |

Notes on rank: rows 4-6 are the ones a *somatic/coverage* parity pass hits first; rows 1-3, 16-18 are
the ones a *gVCF* parity pass hits first. Row 8 (`verify_hc_genotype_priors.py`) is the most fragile
non-measured item: it pins exact PairHMM likelihood strings on a synthetic fixture and nothing in
the test can tell you whether they are still right.

---

## 6. Recommended remediation order

1. **Fix the `ReblockGVCF` / `GenotypeGVCFs` contracts first** (ranked 1-3, 6, 7). They sit in
   the same subsystem as the original incident, they are cheap to re-verify, and two of them
   (`:64`, `:391`) assert a *record that GATK does not emit*. Either (a) make them real oracles by
   running the pinned jar on the same fixture and diffing, or (b) correct the expectations and
   record the divergence in the project catalogue.
2. **Then the six variant-manipulation (c) items** (ranked 16-18 and the two measured (b)s at 19-21).
   All of them were settled with a single short GATK run each; three of them
   (`verify_select_variants.py:820`, `verify_variant_recalibrator.py:84`, `verify_variant_filtration.py:504-505`)
   pin values GATK never produces, so they are the cheapest places to stop locking in divergence.
3. **Convert "jar-as-reader" tests into comparators.** For every `-contract` script that calls
   `SelectVariants`/`IndexFeatureFile` and then asserts literals (`verify_genotype_gvcf.py`,
   `verify_reblock_gvcf.py`, `verify_combine_gvcfs.py`, `verify_select_variants.py`,
   `verify_gather_vcfs.py`, `verify_variant_recalibrator.py`, `verify_apply_vqsr.py`,
   `verify_filter_mutect_calls.py`, `verify_mutect2.py`), add the matching GATK tool run and diff
   the data rows. This is the change that would have caught `verify_indel.py:173` automatically.
4. **Make the oracle guard fail closed.** Replace `if java.exists() and jar.exists():` with a hard
   requirement (or honour the existing `FASTGATK_REQUIRE_GATK_ORACLE=1` convention used by
   `verify_hc_likelihood_filter.py` and by the second registration of `verify_analyze_covariates.py`)
   in the 11+ scripts in §4.3, and add CI enforcement so a missing pinned jar is a failure, not a
   silent downgrade to native-only expectations. 13 of the 18 variant-manipulation scripts are
   "dual-mode" — they run the jar *somewhere* while keeping unprotected native-only assertions.
5. **De-positionalise FORMAT assertions.** Replace `record[9].split(":")[...]` and
   `split(":")[-1]` arithmetic with name-keyed lookups (`dict(zip(fmt, values))`) wherever the field
   is asserted (`verify_reblock_gvcf.py:77,80`, `verify_genotype_gvcf.py:362-363`,
   `verify_select_variants.py:116,790-795`, `verify_left_align.py:111-113,326-327`,
   `verify_combine_gvcfs.py`, `verify_genomicsdb_import.py:372-374`). This removes a whole class of
   latent breakage if native ever adopts GATK's FORMAT ordering.
6. **Quarantine the input-envelope divergences** (§4.1.5-7) into explicitly-named
   `*-native-extension-*` tests, so that a future "match GATK's CLI validation" pass can delete
   them without ambiguity.
7. **Record the new (c) items in the project's divergence catalogue** so the next audit does not
   have to rediscover them: ReblockGVCF `--drop-low-quals` drop-vs-block; ReblockGVCF
   min-PL-genotype reblocking; GenotypeGVCFs spanning-deletion emission; native multi-sample
   ReblockGVCF; FilterMutectCalls MBQ/MMQ arity + `--min-tlod`/`low_tlod`/`FAIL`;
   DepthOfCoverage per-sample columns; Mutect2 sample-column set; negative `-ip`/`-ixp`;
   unindexed-gzip interval inputs; 3-field `.interval_list`;
   `SelectVariants --remove-unused-alternates` ref-only record drop;
   `VariantRecalibrator culprit=full-covariance-gmm`; `VariantFiltration` allele-filter
   `AS_FilterStatus=LowASQD,…`; native JEXL superset; FORMAT-order preservation.

---

## 7. Confirmations performed vs reasoned only

**Confirmed by execution (literal output in §3 and §4.1):**

* `verify_indel.py` — native vs GATK gVCF on the test's own fixture: **byte-identical**, all three
  records (§3.1).
* `verify_reblock_gvcf.py` — fixture A, deletion fixture, triploid fixture, multi-sample fixture,
  all run through both `fastgatk-reblock-gvcf` and GATK `ReblockGVCF` with identical flags; plus
  GATK source reading of `lowQualVariantToGQ0HomRef`, `changeCallToHomRefVersusNonRef`,
  `shouldBeReblocked`, `addRefBlockIfNecessary`, `updateMQAnnotations` (§3.2).
* `verify_genotype_gvcf.py` — the exact `:378-381` spanning-deletion fixture (3 PL variants) and
  the `:329-343` multi-sample fixture through `fastgatk-genotype-gvcf` and GATK `GenotypeGVCFs`
  default / `-all-sites` (§3.3, §3.5).
* `verify_hc_variant_annotations.py`, `verify_hc_likelihood_filter.py`,
  `verify_hc_chr17_independent_regions.py` — HC runs on the bundled chr17 fixture at the default
  and `--dont-use-soft-clipped-bases` thresholds (§3.4).
* `verify_combine_gvcfs.py` — shard-merge, `--break-bands-at-multiples-of 5`, and
  `--convert-to-base-pair-resolution` fixtures through both binaries (§3.6).
* Delegated: `verify_count_reads.py`, `verify_flag_stat.py`, `verify_split_intervals.py`,
  `verify_annotate_intervals.py`, `verify_count_bases_in_reference.py` (GATK executions, §4.1.5-7);
  `verify_filter_mutect_calls.py`, `verify_depth_of_coverage_multisample.py`,
  `verify_mutect2.py` sample columns, `verify_get_pileup_summaries.py`/`verify_calculate_contamination.py`
  (GATK executions, §4.1.1-4); GATK-source readings of `GATKVCFConstants.java`,
  `BaseQualityFilter.java`, `MappingQualityFilter.java`, `CoverageOutputWriter.java`,
  `GatherPileupSummaries.java`, `PileupSummary.java`, `TableWriter.java`.
* Delegated (variant manipulation, §4.4): GATK executions on the scripts' own fixtures for
  `verify_select_variants.py` (main fixture, ref-only record, block `-L`, ~60 JEXL selections),
  `verify_variant_recalibrator.py` (model text + recal rows, plus `grep -c full-covariance-gmm` on
  the jar → 0), `verify_variant_filtration.py` (AS filters + the 4-filter/3-filter crash),
  `verify_left_align.py` (split fixture), `verify_mark_duplicates.py` (paired fixture metrics),
  `verify_gather_vcfs.py` (comment header combinations), `verify_bqsr.py`
  (`#:GATKReport.v1.1:5`, 94 quantized rows), `verify_sort_sam.py` (`@HD`), and
  `verify_variants_to_table.py` (three TSV comparisons).

**Reasoned only (not executed):**

* Classification of the class-(b) literal sets in `verify_collect_allelic_counts.py`,
  `verify_collect_read_counts.py`, `verify_denoise_read_counts.py`,
  `verify_call_copy_ratio_segments.py`, `verify_create_read_count_panel_of_normals.py`,
  `verify_model_segments.py` (extension parts), `verify_genomicsdb_import.py`,
  `verify_select_variants.py:112-124`, `verify_gather_vcfs.py:70,188`, `verify_apply_vqsr.py:93-96`,
  `verify_hc_genotype_priors.py` — these rest on reading the scripts (no comparator for those
  values in-file), not on a GATK run.
* All manifest/telemetry literals (Kokkos execution space, kernel names/lifetimes,
  `pipeline_lifecycle` strings, batch/thread counters, `HDF5-SVD-ReadCountPanelOfNormals-v7`):
  class (b) by construction; none was executed.
* `verify_mutect2.py:647,678,684` (`ALT=="A,T"`, `GT=="0/1/2"`, `DP=="0"`): the delegated GATK
  reproduction produced no data records for that fixture, so these are unresolved.
* `verify_validate_variants.py:52,71` (native report TSV) and `verify_variant_eval.py:100-112`
  (native `##table=` schema): the "GATK has no such artifact" claim is a source/grep reading
  (`ValidateVariants.java` declares no OUTPUT argument), not a GATK execution.
* `verify_sort_sam.py:98-100`, `verify_mark_duplicates.py:144-146` (`@PG` rendering),
  `verify_bqsr.py:763,783`, `verify_gather_vcfs.py:206,219,237`, `verify_optional_boolean_contract.py:77`,
  `verify_apply_vqsr.py:410`, `verify_variant_recalibrator.py:481,501`, and all native manifest
  telemetry: native-only by construction; not executed.
* `verify_hc_multi_input.py`, `verify_flow_hmer_option.py`, `verify_real_assembly_graph.py`,
  `verify_vqsr_scatter_joint.py` (`:99-111` tranche headers were read from `VQSLODTranche`/
  `TruthSensitivityTranche` source, not run), `verify_activity_region_controls.py`,
  `verify_overlapping_quality_correction.py`, `verify_analyze_covariates.py` telemetry:
  classified by reading only.

**Known gaps in this audit:** `verify_hc_genotype_priors.py`'s exact PL/GQ literals were not
re-run against GATK (the fixture is a synthetic SAM that needs its own reference/dictionary);
`verify_flow_hmer_option.py` and `verify_hc_multi_input.py` were read but not executed against GATK.

**Environment limitation:** `python3 -c "import gzip"`-written `.g.vcf.gz` fixtures cannot be
indexed by htsjdk (`MalformedFeatureFile: Input file is not in valid block compressed format`), so
all GATK-side replays of gVCF fixtures were performed on the identical content written as plain
`.vcf` + Tribble `.idx`. Native was additionally run in its normal `.gz` mode and produced the same
records, so the substitution does not affect any classification above.
