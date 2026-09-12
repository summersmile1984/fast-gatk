# Round: the unchecked htslib read idiom (`bcf_read`/`sam_read1`/`hts_getline`)

Scope: the side finding of `.diag/round-silent-drop.md` §6 — the
`while (bcf_read(...) == 0)` idiom that made GenotypeGVCFs silently truncate its
input "still exists in about 19 other `fastgatk-native/src/*.cpp` files".  That
number was a grep count, never measured.  This round censuses the whole class,
proves it on a representative subset, fixes it, and gates it.

Probes: `.diag/audit-work/probe.py` (per-tool, `.diag/audit-work/probe-before.json`
/ `probe-after.json`), `.diag/audit-work/truncgz2.py` (truncated `.vcf.gz`),
`.diag/audit-work/bamtrunc.sh` / `bamflip.py` / `bammulti.py` (corrupt BAM).
New gate: `fastgatk-native/scripts/verify_malformed_input_fail_loud_oracle.py`
(`.diag/audit-work/gate-before.log` exit 1, `.diag/audit-work/gate-after.log` exit 0).
New helper: `fastgatk-core/include/fastgatk/io/hts_read_guard.hpp`.

## 0. Verdict

**(i) The class is real, systemic, and larger than the report: 51 sites in 27
first-party files, not "about 19 files".**  Every one of them ended a traversal
at the first unreadable record as if the input had ended.  For 50 of the 51 that
meant silently discarded records with exit 0; the one exception
(`gather_vcfs_tool.cpp:317`) is flagged in §1a and had an independently guarded
sibling read that made the loss unreachable in practice.

Measured, on the *pre-fix* binaries, per tool (`bad-first` = an unreadable record
at chr1:10 followed by two well-formed records at chr1:20 and chr1:30):

| tool | control exit/rows | bad-first exit/rows | verdict |
| --- | --- | --- | --- |
| variants-to-table | 0 / 2 | **0 / 0** | both well-formed records lost |
| select-variants | 0 / 2 | **0 / 0** | lost |
| left-align-trim | 0 / 2 | **0 / 0** | lost |
| variant-filtration | 0 / 2 | **0 / 0** | lost |
| filter-mutect-calls | 0 / 2 | **0 / 0** | lost |
| reblock-gvcf | 0 / 2 | **0 / 0** | lost |
| combine-gvcfs | 0 / 2 | **0 / 0** | lost |
| validate-variants | 0 / 2 | **0 / 2** | worse: reports a *valid* file after skipping the unreadable record |
| gather-vcfs | 0 / 2 | 2 / 0 | already guarded (class B) — the fix's control |

And the realistic trigger, a **truncated `.vcf.gz`** (an interrupted pipeline
write, a partial copy) — 3000 well-formed records, cut at 90 %:

| tool | pre-fix | post-fix |
| --- | --- | --- |
| variants-to-table | **exit 0, 2057 / 3000 records** | exit 2, `BAD_INPUT: malformed or truncated VCF/BCF record (htslib status -2)` |
| select-variants | **exit 0, 2057 / 3000** | exit 2, same |
| left-align-trim / variant-filtration / gather-vcfs / filter-mutect-calls | **exit 0, 2057 / 3000** | exit 2, same |
| reblock-gvcf / combine-gvcfs | **exit 0, 2664 / 3000** | exit 2, same |
| validate-variants | **exit 0, report on 2 records** | exit 2, same |

htslib printed `[E::bgzf_read_block] Reading GZIP stream failed at offset …` to
stderr in every one of those runs, and native exited **0**.  A downstream
pipeline has no way to tell such a run from a complete one.  This is the same
defect as the previous round's, with a much more ordinary trigger.

**(ii) Two bonus findings, both measured, both fixed by the same guard.**

* A corrupt BAM used to make `sort-sam`, `mark-duplicates` and
  `collect-f1r2-counts` **abort on heap corruption** (`free(): invalid pointer`,
  exit 134) or **SIGSEGV** (exit 139) — they failed loudly but with a crash and
  no diagnostic.  Post-fix they exit **2** with
  `BAD_INPUT: malformed or truncated BAM/CRAM record (htslib status -2): <path>`.
  The crash came from the code that ran *after* the unguarded loop had exited on
  a half-populated record; throwing at the read removes it.
* `validate-variants` on a malformed input reported success.  It validates
  *records*, so a record it never saw is a record it never rejected.

**(iii) Class size detail.** 51 sites: 34 `bcf_read`, 2 `sam_read1`/`sam_itr_next`
in `fastgatk-core`, and 15 `hts_getline` (same two-negative contract, documented
in `htslib/hts.h`).  12 further sites were already correct (class B) and 5 are
not read loops at all (class C) — the census table in §1 lists every one.

## 1. STEP 1 — the census

Method: `grep` for `bcf_read|vcf_read|bcf_read1|bcf_readrec|sam_read1|sam_itr_next|
bam_read1|hts_getline` over `fastgatk-native/src`, `fastgatk-native/include`,
`fastgatk-core/src` and `fastgatk-core/include` (the only first-party trees;
`third_party/`, `gatk-source/` and `gatk-rs-source/` are vendored and excluded),
then **read every hit** together with its loop tail to decide whether a
`status < -1` check follows.

Classes: **A** can silently truncate (no `-2` guard), **B** already guarded,
**C** not a read loop / safe by construction, **D** reads a different stream where
the idiom is correct.

Line numbers are those of the **fixed tree** (the census line numbers differ by at
most +3 because the guard changes are one line each and the include adds one).

### 1a. Class A — `bcf_read` family (36 sites, all fixed)

| file:line | tool / enclosing function | idiom | class | reason |
| --- | --- | --- | --- | --- |
| `apply_vqsr_tool.cpp:425` | ApplyVQSR, `read_recal_vcf` | `while (bcf_read(input, header, record) == 0)` | A | the recal VCF is a *scored* input: losing its tail silently drops VQSLOD assignments |
| `apply_vqsr_tool.cpp:780` | ApplyVQSR, main loop | `while (bcf_read(input, input_header, record) == 0)` | A | primary variant loop |
| `bqsr_tool.cpp:3266` | ApplyBQSR, `decode` lambda | `while (batch.records.size() < N && sam_read1(...) >= 0)` | A | a BAM read failure ends the decode batch as if the BAM had ended |
| `collect_f1r2_counts_tool.cpp:494` | CollectF1R2Counts, unindexed branch | `while (sam_read1(...) >= 0) process_record();` | A | **the indexed branch 5 lines above *is* guarded** (`collect_f1r2_counts_tool.cpp:490`); the report's claim that `:489` covered both branches was wrong |
| `combine_gvcf_tool.cpp:1158` | CombineGVCFs, `read_next_raw` lambda | `return bcf_read(...) == 0;` | A | its caller (`:1290`) maps `false` to `cursor.eof = true`; the lambda cannot report an error |
| `combine_gvcf_tool.cpp:1709` | CombineGVCFs, `run_tool` | `while (bcf_read(input, header, record) == 0)` | A | non-streaming record read |
| `fasta_reference_tool.cpp:433` | FastaAlternateReferenceMaker, `read_variants` | `while (record && bcf_read(file, header, record) == 0)` | A | a lost tail writes the reference base instead of the ALT |
| `filter_mutect_tool.cpp:4003` | FilterMutectCalls, haplotype pass | `while (bcf_read(input, header, record) == 0)` | A | |
| `filter_mutect_tool.cpp:4669` | FilterMutectCalls, orientation learning | `while (bcf_read(input, input_header, record) == 0)` | A | a truncated first pass silently learns a threshold from a partial population |
| `filter_mutect_tool.cpp:4726` | FilterMutectCalls, haplotype learning | `while (bcf_read(input, input_header, record) == 0)` | A | same |
| `filter_mutect_tool.cpp:4819` | FilterMutectCalls, ErrorProbabilities learning | `while (bcf_read(input, input_header, record) == 0)` | A | same |
| `filter_mutect_tool.cpp:5093` | FilterMutectCalls, main filtering loop | `while (bcf_read(input, input_header, record) == 0)` | A | |
| `gather_vcfs_tool.cpp:317` | GatherVcfs, `reorder_inputs_by_first_variant` | `if (bcf_read(file, header, record) == 0)` | A\* | not independently reachable: the copy loop at `:630` is guarded and re-reads the same file, so a `-2` here can never be the *only* signal.  Fixed anyway — treating a parse failure as "empty shard" is wrong on its face |
| `get_pileup_summaries_tool.cpp:482` | GetPileupSummaries, `load_interval_variant_keys` (`-L` VCF) | `while (bcf_read(file, header, record) == 0)` | A | a truncated `-L` VCF silently narrows the site set |
| `get_pileup_summaries_tool.cpp:842` | GetPileupSummaries, `load_sites`, unindexed branch | `while (bcf_read(file, header, record) == 0)` | A | **the indexed branch of the same `if/else` is guarded at `:834`** |
| `get_pileup_summaries_tool.cpp:934` | GetPileupSummaries, `load_intervals_for_sites` | `while (bcf_read(...) == 0)` | A | |
| `left_align_tool.cpp:1005` | LeftAlignAndTrimVariants | `while (bcf_read(input, input_header, record) == 0)` | A | |
| `mark_duplicates_tool.cpp:1113` | MarkDuplicates, first pass | `while (sam_read1(input, header, record) >= 0)` | A | also crashed (exit 134) on a corrupt BAM pre-fix |
| `mark_duplicates_tool.cpp:1350` | MarkDuplicates, second pass | `while (sam_read1(second_input, second_header, record) >= 0)` | A | |
| `mutect2_tool.cpp:320` | Mutect2, `load_mutect_feature_table` | `while (bcf_read(input, header, record) == 0)` | A | a partial `--germline-resource`/`--alleles` table shifts the AF prior silently |
| `reblock_gvcf_tool.cpp:1847` | ReblockGVCF | `while (bcf_read(input, header, record) == 0)` | A | |
| `select_variants_tool.cpp:541` | SelectVariants, `load_comparison_data` | `while (bcf_read(input, header, record) == 0)` | A | a short `--concordance`/`--discordance` file inverts the selection |
| `select_variants_tool.cpp:2358` | SelectVariants, main loop | `while (bcf_read(input, input_header, record) == 0)` | A | |
| `sort_sam_tool.cpp:755` | SortSam, input pass | `while (sam_read1(input, header, record) >= 0)` | A | also crashed (exit 134/139) on a corrupt BAM pre-fix |
| `sort_sam_tool.cpp:838` | SortSam, spill-merge heap | `const int status = sam_read1(...)` | A | a failed re-read of a spill run silently shortens the merge |
| `validate_variants_tool.cpp:556` | ValidateVariants, `load_dbsnp_ids` | `while (bcf_read(file, header, record) >= 0)` | A | a short dbSNP file turns "ID not in dbSNP" errors into passes |
| `validate_variants_tool.cpp:777` | ValidateVariants, main loop | `while (bcf_read(input, header, record) >= 0)` | A | **worst variant of the class**: exit 0 with a *pass* verdict after never examining the record |
| `variant_filtration_tool.cpp:2152` | VariantFiltration, mask VCF | `while (bcf_read(input, header, record) == 0)` | A | a truncated `--mask` silently unmasks the tail |
| `variant_filtration_tool.cpp:2209` | VariantFiltration, clustering VCF | `while (bcf_read(input, header, record) == 0)` | A | |
| `variant_filtration_tool.cpp:2521` | VariantFiltration, main loop | `while (bcf_read(input, input_header, record) == 0)` | A | |
| `variant_recalibrator_tool.cpp:652` | VariantRecalibrator, `read_resource_keys` | `while (bcf_read(input, header, record) == 0)` | A | a partial resource silently relabels training/truth membership |
| `variant_recalibrator_tool.cpp:2396` | VariantRecalibrator, `read_entries` | `while (bcf_read(input, header, record) == 0)` | A | a partial training set silently changes the GMM |
| `variant_recalibrator_tool.cpp:2709` | VariantRecalibrator, `write_recal_vcf` | `while (bcf_read(input, input_header, record) == 0)` | A | |
| `variants_to_table_tool.cpp:1117` | VariantsToTable | `while (bcf_read(input, header, record) == 0)` | A | |
| `hts_reader.cpp:470` | `HtsReader::Impl::next_record`, unindexed | `const auto status = sam_read1(...); return status >= 0;` | A | the shared Host reader used by the calling tools |
| `hts_reader.cpp:486` | `HtsReader::Impl::next_record`, indexed | `const auto status = sam_itr_next(...)` | A | `-1` advances to the next interval, `< -1` is an error — both were swallowed |

### 1b. Class A — `hts_getline` (15 sites, all fixed)

`hts_getline` documents `-1 on end-of-file; <= -2 on error`
(`third_party/htslib-build/htslib-src/htslib/hts.h:680-681`), and htslib's own
caller honours it (`hts.c:2060` `while ((ret = bgzf_getline(...)) >= 0)` …
`hts.c:2072` `if (ret < -1) // Read error`).  For a BGZF/gzip stream `-2` is what a
truncated file produces, so these share the defect exactly; for a plain text file
`-2` needs a real I/O error (`hts.c:2031`).

| file:line | tool / function | read |
| --- | --- | --- |
| `check_reference_compatibility_tool.cpp:100` | CheckReferenceCompatibility, `parse_vcf_text_contigs` | a text VCF `.dict` source |
| `apply_vqsr_tool.cpp:294` | ApplyVQSR, `for_each_tranche_line` | tranches CSV |
| `preprocess_intervals_tool.cpp:182` | PreprocessIntervals, `read_lines` | interval list |
| `variants_to_table_tool.cpp:276` | VariantsToTable, `RawInfoReader` | **the raw AF/INFO re-read of the input VCF** — a truncated `.vcf.gz` here loses the high-precision INFO text |
| `split_intervals_tool.cpp:324` | SplitIntervals, compressed interval branch | interval list |
| `count_bases_in_reference_tool.cpp:183` | CountBasesInReference, compressed interval branch | interval list |
| `get_pileup_summaries_tool.cpp:571` | GetPileupSummaries, `load_raw_vcf_af` | **the population VCF re-read that preserves `AF` digits** |
| `fasta_reference_tool.cpp:175` | FastaReferenceMaker/Alterer, `read_text_lines` | `.dict` / FASTA index text |
| `genotype_gvcf_tool.cpp:4670` | GenotypeGVCFs, `input_pass_filter_line` | input VCF header scan |
| `variant_recalibrator_tool.cpp:895` | VariantRecalibrator, `for_each_report_line` | GATKReport model text |
| `filter_intervals_tool.cpp:296` | FilterIntervals, `read_lines` | interval list |
| `read_metrics_tool.cpp:779` | CollectAlignmentSummaryMetrics, compressed interval branch | interval list |
| `annotate_intervals_tool.cpp:215` | AnnotateIntervals, `read_lines` | interval list |
| `hts_reader.cpp:559` | `HtsReader::Impl::append_file` | interval list |
| `intervals.hpp:97` | `fastgatk::io::read_interval_lines` (shared by every interval-taking tool) | interval list |

### 1c. Class B — already guarded (12 sites, **not** touched)

| file:line | guard, quoted |
| --- | --- |
| `bqsr_tool.cpp:589` | `const auto status = bcf_read(...); if (status < 0) { if (status != -1) { … throw "BAD_INPUT: failed reading known-sites VCF/BCF" } break; }` |
| `collect_f1r2_counts_tool.cpp:487` | `while ((status = sam_itr_next(...)) >= 0) …; if (status < -1) throw "BAD_INPUT: indexed BAM/CRAM iterator failed"` |
| `gather_vcfs_tool.cpp:629/670` | `while ((read_status = bcf_read(...)) == 0) { … } if (read_status < -1) throw "BAD_INPUT: GatherVcfs failed while reading shard"` |
| `genomicsdb_import_tool.cpp:841/864` | `if (status < -1) { … throw "BAD_INPUT: failed scanning native GenomicsDB input" }` |
| `genomicsdb_import_tool.cpp:1192/1211` | `if (status < -1 \|\| output_status != 0 \|\| input_status != 0) { … throw … }` |
| `genotype_gvcf_tool.cpp:4927` | previous round: `if (status < -1) throw "BAD_INPUT: streaming GenotypeGVCFs VCF record parse failed"` |
| `genotype_gvcf_tool.cpp:5605/5635` | previous round: `if (probe_status < -1) { … throw … }` |
| `genotype_gvcf_tool.cpp:6535/6540` | previous round: `if (read_status < -1) { … throw "BAD_INPUT: GenotypeGVCFs input record parse failed" }` |
| `hc_call.cpp:2252` | `const auto status = bcf_read(...); if (status < 0) { if (status != -1) throw "BAD_INPUT: cannot read --alleles VCF/BCF record"; break; }` |
| `read_metrics_tool.cpp:1276` | `status = sam_itr_next(...); … if (status < -1) throw "BAD_INPUT: HTSlib failed while reading indexed interval"` |
| `read_metrics_tool.cpp:1287` | `status = sam_read1(...); … if (status < -1) throw "BAD_INPUT: HTSlib failed while reading records"` |
| `variant_eval_tool.cpp:893/901` | `while ((read_status = bcf_read(file, header, raw)) >= 0) { … } if (read_status < -1) throw "BAD_INPUT: failed reading variant input"` |

These prove the idiom was already the house convention — the previous round's
GenotypeGVCFs fix copied it, and this round applies it to the rest.

### 1d. Class C — not a read loop, or safe by construction (5 sites, not touched)

| file:line | why it is safe |
| --- | --- |
| `bqsr_tool.cpp:3217` | `if (sam_read1(prefix_input, prefix_header, prefix_record) < 0 \|\| sam_write1(...) < 0) throw "BAD_INPUT: ApplyBQSR checkpoint output is shorter than checkpoint"` — *any* negative throws; there is no silent path |
| `sort_sam_tool.cpp:817` | `if (!state.header \|\| !state.record \|\| sam_read1(...) < 0) { … throw "cannot initialize spill merge" }` — same, and the file is a spill run this process wrote itself |
| `variants_to_table_tool.cpp:284` | a comment naming `bcf_read` for alignment; not a call |
| `get_pileup_summaries_tool.cpp:557` | a comment; not a call |
| `genotype_gvcf_tool.cpp:4914-4915` | comments; not calls |

### 1e. Class D — a different stream where the idiom is correct

The indexed iterators `tbx_itr_next` / `bcf_itr_next` at
`get_pileup_summaries_tool.cpp:812` and `:830` (guarded once for the whole
`if/else` by `if (status < -1)` at `:834`), `genotype_gvcf_tool.cpp:6501`/`:6520`
(guarded), `genotype_gvcf_tool.cpp:4938`/`:4947` (throws on `< -1`) and
`combine_gvcf_tool.cpp:1165`/`:1175` (guarded at `:1178`) are already correct:
for these, `-1` legitimately means "this iterator is exhausted, move to the next
interval", which is exactly why the check must be `< -1` and not `< 0`.

## 2. STEP 1 (proof) — what was measured and what was only inferred

**Proved by measurement (9 tools, in the gate):** variants-to-table,
select-variants, left-align-trim, variant-filtration, validate-variants,
reblock-gvcf, combine-gvcfs, filter-mutect-calls,
fasta-alternate-reference-maker — plus gather-vcfs as an already-correct control.
For each: pre-fix exit 0 with records missing (§0), post-fix exit 2 with
`BAD_INPUT: … (htslib status -2)` and the controls unchanged.
Literal pre-fix numbers are in `.diag/audit-work/probe-before.json` and
`.diag/audit-work/gate-before.log` (53 violations, exit 1).

**Proved by measurement (3 tools, BAM):** sort-sam, mark-duplicates,
collect-f1r2-counts on a 50 %-truncated BAM — pre-fix exit 134
(`free(): invalid pointer`) / 139 (SIGSEGV) with no diagnostic, post-fix exit 2
with `BAD_INPUT: malformed or truncated BAM/CRAM record (htslib status -2)`.

**Inferred, not measured (each is the same one-line idiom on the same reader):**
ApplyVQSR's recal-VCF reader and main loop, BaseRecalibrator's known-sites
reader (`bqsr_tool.cpp:589` is class B, its sibling reader at `:3266` is class A),
Mutect2's feature table, VariantRecalibrator's resource/entry/recal readers,
SelectVariants' comparison file, ValidateVariants' dbSNP loader, MarkDuplicates'
second pass, SortSam's spill merge, `HtsReader`, all 15 `hts_getline` sites, and
GatherVcfs' reorder scan.  The inference is strong — identical idiom, identical
reader, measured effect on the same reader elsewhere — but it is an inference and
is labelled as such.

## 3. The htslib rule, with file:line

Native is not entitled to invent its own contract: htslib states the rule in its
headers and follows it in its own code.

* `third_party/htslib-build/htslib-src/htslib/vcf.h:403` — `bcf_read` doc:
  `@return 0 on success; -1 on end of file; < -1 on critical error`.
* `third_party/htslib-build/htslib-src/vcf.c:2019-2026` — `bcf_read` dispatches to
  `vcf_read` (`vcf.c:3930-3936`, which returns `hts_getline`'s value when it is
  negative) or to `bcf_read1_core` (`vcf.c:1708`).
* `third_party/htslib-build/htslib-src/vcf.c:3749` — `vcf_parse` sets
  `int ret = -2` at its `err:` label, so a rejected record leaves `bcf_read`
  returning **-2**, never -1.  `vcf_parse1` is a macro for `vcf_parse`
  (`htslib/vcf.h:274`).
* `third_party/htslib-build/htslib-src/htslib/sam.h:974-976` — `bam_read1` doc:
  `-1 at end of file; < -1 on failure` (the contract `sam_read1` inherits).
* `third_party/htslib-build/htslib-src/htslib/hts.h:680-681` — `hts_getline` doc:
  `-1 on end-of-file; <= -2 on error`; `hts.c:2031` sets `-2` when the underlying
  read failed.
* `third_party/htslib-build/htslib-src/hts.c:2060` + `:2072` — htslib's own
  `hts_readlist` loops on `>= 0` and then tests `if (ret < -1) // Read error`.
  **The correct idiom is in the library we link**, and native's loops did not
  use it.

GATK's side: GATK rejects an unreadable record rather than continuing.  The
previous round established the mechanism (`htsjdk.tribble.TribbleException` out
of `AbstractVCFCodec.createGenotypeMap`, GATK exit 3) with file:line; the same
applies here, and the tools' own error convention is
`throw std::runtime_error("BAD_INPUT: ...")` → the tool's name + message on
stderr, exit 2.

## 4. STEP 2 — the gate

`fastgatk-native/scripts/verify_malformed_input_fail_loud_oracle.py`
(`main() -> int`, argparse, `tempfile.TemporaryDirectory`, final JSON status line,
non-zero on any mismatch).

10 tools × 5 cases = 50 cases.  Contract per case:

* `malformed-record-first` — an unreadable record at chr1:10 followed by two
  well-formed records: must exit non-zero, print `BAD_INPUT:`, and write **no**
  data record (nothing well-formed precedes the failure).
* `malformed-record-last` — the mirror image: must exit non-zero, print
  `BAD_INPUT:`, and the output may contain **only** the two-record prefix that a
  streaming writer legitimately flushed before it noticed (this is exactly what
  pinned GATK does before it exits 3), never the unreadable record and never
  anything beyond.
* `truncated-gzip` — 3000 well-formed records compressed and cut at 90 %: must
  exit non-zero, print `BAD_INPUT:`, and the output must be **provably
  incomplete** (fewer records than the 3000 the stream can hold).  The fixture is
  deliberately larger than one BGZF decompression buffer so the cut lands in the
  record stream rather than in the header read.
* `control-plain` / `control-gzip` — the same two well-formed records, plain and
  gzipped: must still exit 0 and emit exactly 2 records.  This is the
  over-rejection control the task asks for.

Exit status: **before the fix, 1** (`status: divergence`, 53 violations across 9
of the 10 tools — `.diag/audit-work/gate-before.log`).  **After the fix, 0**
(`status: pass`, `violations: []` — `.diag/audit-work/gate-after.log`).

Registration: the gate is **not** registered in CTest — hard constraint 1 forbids
editing `fastgatk-native/CMakeLists.txt`.  The orchestrator should add, next to
the existing malformed oracle:

```cmake
    add_test(NAME fastgatk-malformed-input-fail-loud-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_malformed_input_fail_loud_oracle.py")
    set_tests_properties(fastgatk-malformed-input-fail-loud-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
```

## 5. STEP 3 — the fix

One new header-only helper plus one line changed per call site: **+110 / −52
across 27 files**, no restructuring, no Kokkos ABI change, no vendored-source
change, no behaviour change for well-formed input.

`fastgatk-core/include/fastgatk/io/hts_read_guard.hpp` (new) — four inline
wrappers that perform the read and raise `BAD_INPUT` on `< -1`, so a loop
condition can only observe `0` (a record) or `-1` (clean end of input):

```cpp
inline int read_variant_record(htsFile* file, const bcf_hdr_t* header,
                               bcf1_t* record, const std::string& path) {
    const int status = bcf_read(file, header, record);
    if (status < -1)
        throw_read_failure("malformed or truncated VCF/BCF record", status, path);
    return status;
}
```

plus `read_alignment_record` (`sam_read1`), `read_indexed_alignment_record`
(`sam_itr_next`) and `read_text_line` (`hts_getline`), all raising

```
BAD_INPUT: <what> (htslib status <n>): <path>
```

htslib exposes no per-file error-string accessor, so the diagnostic carries the
status and the path; htslib has already logged its own `[E::vcf_parse_*]` /
`[E::bgzf_read]` line naming the offending record.  The helper lives in
`fastgatk-core/include` (not `fastgatk-native/include`) because
`fastgatk-core/src/hts_reader.cpp` needs it and only
`${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include` is on that target's
include path — and `fastgatk-hts` exports that directory `PUBLIC`, so every tool
sees it transitively.  No CMake change was needed or made.

Per-file diff (`added/deleted`):

| file | ± | sites |
| --- | --- | --- |
| `fastgatk-core/src/hts_reader.cpp` | 6/3 | `:470` sam_read1, `:486` sam_itr_next, `:559` hts_getline |
| `fastgatk-native/include/fastgatk/io/intervals.hpp` | 2/1 | `:97` hts_getline |
| `annotate_intervals_tool.cpp` | 2/1 | `:215` |
| `apply_vqsr_tool.cpp` | 5/3 | `:294` hts_getline, `:425`, `:780` bcf_read |
| `bqsr_tool.cpp` | 4/1 | `:3266` sam_read1 |
| `check_reference_compatibility_tool.cpp` | 2/1 | `:100` hts_getline |
| `collect_f1r2_counts_tool.cpp` | 3/1 | `:494` sam_read1 |
| `combine_gvcf_tool.cpp` | 6/2 | `:1158`, `:1709` bcf_read |
| `count_bases_in_reference_tool.cpp` | 2/1 | `:183` hts_getline |
| `fasta_reference_tool.cpp` | 3/2 | `:175` hts_getline, `:433` bcf_read |
| `filter_intervals_tool.cpp` | 2/1 | `:296` hts_getline |
| `filter_mutect_tool.cpp` | 11/5 | `:4003`, `:4669`, `:4726`, `:4819`, `:5093` bcf_read |
| `gather_vcfs_tool.cpp` | 2/1 | `:317` bcf_read |
| `genotype_gvcf_tool.cpp` | 2/1 | `:4670` hts_getline |
| `get_pileup_summaries_tool.cpp` | 8/4 | `:482`, `:842`, `:934` bcf_read, `:571` hts_getline |
| `left_align_tool.cpp` | 3/1 | `:1005` |
| `mark_duplicates_tool.cpp` | 5/2 | `:1113`, `:1350` sam_read1 |
| `mutect2_tool.cpp` | 2/1 | `:320` |
| `preprocess_intervals_tool.cpp` | 2/1 | `:182` hts_getline |
| `read_metrics_tool.cpp` | 2/1 | `:779` hts_getline |
| `reblock_gvcf_tool.cpp` | 3/1 | `:1847` |
| `select_variants_tool.cpp` | 4/2 | `:541`, `:2358` |
| `sort_sam_tool.cpp` | 5/2 | `:755`, `:838` sam_read1 |
| `split_intervals_tool.cpp` | 2/1 | `:324` hts_getline |
| `validate_variants_tool.cpp` | 4/2 | `:556`, `:777` |
| `variant_filtration_tool.cpp` | 5/3 | `:2152`, `:2209`, `:2521` |
| `variant_recalibrator_tool.cpp` | 7/4 | `:652`, `:2396`, `:2709` bcf_read, `:895` hts_getline |
| `variants_to_table_tool.cpp` | 6/3 | `:276` hts_getline, `:1117` bcf_read |

Two call sites needed one extra line of context rather than a bare substitution:

* `sort_sam_tool.cpp:838` — the merge loop is *outside* the `for (const auto& path
  : runs)` scope, so the path is `runs[node.run].string()` (`merge_runs` is pushed
  in `runs` order).
* `variants_to_table_tool.cpp:276` — `RawInfoReader::next()` is a member function
  and the constructor's `path` is not stored, so the class gained a
  `std::string source_path_` member initialised from the constructor argument.

## 6. Sites judged intentional / left alone, with the argument

* **Class B (12 sites)** — already correct; changing them would be a no-op diff.
  Quoted above.  They are the reason the class is smaller than the grep suggested:
  the idiom *already* had a correct form in the tree.
* **Class C (5 sites)** — `bqsr_tool.cpp:3217` and `sort_sam_tool.cpp:817` throw on
  *any* negative `sam_read1`, i.e. they already refuse to distinguish; the other
  three are comments.  Nothing to fix.
* **Class D (indexed iterators)** — `<0` would be wrong there: `-1` means "this
  interval's iterator is done".  Left as `< -1`, which is what the already-correct
  sites do.
* **`gather_vcfs_tool.cpp:317`** — fixed, but honestly labelled: the main copy loop
  is guarded and re-reads the same file, so this site was never *independently*
  able to cause silent loss.  It was fixed because treating a parse failure as
  "this shard is empty" is wrong for the reorder scan.
* **SAM text input** — deliberately *not* fixable here.  For SAM (not BAM),
  htslib's `sam_read1_sam` (`sam.c:4270`) returns `sam_parse1`'s `-1` for a
  malformed line, i.e. htslib itself makes a bad SAM line indistinguishable from
  EOF (`sam.c:4363-4369`: it only logs `Parse error at line N` and returns the
  negative value).  A `status < -1` guard cannot catch that without patching
  htslib, which is forbidden.  Native's BAM/CRAM path is covered; the SAM *text*
  path inherits htslib's own collapse.  This is the one part of the class that
  remains open, and it is htslib's contract, not native's loop.

## 7. STEP 4 — stale assertions

**None.  No test line was corrected.**

Two independent checks:

1. Structural: every script mentioning "malformed"/"truncated" was read for an
   assertion on a *malformed input stream*:
   * `verify_get_pileup_summaries.py:105` asserts `returncode == 0` for its
     "malformed-cigar" SAM — but that record is a zero-length CIGAR element that
     htslib *accepts* (`sam_read1` returns 0, the record is retained and
     `GoodCigarReadFilter` rejects it downstream).  It never reaches the `-2`
     path, so the assertion is correct and untouched.
   * `verify_count_bases_in_reference.py:116`, `verify_count_reads.py:222`,
     `verify_flag_stat.py:165` assert the *opposite* (non-zero for a malformed
     interval file) — consistent with this round's contract.
   * every other "truncated" hit refers to
     `testdata/real/cnv_somatic/human_g1k_v37.chr-20.truncated.fasta`, a
     well-formed 1 Mb reference, or to a truncated *tag*.
2. Empirical: the full 301-test dual-backend suite passes after the change (§8d),
   so no registered test depended on the truncation.

## 8. STEP 5 — gate results (mandatory order)

| step | command | result |
| --- | --- | --- |
| a | `python3 fastgatk-native/scripts/verify_malformed_input_fail_loud_oracle.py` **before the fix** | **exit 1** — 53 violations, 9/10 tools (`.diag/audit-work/gate-before.log`) |
| a | same gate, final tree | **exit 0** — `{"status": "pass", "cases": 50, "violations": []}` (`.diag/audit-work/gate-after.log`) |
| b | `ctest --test-dir fastgatk-native/build -j 8 -R 'genotype-gvcf\|mutect2\|vqsr\|variant-recalibrator' -V` | **exit 0** — **64/64 passed** (`.diag/audit-work/ctest-subset2.log`) |
| c | the 18-name strict-gate filter of the task, `-j 8 -V` | **exit 0** — **20/20 passed** (`.diag/audit-work/ctest-strict.log`) |
| d | `fastgatk-native/scripts/run_regression.sh --label bcfread-audit` | **exit 0** — **301/301 omp** (1390.14 s) and **301/301 serial** (1500.72 s) (`.diag/regression/20260912-154012/{omp,serial}.log`, `summary.txt`) |

`-j 8` is the house convention for the same suites
(`run_regression.sh:95` builds `ctest … -j ${JOBS} --output-on-failure` with
`JOBS=8` by default); the filters and the pass criteria are unchanged, and the
serial run of (b) was abandoned after 4 tests only because eight-way parallelism
is what the repository's own runner uses.

**Nothing was edited after the suite run.**  The 19 changed tools were rebuilt in
both trees immediately before it (`cmake --build`, exit 0 each) and their md5s
are byte-identical before and after the suite
`.diag/audit-work/md5-before-suite.txt` vs `-after-suite.txt` (diff empty), so no
re-run was needed.  No file under `fastgatk-native/src`, `fastgatk-native/include`,
`fastgatk-core/src`, `fastgatk-core/include` or `fastgatk-native/scripts` is newer
than the newest binary in either tree.

## 9. Tree state

* 28 modified files + 2 new files (the helper header and the gate script);
  `git diff --numstat` totals **110 added / 52 deleted**.  No commit, no branch,
  changes left in the working tree.
* Both trees rebuilt (`cmake --build fastgatk-native/build` and
  `fastgatk-native/build-serial`), exit 0.  Staleness check: no file under
  `fastgatk-native/src`, `fastgatk-native/include`, `fastgatk-core/src`,
  `fastgatk-core/include` is newer than the newest binary in either build tree.
* NOT touched: `fastgatk-native/CMakeLists.txt`, every root `*.md`, every test
  script, `third_party/`, the Kokkos kernel ABI.
* Artifacts: probes `.diag/audit-work/probe.py` (+ `probe-before.json`,
  `probe-after.json`), `.diag/audit-work/truncgz2.py`, `.diag/audit-work/bamtrunc.sh`,
  `.diag/audit-work/bamflip.py`, `.diag/audit-work/bammulti.py`; the census
  working notes `.diag/audit-work/sites.txt`, `tail.txt`, `post.txt`; the applied
  edit scripts `.diag/audit-work/apply_fix.py`, `apply_fix_getline.py`; the full
  diff `.diag/audit-work/full.diff`; gate logs `gate-before.log` / `gate-after.log`;
  test logs `ctest-subset2.log`, `ctest-strict.log`, `regression.log`;
  binary hashes `md5-before-suite.txt`, `md5-after-suite.txt`.
* Gate script md5 `406ed1bb7887744232409594a31e2538`,
  helper header md5 `1da9b3a0f34fae8c2e758423e90333d8`.

## 10. What remains unproven

* The gate is **not registered** in CTest (constraint 1); it was run by hand.
* 14 of the 51 fixed sites are directly exercised by the measurements above
  (`variants_to_table:276`, `:1117`, `select_variants:541`, `:2358`,
  `left_align:1005`, `variant_filtration:2521`, `validate_variants:777`,
  `reblock_gvcf:1847`, `combine_gvcf:1709`, `filter_mutect:5093`,
  `fasta_reference:433`, `sort_sam:755`, `mark_duplicates:1113`,
  `collect_f1r2_counts:494`).  The remaining **37 sites are fixed by inference,
  not by measurement**, and are labelled as such in §1.
* The BAM half was measured with a *truncated* BAM (a `-2` from `bgzf_read`); a
  BAM corrupted in a way that produces a *parse* -2 inside a valid BGZF block was
  not constructed.  The SAM **text** path cannot be fixed at all (§6).
* `hts_getline` sites were fixed by inspection of the same contract; none was
  exercised with a truncated `.gz` interval list or GATKReport.
* Multi-input (`-V a -V b`) behaviour: pre-fix the rest of the affected input was
  lost while the other inputs continued; post-fix the run aborts.  Only
  single-input cases were measured.
* The `gendb://` GenomicsDB bridge path reads through a separate process and was
  not exercised; `genomicsdb_import_tool.cpp` was already class B.
* Residual contract boundary, stated plainly: for inputs the tool cannot read,
  native leaves whatever it had already flushed (exit non-zero), exactly as pinned
  GATK does.  "Fails loudly, output provably incomplete" is the asserted
  contract; "no output file at all" is **not** asserted for streaming writers,
  because achieving it would require restructuring every tool's writer.
