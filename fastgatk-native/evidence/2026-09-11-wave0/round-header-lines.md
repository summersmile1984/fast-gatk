# Round: the output header's CONTENT lines (GATK's own declarations, htsjdk re-quoting, MLEAC/MLEAF wording)

Scope: the divergence class the previous two rounds *measured but did not fix*
(`.diag/round-header-filter.md` §1.4, `.diag/round-filter-pass-line.md` §1.3):
the existing gates compared DATA ROWS only, so an entire class of **header
content** divergences was invisible.  On the fixture of
`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`, pinned GATK
4.6.2.0 GenotypeGVCFs writes `##INFO`/`##FORMAT` declarations native omitted,
re-quotes every `Description` native kept unquoted, and words `MLEAC`/`MLEAF`
differently.  Gate first, fix second.

## 0. Verdict

* **GATK truth (measured, pinned 4.6.2.0):** for the gate fixture, GATK's output
  header has **28 content lines** (`##` lines minus the two producer-identity
  ones); native had 23 of them.  Exactly:

  ```
  (a) missing in native (15)      (b) extra in native (10)     (c) 10 keys with different text
  ```

  §1.1 lists every literal line.  They fall into three families:
  1. **declarations GATK adds and native never wrote** (7 lines): the three
     rank-sum INFO lines, the standard `##INFO=<ID=DP,...>` line, the standard
     `##FORMAT=<ID=AD,...>` line, and the long-wording `MLEAC`/`MLEAF` pair;
  2. **htsjdk re-quoting** (7 lines): the input's own `##ALT`, `##INFO=<ID=DP`,
     `##INFO=<ID=AD`, `##FORMAT=<ID=AD|DP|GQ|GT|PL` lines come back with
     `Description="..."` because htsjdk writes a header out of parsed objects;
  3. **native's own shorter spellings** of `MLEAC`/`MLEAF` that GATK never
     writes.
* **GATK rule with source:** `GenotypeGVCFsEngine.setupVCFWriter()` seeds its
  header set from the input (`new LinkedHashSet<>(
  inputVCFHeader.getMetaDataInInputOrder())`, `GenotypeGVCFsEngine.java:395`),
  strips only `GVCFBlock` (`:398-399`), then adds
  `annotationEngine.getVCFAnnotationDescriptions(false)` (`:401`),
  `genotypingEngine.getAppropriateVCFInfoHeaders()` (`:402`), MLEAC/MLEAF/RGQ
  (`:404-406`), `VCFStandardHeaderLines.getInfoLine(DEPTH_KEY)` (`:407`) and the
  LowQual filter (`:416`).  The annotation set is the unconditional
  `StandardAnnotation` group (`GenotypeGVCFs.java:255-257`).  The set
  de-duplicates by the **whole line**, not by ID — measured, GATK's header
  carries two `##INFO=<ID=DP,...>` and two `##FORMAT=<ID=AD,...>` lines with
  different descriptions.  `getDescriptions()` resolves through
  `GATKVCFHeaderLines.getInfoLine/getFormatLine(key, true)`
  (`VariantAnnotation.java:22-33`, `GATKVCFHeaderLines.java:18-43`), which falls
  through to htsjdk's `VCFStandardHeaderLines` for `AD`/`DP`.
* **Fix kept in the tree:** one local change in
  `fastgatk-native/src/genotype_gvcf_tool.cpp` (`+165/-14`): a
  `gatk_canonical_header_line()` re-quoter, two GATK declaration lists injected
  into the existing INFO/FORMAT groups by exact text, and three compat-gated
  append texts (DP/MLEAC/MLEAF, plus NDA) so native no longer writes its own
  wording into the GATK-compatibility header.  No restructuring, no kernel ABI,
  no `CMakeLists.txt`, no Mutect2, no other tool.
* **Gate:** `verify_genotype_gvcf_spandel_gatk_oracle.py` extended with a
  whole-content assertion for every gated case.  **exit 1 before** the fix
  (29 violations, one per gated case, all of them that one assertion — §2.1),
  **exit 0 after** (31 cases, 0 violations).
* Step 5: (a) exit 0; (b1) exit 0; (b2) 15/15 exit 0; (c) 19/19 exit 0;
  (d) see §6.  Change **kept in the tree**.

## 1. STEP 1 — the full grouped header diff, byte for byte

Reproducer: `.diag/header_lines_probe.py` — log/JSON
`.diag/header-lines-probe-before.{log,json}` (before) and
`.diag/header-lines-probe-after.{log,json}` (after, final).  It runs both tools
on the oracle's own fixture (its 100 bp `chr1` reference, `IndexFeatureFile`,
`--gatk-compatible-annotations` on the native side) and classifies every
difference between the two **full header line lists** into (a) missing in
native, (b) extra in native, (c) same key + different text, (d) order-only.
`##GATKCommandLine` is excluded from the diff as producer identity (§4);
`##source=GenotypeGVCFs` is byte-identical on both sides so it never appears.

### 1.1 The fixture that the gate pins (`gate-fixture-star`)

GATK 31 lines (30 + `#CHROM`), native 25 before / 30 after:

```
                      BEFORE                                    AFTER
(a) missing in native  15                                        0
(b) extra in native    10                                        0
(c) same key, diffs 10 keys                                      0
(d) order-only          0                                        14   (unchanged divergence, §3.3)
content lines         GATK 28  NATIVE 23                        GATK 28  NATIVE 28
```

**(a) missing in native (15)** — every literal line, GATK's exact bytes:

```
##ALT=<ID=NON_REF,Description="Represents any possible alternate allele">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
##INFO=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
```

**(b) extra in native (10)** — native's own bytes for the same declarations:

```
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=Maximum likelihood allele count>
##INFO=<ID=MLEAF,Number=A,Type=Float,Description=Maximum likelihood allele frequency>
##INFO=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
```

**(c) same key, different text (10 keys):** `ALT/NON_REF`, `FORMAT/AD`,
`FORMAT/DP`, `FORMAT/GQ`, `FORMAT/GT`, `FORMAT/PL`, `INFO/AD`, `INFO/DP`,
`INFO/MLEAC`, `INFO/MLEAF` — i.e. exactly the (a)/(b) pairs above: 7 of them are
"unquoted vs quoted", 1 (`INFO/DP`) is "quoted input line vs quoted standard
line" (both present in GATK, native had only the input's), and 2
(`INFO/MLEAC|MLEAF`) are the wording divergence, where native's line is a
*substitution* rather than a missing declaration.

**(d) order-only: 0 before** — the two multisets differed, so the order
comparison is meaningless before the fix.  **After the fix: 14 positions**, all
of them the pre-existing whole-header ordering divergence (§3.3), not touched.

### 1.2 Every other fixture measured (after the fix)

`(a)`/`(b)`/`(c)` are **0 for all of them**; only `(d)` order-only remains:

| fixture | what it isolates | after: (a)/(b)/(c) | after: (d) |
| --- | --- | --- | --- |
| `gate-fixture-star` | the gate's own fixture, default mode | 0/0/0 | 14 |
| `gate-fixture-star-dense` | `--include-non-variant-sites` | 0/0/0 | 14 |
| `ordinary-variant-dense` | the ordinary variant path | 0/0/0 | 14 |
| `quoted-input-header` | input Descriptions already quoted: isolates "missing declarations" from "re-quoting" | 0/0/0 | 14 |
| `input-declares-gatk-mleac` | input declares GATK's long MLEAC/MLEAF (the HaplotypeCaller shape): pins that native must NOT duplicate them | 0/0/0 | 14 |
| `description-with-comma` | `Description="Read depth, filtered"`: quoting inside a value | 0/0/0 | 14 |
| `input-without-info-dp` | no `INFO/DP` and no `MLEAC` in the input: exercises native's own append path | 0/0/0 | 12 |
| `num-discovered-alleles-flag` | `--annotate-with-num-discovered-alleles` (NDA) | 0/0/0 | 18 |
| `reordered-attribute-keys` | input with `ID`,`Description`,`Type`,`Number` order | *GATK rejects the input* (exit 3, *"Tag Description in wrong order (was #2, expected #4)"*) — no parity obligation | — |

Before the fix, the same probe reported, for `quoted-input-header`,
(a)=7 (only the declarations, no re-quoting) and for `reordered-attribute-keys`
the same rejection — the two controls that separate the halves:
`.diag/header-lines-probe-before.log`.

### 1.3 Does the same gap exist for other tools?

Not measured this round, and not changed.  The two functions touched are in
`genotype_gvcf_tool.cpp`: `gatk_compatible_header_text()` (the
`--gatk-compatible-annotations` boundary, used by **both** GenotypeGVCFs
writers: streaming `:5817`, aggregate `:6648`) and the two header-append helpers
of the same tool.  `grep -rn "BaseQRankSum" fastgatk-native/src/` shows the
other tools declare the rank-sum lines themselves.  No other tool was touched.

## 2. STEP 2 — the gate, extended and run BEFORE the fix

**Extended** the already-registered
`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` rather
than adding a new script, for the same reason as the previous two rounds: it
already runs pinned GATK and native on the same fixtures with the same
arguments, and it is already registered in `fastgatk-native/CMakeLists.txt` under
the strict-gate filter of step 5(c).  `CMakeLists.txt` must not be edited this
round, so a brand-new file could not have been registered as a ctest gate.

Added (oracle diff `+224/-2`):

* a docstring section *"The output header's other content lines"* quoting the
  GATK `file:line`s of every added line (§3.1) and the measured rule that GATK's
  set de-duplicates by the whole line;
* `HEADER_QUOTED_DESCRIPTIONS` and `HEADER_WITH_GATK_MLE_LINES` fixtures;
* `content_header_lines()` — every `##` line except the two producer-identity
  ones (`##GATKCommandLine`, `##source`) — and a **gated** assertion that the two
  outputs' content lines are equal **as a multiset** (order aside, duplicates
  included);
* `GATK_DECLARED_HEADER_LINES` (the 19 lines measured to be GATK's own) and a
  separate violation if GATK stops declaring one of them, so the comparison
  cannot pass for the wrong reason;
* three new **GATED** cases: `header-content-with-unquoted-descriptions`,
  `header-content-with-quoted-descriptions`, `header-content-with-gatk-mle-wording`;
* the JSON payload and the printed report now carry
  `gatk_content_header_lines`, `native_content_header_lines`,
  `header_content_missing_in_native` and `header_content_extra_in_native`.

The previous rounds' `##FILTER` group assertion is unchanged and still gated
**in order**; order-only differences elsewhere stay observations.

### 2.1 Exit status and the literal lines, BEFORE the fix

**exit 1**, `"status": "divergence"`, 31 cases, **29 violations** — one per gated
case, every one of them the new content assertion and nothing else
(`.diag/header-lines-gate-before.log`), e.g.:

```
[default] VIOLATION: header content lines are not byte-identical (as a set):
  missing in native=['##ALT=<ID=NON_REF,Description="Represents any possible alternate allele">',
   '##FORMAT=<ID=AD,...Description="Allele depths">',
   '##FORMAT=<ID=AD,...Description="Allelic depths for the ref and alt alleles in the order listed">',
   ... 15 lines ...]
  extra in native=['##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>',
   ... '##INFO=<ID=MLEAC,...Description=Maximum likelihood allele count>', ... 10 lines ...]

[header-content-with-quoted-descriptions] VIOLATION: ... missing in native=
  ['##FORMAT=<ID=AD,...Description="Allelic depths for the ref and alt alleles in the order listed">',
   '##INFO=<ID=BaseQRankSum,...', '##INFO=<ID=DP,...Approximate read depth; some reads may have been filtered',
   '##INFO=<ID=MLEAC,...long wording...', ...]  extra in native=['##INFO=<ID=MLEAC,...Maximum likelihood allele count>', ...]
```

The three new cases behave exactly as designed before the fix: the
unquoted-description case reports both halves, the **quoted**-description case
reports *only* the missing declarations (its diff contains no re-quoting pair),
and the MLEAC case reports the wording substitution.

**AFTER the fix: exit 0**, `"status": "pass"`, 31 cases, 0 violations
(`.diag/header-lines-gate-after.log`).

## 3. STEP 3 — the fix

`git diff --stat` = **2 files changed, 389 insertions(+), 16 deletions(-)**
(oracle `+224/-2`, tool `+165/-14`; the line counts include the comment blocks
that carry the `file:line` provenance).

### 3.1 What was changed, and why it is the correct minimal shape

All of it is header construction in `genotype_gvcf_tool.cpp`:

1. **`gatk_canonical_header_line()`** (new, ~20 lines of code): for a
   `##INFO=<`/`##FORMAT=<`/`##FILTER=<`/`##ALT=<` line whose `Description=`
   value is not quoted, wrap the value in quotes.  This is htsjdk's
   re-serialization (`VCFHeaderLineTranslator.parseLine()` builds a
   `VCFCompoundHeaderLine`; the writer emits the line from that object).  It is
   applied to every retained line inside `gatk_compatible_header_text()` —
   including the input's PASS line re-read from the raw input text by the
   previous round's `input_pass_filter_line()`, which is canonicalized at the
   point of insertion — and it is a no-op for lines that are already quoted.
   *The attribute order inside `<>` is not touched*, because htsjdk refuses an
   input that does not already use its order (measured, §1.2), so no valid input
   needs re-ordering.
2. **`gatk_declared_info_lines()` / `gatk_declared_format_lines()`** (new, data):
   GATK's own declarations, measured byte for byte and cited per line.  They are
   appended to the existing `info_lines`/`format_lines` groups **only when their
   exact text is absent**, before the existing sorts — i.e. the same
   collect/sort/re-insert shape the function already had, with GATK's
   whole-line de-duplication rule instead of native's "is the ID present"
   rule.  This is what reproduces GATK's duplicate `INFO/DP` and `FORMAT/AD`
   pairs (the input's line is a different line and is kept alongside) while
   keeping the canonical HaplotypeCaller `MLEAC`/`MLEAF` pair single.
3. **Three compat-gated append texts.**  Native produced `MLEAC`, `MLEAF`, `DP`
   and `NDA` declarations with descriptions of its own from the two header
   builders that feed *both* output profiles.  In the GATK-compatibility profile
   they now use GATK's text (`kGatkMleacInfoLine`/`kGatkMleafInfoLine`/
   `kGatkStandardDpInfoLine`/`kGatkNdaInfoLine`); the diagnostic profile keeps
   the exact bytes it wrote before.  The gate is `options.
   gatk_annotation_compatibility`, the same pattern the file already uses for
   the `GP`/`PG` appends, so the non-parity diagnostic header is byte-unchanged.

Provenance of every added declaration, `file:line`:

| line | GATK source |
| --- | --- |
| `##INFO=<ID=BaseQRankSum,...>` | `BaseQualityRankSumTest` (StandardAnnotation) + `GATKVCFHeaderLines.java:151` |
| `##INFO=<ID=MQRankSum,...>` | `MappingQualityRankSumTest` + `GATKVCFHeaderLines.java:168` |
| `##INFO=<ID=ReadPosRankSum,...>` | `ReadPosRankSumTest` + `GATKVCFHeaderLines.java:195` |
| `##INFO=<ID=DP,..."Approximate read depth; some reads may have been filtered">` | `VCFStandardHeaderLines.getInfoLine(DEPTH_KEY)` via `GenotypeGVCFsEngine.java:407` |
| `##INFO=<ID=MLEAC,...>` / `MLEAF` | `GATKVCFHeaderLines.java:148-149` via `GenotypeGVCFsEngine.java:404-405` |
| `##FORMAT=<ID=AD,..."Allelic depths for the ref and alt alleles in the order listed">` | `DepthPerAlleleBySample` via `VariantAnnotation.java:22-33` -> `GATKVCFHeaderLines.java:32-43` -> htsjdk `VCFStandardHeaderLines` |
| `##FORMAT=<ID=RGQ,...>` | `GATKVCFHeaderLines.java:133` via `GenotypeGVCFsEngine.java:406` |
| `##INFO=<ID=NDA,...>` | `GATKVCFHeaderLines.java:205` via `GenotypingEngine.getAppropriateVCFInfoHeaders()` (`GenotypingEngine.java:90-94`) |
| `##INFO=<ID=AC/AF/AN/ExcessHet/FS/InbreedingCoeff/MQ/QD/SOR,...>` | the same `StandardAnnotation` group; native already wrote these bytes, so the injection is a no-op for them and only guarantees the rule |

Not restructured: the genotyping engine, the FILTER column handling, the AF/PL
kernels and their ABI, the `##INFO`/`##FORMAT` **ordering** logic, the
GVCFBlock/PASS handling, the `##source`/lowqual group logic, the
non-compatibility diagnostic profile, and every other tool.
`fastgatk-native/CMakeLists.txt` and every root `*.md` are unmodified.

### 3.2 Collateral, explicitly

* The two changed functions are called only from GenotypeGVCFs itself
  (`gatk_compatible_header_text` from `:5817` and `:6648`;
  `add_genotype_output_header_fields` from one streaming call site, plus the
  aggregate block).
* The three append-text gates add three `const std::string` locals in
  `add_genotype_output_header_fields` and change two `bcf_hdr_append` arguments
  in the aggregate path; in non-compatibility mode the appended bytes are
  character-identical to before, so the diagnostic header is unchanged.
* `output.reserve(...)` and the `##FILTER` logic are untouched.

### 3.3 What is deliberately NOT fixed (order, identity)

* **Order.** `(d) order-only` is 14 positions on the gate fixture and 12-18 on
  the others, unchanged by this round: native preserves the input's line order
  and its own INFO rank list, GATK emits htsjdk's fully sorted header.  The
  visible consequences are the `##contig` position (native 1, GATK 28), the
  `INFO/AD` position (native last, GATK right after `AC`) and the order of the
  two `INFO/DP` lines.  Out of scope by the task statement, reported only.
* The `##FILTER` group's absolute index (native 3, GATK 2) — same cause, same
  status (unchanged from the previous rounds).

## 4. Lines judged PRODUCER IDENTITY and deliberately left

| line | why it is identity, not content |
| --- | --- |
| `##GATKCommandLine=<ID=GenotypeGVCFs,CommandLine="...",Version="4.6.2.0",Date="...">` | It records the **command line, version and wall-clock date of the GATK process that wrote this file**.  Native is a different program with a different CLI, and it has no `--add-output-vcf-command-line` equivalent on this surface; reproducing the line would mean writing a false provenance record.  Excluded from the probe diff and from the gate. |
| `##source=GenotypeGVCFs` | The program that wrote the file.  Measured: **both sides write the identical byte string** (`##source=GenotypeGVCFs`), so there is no divergence to fix; it is excluded from the gated set because naming the writer is identity, and native is a different program. |

Everything else in the header is a declaration about the file's *contents*
(sequence dictionary, ALT/FILTER/FORMAT/INFO schemas, `##fileformat`) and is
therefore content: it is gated byte for byte (order aside).

## 5. STEP 4 — stale assertions

**No registered test pinned the old behaviour; no test line needed correcting.**
The lines checked, and why each was or was not stale:

| File / check | Why it could have been stale | Outcome |
| --- | --- | --- |
| `verify_gatk_genotype_gvcf.py:234-236` (`semantic_header`, a **sorted set** of every `##` line minus provenance) | compares GATK's and native's header sets | **not stale**: it compares the two tools, so this round can only make it *more* equal; its input is a GATK HaplotypeCaller gVCF, whose header already carries every declaration native now adds (so the additions are no-ops) and whose Descriptions are already htsjdk-canonical (so re-quoting is a no-op).  Ran green in step 5(b1) |
| `verify_gatk_genotype_gvcf_legacy_qual.py:105`, `verify_gatk_genotype_gvcf_multisample.py:116` (`assert native_header == gatk_header`) | same, order-insensitive `sorted({...})` | **not stale**, same reasoning; green in step 5(b2) |
| `verify_hc_dense_gvcf_genotype_gatk_oracle.py:128` | the strictest one: every non-ignored header line **in order** | **not stale and the sharpest evidence for "no over-declaration"**: its input gVCF is produced by HaplotypeCaller at run time, which already declares the rank-sum/MLEAC/MLEAF/AD/DP lines with GATK's own bytes, so "add only when the exact text is absent" must add nothing and the order must not change.  Green in 5(b2)/5(c) |
| `verify_genotype_gvcf_assignment_gatk_oracle.py:27-40,116-138` | parses header lines into a `(kind, id)` set, asserts `required <= ids` | **not stale**: subset assertions, which adding lines cannot break |
| `verify_gatk_genotype_gvcf_inbreeding.py:137-149`, `verify_gatk_genotype_gvcf_multiallelic.py:123-134` | header **ID** sets | **not stale**: subset assertions |
| `verify_genotype_gvcf_max_alternate_alleles_gatk_oracle.py:94` | asserts `NDA=7` in a data row | **not stale**: a row assertion; the NDA *header* text it never reads changed (step 3.1 item 3), and the oracle compares no headers.  Green in 5(b2) |
| `grep -rn 'Maximum likelihood allele count\|Number of discovered alternate alleles\|Approximate read depth' fastgatk-native/scripts/` | could pin native's own wording | only the source file itself (before this round) and this round's oracle text; no test script pins them |
| `fastgatk-native/scripts/verify_genotype_gvcf_starts_in_intervals_gatk_oracle.py:46-49`, `verify_reblock_*.py` | contain unquoted `Description=Approximate read depth` lines | **not stale**: those are **fixture inputs** of GenotypeGVCFs/reblock, not output-header expectations |

The only test-script change in this round is the oracle itself (§2).

## 6. STEP 5 — gate results (mandatory order)

| Step | Command | Result |
| --- | --- | --- |
| a | `python3 fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` | **exit 0** (31 cases, 0 violations); **exit 1 before** the fix (29 violations, §2.1) |
| b1 | `python3 fastgatk-native/scripts/verify_genotype_gvcf.py` | `B1_EXIT=0` — `{"status": "pass", "output_records": 1}` (`.diag/header-lines-verify-gvcf.log`) |
| b2 | `ctest --test-dir fastgatk-native/build -R 'genotype-gvcf' -V` | `B2_EXIT=0` — **15/15 passed** (`.diag/header-lines-ctest-genotype-gvcf.log`) |
| c | the 17-name strict-gate filter from the task | `C_EXIT=0` — **19/19 passed** (`.diag/header-lines-strict-gates.log`) |
| d | `fastgatk-native/scripts/run_regression.sh --label header-lines` | `REG_EXIT=0` — **300/300 on BOTH backends** (§6.1).  Both trees were rebuilt after the last source edit and before every gate above (`.diag/header-lines-build-{omp,serial}.log`) |

### 6.1 Dual-backend regression (d)

```
REG_EXIT=0
| 后端   | 构建目录                              | 结果 | 通过/总数 | 耗时        |
| omp    | OpenMP (fastgatk-native/build)        | 通过 | 300/300   | 1319.02 sec |
| serial | Serial (fastgatk-native/build-serial) | 通过 | 300/300   | 1302.89 sec |
```

Evidence block `.diag/regression/20260912-102149/{omp.log,serial.log}` (+
`.diag/header-lines-regression.out` for `REG_EXIT`), run window
10:21:49 → 10:43:48, label `header-lines`, git `7c021bf` + the two uncommitted
files of §8.  Both backends ran the full 300-test suite with
`FASTGATK_REQUIRE_GATK_ORACLE=1` (the runner's default) and ctest parallelism 8;
the runner's staleness check reported no warning (`grep -c 警告` = 0 in both
logs), i.e. no source file is newer than the binaries it tested.

**Nothing was edited after the builds.**  Provenance (timestamps):

```
09:35:26  fastgatk-native/src/genotype_gvcf_tool.cpp
09:35:56  fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py
09:36:09  fastgatk-native/build/fastgatk-genotype-gvcf          (omp, rebuilt)
09:38:19  fastgatk-native/build-serial/fastgatk-genotype-gvcf   (serial, rebuilt)
10:21:49 → 10:43:48   regression window (both backends)
```

`md5sum` of the two changed test inputs and the two binaries:

```
3f454bad4ed76d43df70258cf9c06ee1  fastgatk-native/src/genotype_gvcf_tool.cpp
9da55dd2c359fce7d0898265c40a93e3  fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py
fc4f8024ebc81f26eeddd41b748603c9  fastgatk-native/build/fastgatk-genotype-gvcf        (omp)
51bf95f0870908bcd9c41f65899f1748  fastgatk-native/build-serial/fastgatk-genotype-gvcf (serial)
```

`.diag/header-lines-binaries.md5` re-confirms the binary hashes.  The only files
written after the builds are this report, that hash file and the read-only gate
logs of step 5 — no source, test script, header or build input changed after the
binaries were built, so no rebuild and no re-run were required.

## 7. What remains unproven

* The gate pins the header content for the fixtures it runs (31 cases, 28 gated).
  §1.2's probe covers 9 fixture/mode combinations, and the source rule is read,
  but it is a point measurement, not a proof for every input shape.  Untested:
  `.bcf`/`.vcf.gz` inputs, `--dbsnp`, `--keep-combined-raw-annotations`, a
  GenomicsDB source, `--annotations-to-exclude`/`-A` (native has no handling for
  those flags at all, so a run that excludes an annotation would make native
  over-declare), and inputs whose lines are already quoted but with a
  *different* `Number`/`Type` for the same ID.
* The de-duplication rule was inferred as "htsjdk's `VCFHeaderLine` set compares
  the whole line (ID + Number + Type + Description)" from the measured duplicate
  `INFO/DP` and `FORMAT/AD` pairs and from the canonical-MLEAC control (§1.2).
  I did not read htsjdk's `equals()` out of the shaded jar this round.
* `gatk_canonical_header_line()` quotes the value from `Description=` up to the
  line's final `>`.  htsjdk's escaping of an embedded `"` inside a description
  is not measured, so a description containing a quote character is unverified
  (it is returned unchanged, i.e. the previous behaviour).
* The NDA header wording is fixed and verified by the probe
  (`num-discovered-alleles-flag`: 0/0/0), but it is **not covered by the
  oracle**: adding that fixture as a *gated* case failed for an unrelated reason
  — with `--annotate-with-num-discovered-alleles` GATK writes `NDA=1` in the data
  row while native writes no `NDA` at all (measured during this round:
  `GATK chr1 2 . AA A 92.60 . AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;NDA=1;QD=4.63 ...`
  vs native without `NDA`), so the case cannot be gated on rows without fixing a
  **data-row** divergence that is outside this round's scope.  That data-row
  divergence is measured here but not repaired.
* `(d) order-only` remains 12-18 positions per fixture, and the `##FILTER`
  group's absolute index remains off by one; both are the pre-existing
  whole-header ordering divergence, measured but not repaired, and this round
  deliberately does not gate them.
* "No registered test pinned the old behaviour" (§5) is a grep plus the step-5
  results, not a proof: a test could construct such an expectation dynamically.
  Step 5(d) is the empirical check (300/300 on both backends with
  `FASTGATK_REQUIRE_GATK_ORACLE=1`), which is evidence, not proof.

## 8. Tree state

* Change **kept in the tree** (no commit, no branch, no branch switch).
* Modified — exactly two files: `fastgatk-native/src/genotype_gvcf_tool.cpp`
  and `fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`.
* NOT touched: `fastgatk-native/CMakeLists.txt`, every root `*.md`, Mutect2,
  every other tool, and every registered test script other than the oracle.
* Artifacts: probe `.diag/header_lines_probe.py` (+
  `.diag/header-lines-probe-{before,after}.{log,json}`), gate logs
  `.diag/header-lines-gate-{before,after}.log`, build logs
  `.diag/header-lines-build-{omp,serial}.log`, step-5 logs
  `.diag/header-lines-{verify-gvcf,ctest-genotype-gvcf,strict-gates,regression}.log|out`,
  regression evidence `.diag/regression/20260912-102149/`, hashes
  `.diag/header-lines-binaries.md5`.
