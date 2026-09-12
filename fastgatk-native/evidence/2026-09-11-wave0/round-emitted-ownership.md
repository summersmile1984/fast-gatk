# Round: the recorded-deletion set must be built from EMITTED alleles

Scope: the REPORTED-ONLY case
`unemitted-upstream-deletion-star-plus-concrete-alt` in
`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`, i.e. the
difference between GATK's `recordDeletions()` (which sees only what a locus
actually emitted) and native's input-derived span list, which decides whether a
symbolic spanning deletion `*` is orphaned.

## 0. Verdict

* The hand-off hypothesis is **CONFIRMED and now fixed**. It is *not* a
  traversal/buffering problem and *not* a stage-boundary problem: every consumer
  of the ownership answer already lived inside `apply_gatk_output_allele_subset()`
  on the pipeline's single-threaded, in-order compute worker, while the value was
  produced earlier, in the union/merge step, from a pre-pass over the **input**
  records. Moving the answer to where GATK computes it (`GenotypingEngine.java:314`,
  inside `calculateOutputAlleleSubset()`) and recording what the locus emitted at
  GATK's `:178-179` point is **state bookkeeping** — no reorder, no buffer, no
  new stage.
* Three measured differences, in **both** directions, all now gated:
  1. native counted deletions from records GATK dropped (the reported case);
  2. native's interval ended at `start + refLen` instead of
     `start + (refLen - altLen)`, so an emitted ALT of two or more bases
     over-covered — a **second, previously unmeasured** divergence;
  3. `recordDeletions()` records an interval for an emitted `*` itself (its
     length is 1), which native's concrete-deletion scan skipped — the
     **opposite direction**, native pruning a `*` GATK keeps.
* Production change: one translation unit
  (`fastgatk-native/src/genotype_gvcf_tool.cpp`, +130/-106 including comments).
  No kernel, no ABI, no `CMakeLists.txt`, no other tool, no root `.md`.
* Oracle: **exit 1 before** the fix (3 violations, §3), **exit 0 after**
  (41 gated cases, 0 violations).
* Step 5: (a) exit 0, (b) exit 0 + 18/18, (c) exit 0 + 23/23,
  (d) 304/304 on **both** backends (`REG_EXIT=0`, §5.1).

## 1. STEP 1 — the GATK order, with `file:line`

`gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/genotyper/GenotypingEngine.java`:

| line | what it does |
| --- | --- |
| `:155` | `calculateOutputAlleleSubset(AFresult, vc, forcedAlleles)` — the output-ALT subset is **decided here**, over the AF result of the merged locus |
| `:167-169` | `return null` for a locus that does not pass the emit threshold → the traversal writes nothing *yet* |
| `:173-175` | `return null` when the only surviving ALT is `*` and the traversal is not `EMIT_ALL_ACTIVE_SITES` |
| `:178` | `outputAlleles = outputAlternativeAlleles.outputAlleles(vc.getReference())` — the concrete emitted allele list (reference first) |
| `:179` | `recordDeletions(vc, outputAlleles)` — **only now** are this locus' deletions recorded, and only for the alleles that survived `:155` |
| `:296-328` | `calculateOutputAlleleSubset()` itself |
| `:312` | `isPlausible = afCalculationResult.passesThreshold(allele, standardConfidenceForCalling)` |
| `:314` | `isSpuriousSpanningDeletion = isSpanningDeletion(allele) && !isVcCoveredByDeletion(vc)` — the ownership test, consulted **before** this locus records anything |
| `:316` | `toOutput = (plausible \|\| …) && !isSpuriousSpanningDeletion` — a spurious `*` is dropped from the output ALT list |
| `:318` | `siteIsMonomorphic &= !(isPlausible && !isSpuriousSpanningDeletion)` |
| `:330-332` | `clearUpstreamDeletionsLoc()` — `@VisibleForTesting`; **no production caller** (`grep -rn clearUpstreamDeletionsLoc gatk-source/src` → only `GenotypingEngineUnitTest.java:68`) |
| `:343-357` | `recordDeletions(vc, emittedAlleles)`: lazily cull the head of `upstreamDeletionsLoc` (`:344-346`), then for each emitted allele with `deletionSize = vc.getReference().length() - allele.length() > 0` add `SimpleInterval(contig, vc.getStart(), vc.getStart() + deletionSize)` (`:349-355`) |
| `:365-371` | `isVcCoveredByDeletion(vc)`: `!upstreamDeletionsLoc.isEmpty() && …anyMatch(loc -> loc.getContig().equals(vc.getContig()) && loc.getStart() < vc.getStart() && vc.getStart() <= loc.getEnd())` |
| `:52` | the state: `private final PriorityQueue<Locatable> upstreamDeletionsLoc = new PriorityQueue<>(Comparator.comparingInt(Locatable::getEnd))` — one instance per `GenotypingEngine`, i.e. per traversal, **accumulated across the whole ordered walk** and **never cleared** in production |

So the order is: **subset (`:155`) → ownership test (`:314`, inside it) → record
(`:179`)**, with the ownership test reading a set that can only contain deletions
recorded by strictly earlier loci. Two consequences:
a locus' own deletion can never own its own `*` (already fixed in the previous
round by making the start test strict), and a deletion allele that was pruned at
`:312-316` or whose *record* was dropped at `:167-175` is never recorded at all.

### 1.1 What native computed instead

Before this round (`genotype_gvcf_tool.cpp`, committed HEAD `dce568f`):
* `:1256-1264` `spanning_deletion_supported_at()` — `span.begin < record.pos &&
  record.pos < span.end`, over
* `:1233-1245` `record_has_concrete_deletion()` and a **pre-pass** that pushed
  `VariantSpan{rid, pos, record_span_end(...)}` for every input record carrying a
  concrete deletion allele (`:6750-6755` aggregate, `:5789-5823` streaming probe),
* consumed at the union/merge step (`:6770-6795` aggregate, `:6097-6114`
  streaming), which set `record.orphan_spanning_deletion` and the
  `orphan_spanning_deletion_loci` telemetry counter.

`record_span_end()` is `pos + strlen(REF)` (or `INFO/END`), i.e. the input
record's own end, and the set was built from inputs, unordered and global. All
three measured differences follow from that.

### 1.2 Measured directions (pinned GATK 4.6.2.0, `.diag/emitted-ownership-probe.log`)

| # | fixture | GATK 4.6.2.0 | native BEFORE | native AFTER |
| --- | --- | --- | --- | --- |
| R | `2 AA A,<NON_REF>` implausible + `3 A *,G,<NON_REF>` plausible | `chr1 3 . A G 82.26 . AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11 GT:AD:DP:PL 0/1:0,20:20:0,0,0` | `chr1 3 . A *,G … GT:AD:DP:GQ:PL 1/2:0,0,20:20:99:100,100,100,100,0,100` | byte-identical to GATK |
| O | `2 AAAA AAC,<NON_REF>` (emitted; GATK records 2-3) + `4 A *,G` | `chr1 2 . AAAA AAC 92.60 …` and `chr1 4 . A G 82.26 … 0/1:0,20:20:0,0,0` | row 1 identical, row 2 = `chr1 4 . A *,G … 1/2` | byte-identical to GATK |
| O2 | as O but the `*` at 3 (inside 2-3) | `chr1 3 . A *,G … 1/2 …` | identical | identical (control) |
| T | `1 AA A` + `2 AA *,AC` (emitted `*`, no concrete deletion allele) + `3 A *,G` | 3 rows, `chr1 3 . A *,G … 1/2` | `chr1 3 . A G … 0/1` (**native pruned what GATK keeps**) | byte-identical to GATK |
| T2 | as T but locus 2 carries `*,A` | 3 rows | identical | identical (control) |
| P / P2 | covering upstream deletions (`2 AAAA AAC` + `3 A *,G`; `2 AAA A` + `4 AA *,A`) | keeps `*` | identical | identical (positive controls) |
| S | dense mode `2 AAA A` + `4 AA *,<NON_REF>` + `5 A *,G` | 4 rows (incl. materialized `3 A *`, `4 A *`) | 3 rows; `5 A G` | 3 rows; `5 A *,G` (row now matches GATK; the missing rows are a separate materialization gap, §7) |

R and O and T are the three divergences; O/T were found **by this round's
measurement** (O was not in the previous round's residual, T is the opposite
direction).

Reachability of the whole predicate (read, not assumed): both call sites of the
ownership answer are inside `apply_gatk_output_allele_subset()`, which returns
immediately unless `--gatk-compatible-annotations` is set, so the native
diagnostic profile is structurally untouched.

## 2. STEP 2 — blast radius, assessed before changing anything

**Is it bookkeeping or restructuring? Bookkeeping.** Measured reasons:

1. The ownership answer is consumed **only** in `apply_gatk_output_allele_subset()`
   (three uses: the `owned` test at `:3047`, the orphan-PL normalisation at
   `:3075`, the MLEAF denominator at `:3184` — all in that one function), which
   runs in the compute stage of `ThreeStagePipeline`. `keep_spanning_deletion`
   itself was used at the union step for **nothing but** the flag and the
   telemetry counter (the `*` stays in the merged ALT union either way — see the
   union loops).
2. That compute stage is a **single worker thread consuming a FIFO queue**
   (`fastgatk-runtime/include/fastgatk/runtime/pipeline.hpp:118-137`; the decoder
   and encoder are the other two threads), so per-locus state maintained inside
   it is in traversal order by construction, in **both** the aggregate and the
   `--stream-by-locus` path (the same `apply_gatk_output_allele_subset()` call at
   `:6238` and `:6958`).
3. GATK's own structure puts the test inside the subset function, so this is a
   relocation **towards** the reference implementation, not around it.
4. What is needed is exactly "which deletions did the earlier loci emit", i.e. a
   replay of `recordDeletions()`; no record has to be buffered, no ordering has
   to change, and no stage boundary moves.

Sensitivity of registered tests (named and counted before the change):

* `ctest -R 'genotype-gvcf'` matches **18** tests. Fixtures that contain **both**
  a `*` ALT and a deletion-ish upstream allele — the only shape whose answer can
  change — are three of them:
  * `fastgatk-genotype-gvcf-spanning-deletion-gatk-oracle` — **the most
    sensitive**: it asserts full rows *and*
    `telemetry.orphan_spanning_deletion_loci == (1 if case=='orphan' else 0)`
    **and** that the aggregate and stream-by-locus counts agree
    (`verify_genotype_gvcf_spanning_deletion_gatk_oracle.py:123-126`). The
    telemetry counter therefore had to keep firing exactly once per orphaned
    `*` locus after moving into the compute stage — it does (§5);
  * `fastgatk-genotype-gvcf-lowqual-gatk-oracle` — `STAR_ONLY_COVERED_RECORD`
    (`2 AA A` + `3 A *,<NON_REF>`, deletion emitted) must still be *owned*;
  * `fastgatk-genotype-gvcf-spandel-gatk-oracle` — the gate this round edits.
* Every other `*` fixture in the registered corpus has **no** deletion allele at
  all (`verify_genotype_gvcf.py:386`, `verify_genotype_gvcf_header_order_gatk_oracle.py:160`,
  one `verify_genotype_gvcf_lowqual_gatk_oracle.py:173` case, and the
  `*`-only star cases), so the recorded set is empty either way and no answer can
  move. `grep -rn '\\t\*[,\\t]' fastgatk-native/scripts/*.py` was used to
  enumerate them; the SAM-`*` hits are unrelated fields.
* No C++ unit test references `spanning_deletion_supported_at`,
  `deletion_spans` or `orphan_spanning_deletion` (`grep` over `*.cpp/*.hpp/*.txt/*.cmake`:
  `genotype_gvcf_tool.cpp` only).
* The cheap ones were run **before** the change: the full `genotype-gvcf` ctest
  subset was already green on the committed HEAD (18/18, exit 0, 904.9 s,
  `.diag/emitted-ownership-ctest-baseline.log`).

## 3. STEP 3 — the gate

`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`
(+176/-26; case census 45 = **41 gated** + 4 reported-only, was 41 = 40 + 1):

* `unemitted-upstream-deletion-star-plus-concrete-alt` **promoted REPORTED-ONLY →
  GATED** with the measured GATK row `GATK_G_DOWNSTREAM_ROW`, its `why` rewritten
  around `:155` / `:178-179` / `:314` / `:343-357` / `:365-371`;
* **4 new GATED isolating cases**, each with its measured rows in `expect`:
  * `emitted-deletion-interval-end-is-deletion-size` (direction 2, O) and its
    control `emitted-deletion-interval-end-covered-star` (O2);
  * `emitted-star-records-its-own-interval` (direction 3, T) and its control
    `emitted-star-interval-control` (T2);
* the module docstring's "deletion-ownership cases" section rewritten to state
  the rule, the three native defects and the `EmittedDeletions` implementation.

**Before the fix: exit 1**, `"status": "divergence"`, exactly 3 violations
(`.diag/emitted-ownership-oracle-before.log`):

```
[unemitted-upstream-deletion-star-plus-concrete-alt] row 0 is not byte-identical:
  GATK  ='chr1\t3\t.\tA\tG\t82.26\t.\tAC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11\tGT:AD:DP:PL\t0/1:0,20:20:0,0,0'
  NATIVE='chr1\t3\t.\tA\t*,G\t82.26\t.\tAC=1,1;AF=0.500,0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1,1;MLEAF=0.500,0.500;QD=4.11\tGT:AD:DP:GQ:PL\t1/2:0,0,20:20:99:100,100,100,100,0,100'
[emitted-deletion-interval-end-is-deletion-size] row 1 is not byte-identical:
  GATK  ='chr1\t4\t.\tA\tG\t82.26\t.\t...\tGT:AD:DP:PL\t0/1:0,20:20:0,0,0'
  NATIVE='chr1\t4\t.\tA\t*,G\t82.26\t.\t...\tGT:AD:DP:GQ:PL\t1/2:0,0,20:20:99:100,100,100,100,0,100'
[emitted-star-records-its-own-interval] row 2 is not byte-identical:
  GATK  ='chr1\t3\t.\tA\t*,G\t82.26\t.\t...\tGT:AD:DP:GQ:PL\t1/2:0,0,20:20:99:100,100,100,100,0,100'
  NATIVE='chr1\t3\t.\tA\tG\t82.26\t.\t...\tGT:AD:DP:PL\t0/1:0,20:20:0,0,0'
```

The two controls and every pre-existing gated case were green **before** the
fix, so the three violations isolate exactly the recorded-set rule.

## 4. STEP 4 — the fix (state bookkeeping, implemented)

`fastgatk-native/src/genotype_gvcf_tool.cpp` only (full patch
`.diag/emitted-ownership-fix.diff`, +130/-106; the patch file was captured before the
post-suite comment correction, so its counts read +131/-105):

1. `struct EmittedDeletions` (`:1262-1307`, preceded by the source
   citation comment at `:1233-1261`), replacing
   `record_has_concrete_deletion` / `spanning_deletion_supported_at` /
   `group_has_supported_spanning_deletion` /
   `group_range_has_supported_spanning_deletion`. It mirrors
   `upstreamDeletionsLoc`: a min-heap by interval end, `cull()` = `:344-346`,
   `covers()` = `:365-371` (`begin < pos && pos < end`, strict start),
   `record()` = `:343-357` over the **emitted** allele list with
   `end = pos + (refLen - altLen) + 1` (native's half-open spelling of GATK's
   inclusive `start + deletionSize`). `deletionSize` is invariant under
   `reverseTrimAlleles()` (which only clips trailing bases), so computing the
   interval on native's untrimmed merged record gives GATK's interval.
2. `apply_gatk_output_allele_subset()` (`:2986`) now takes
   `EmittedDeletions& upstream_deletions` and answers the ownership question at
   the top (`:3002-3016`) — exactly where GATK asks it (`:314`) — instead of
   reading a flag computed upstream; `record()` is called at the three points
   where a locus is actually emitted: the no-subset guard (`:3017`), the
   `pruned == 0` path (`:3052`) and the main path (`:3077`, just before the
   Number=G remap, mirroring `:178-179`). The record-destroyed path deliberately
   records nothing (GATK returns `null` at `:167-169`, before `:179`).
3. The union/merge steps no longer compute the flag: the aggregate loop
   (`:6798-6803`) and the streaming producer (`:6135-6140`) keep `*` in the merged
   ALT union and leave the decision to the compute stage; the input-derived
   `deletion_spans` pre-passes were removed from both paths (aggregate
   `:6783-6788`, streaming probe `:5845-5860` and `:5738`).
4. The telemetry counter moved with the decision (incremented once per locus that
   carries a `*` and is not covered), so the counters the registered
   spanning-deletion oracle asserts are unchanged.

No kernel, no ABIs, no `CMakeLists.txt`, no root markdown, no other tool, no
existing test weakened.

## 5. STEP 5 — gate results (mandatory order)

| step | command | result |
| --- | --- | --- |
| a | `python3 fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` | **exit 0**, `"status": "pass"`, 41 gated / 0 violations; **exit 1 with 3 violations before the fix** |
| b1 | `python3 fastgatk-native/scripts/verify_genotype_gvcf.py` | exit 0 — `{"status": "pass", "output_records": 1}` |
| b2 | `ctest --test-dir fastgatk-native/build -R 'genotype-gvcf' -V` | exit 0 — **18/18 passed**, 933.8 s (baseline before the change: 18/18, 904.9 s) |
| c | the 21-pattern strict-gate filter | exit 0 — **23/23 matched tests passed**, 1856.2 s |
| d | `run_regression.sh --label emitted-ownership` | see §5.1 |

(c) matched 23 tests: the 21 patterns plus their siblings —
window-invariance, ploidy-window-invariance, alleles-overlap, span-del-qual,
gvcf-symbolic-prior, arbitrary-ploidy-span-del-prior,
polyploid-gvcf-span-del-prior, spanning-prior-genotype-gq, af-zero-format,
multialt-owner-annotation, mutect2-recheck ×2, gvcf-indel-end, culprit,
select-variants-refonly, asfilterstatus, variant-filtration-flag-only,
droplowqual, spandel, malformed, malformed-input-fail-loud, header-order,
lowqual.

### 5.1 Full double-backend regression

`fastgatk-native/scripts/run_regression.sh --label emitted-ownership` (both
backends in parallel, `FASTGATK_REQUIRE_GATK_ORACLE=1`, ctest parallelism 8):

```
| 后端   | 构建目录                               | 结果 | 通过/总数 | 耗时        |
| omp    | OpenMP (fastgatk-native/build)         | 通过 | 304/304   | 1612.23 sec |
| serial | Serial (fastgatk-native/build-serial)  | 通过 | 304/304   | 1592.66 sec |
```

`REG_EXIT=0`; the runner's staleness check reported no warning for either build
directory. Evidence block
`.diag/regression/20260912-211611/{evidence.md,omp.log,serial.log}` (suite log
`.diag/emitted-ownership-regression.log`).

After the suite, one **comment-only** correction was made in
`genotype_gvcf_tool.cpp` ("Two consequences" → "Three consequences", three
bullets follow) and both trees were rebuilt: the binaries are **byte-identical**
to the ones the suite ran on — omp `88e7bad4a7daaaf51e1b480626fa1444`, serial
`910962353ecdbffb5aa93cfa93db39da` before and after (`md5sum`, also recorded in
`.diag/emitted-ownership-binaries.md5`) — so the suite was **not** re-run. No
behavioural, test or script change was made after the suite.

Logs: `.diag/emitted-ownership-oracle-before.log`,
`.diag/emitted-ownership-oracle-after.log`,
`.diag/emitted-ownership-ctest-baseline.log`,
`.diag/emitted-ownership-ctest-genotype-gvcf.log`,
`.diag/emitted-ownership-strict-gates.log`,
`.diag/emitted-ownership-verify-gvcf.log`,
`.diag/emitted-ownership-regression.log`.

### 5.2 Extra (scratch) verification of the `--stream-by-locus` path

The new gated cases run the aggregate traversal only. `.diag/emitted_ownership_stream_probe.py`
runs R, O, T and the P2 control through GATK, native aggregate and native
`--stream-by-locus`: all three are byte-identical for all four fixtures
(`{"status": "pass", "failures": 0}`, `.diag/emitted-ownership-stream-probe.log`).
The registered oracle that does assert both traversals (rows **and** telemetry)
is `fastgatk-genotype-gvcf-spanning-deletion-gatk-oracle`, green in (b2).

## 6. Tree state

* Change **kept in the tree** (no commit, no branch).
* Modified: `fastgatk-native/src/genotype_gvcf_tool.cpp`,
  `fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`. Scratch
  (`.diag/`, git-ignored) only otherwise.
* Not touched: `fastgatk-native/CMakeLists.txt`, root `*.md`, `third_party/`,
  Mutect2, every other tool, and every registered test script other than the
  gate this round owns.
* Both trees were rebuilt after the last **behavioural** source edit and before
  every gate above (omp `build/fastgatk-genotype-gvcf`, serial
  `build-serial/fastgatk-genotype-gvcf`); the oracle script was frozen before
  them. Binary digests before/after the post-suite comment correction:
  `.diag/emitted-ownership-binaries.md5` (§5.1).
* `git status --short`: exactly the two modified files above and nothing else
  (`.diag/` is git-ignored, so the report, probes and logs do not appear).
* `git diff --stat`: `genotype_gvcf_tool.cpp` +130/-106 (comments included),
  `verify_genotype_gvcf_spandel_gatk_oracle.py` +176/-26.

## 7. What remains unproven / open

* **Dense-mode record materialization is still divergent** (measured, not
  gated): with `--include-non-variant-sites` GATK emits a spanning-deletion row
  for positions an emitted deletion covers (`2 AAA A` + `4 AA *,<NON_REF>` +
  `5 A *,G` → GATK 4 rows incl. `3 A *` and `4 A *`, native 3) and reverse-trims
  merged records (`4 AA *` native vs `4 A *` GATK). This is a
  record-materialization gap, not ownership: it is identical before and after
  this round, except that the `*` at 5 now matches GATK. NOTE: *speculation* —
  fixing it is a different, larger change.
* The new gated cases exercise the aggregate traversal; the streaming traversal
  is verified by the scratch probe in §5.2 and by the registered oracle, but no
  **gated** streaming case was added for R/O/T.
* Fixtures with **two records starting at one locus** (the cross-sample
  overlapping-record shape, already REPORTED ONLY) are unmeasured for this rule;
  in that shape native drops the overlapping record's alleles before the compute
  stage, so a deletion allele contributed only by the dropped record could still
  be missing from the state. *Speculation*: unreachable through the documented
  pipeline (GATK rejects or splits that input), but not demonstrated.
* Multi-shard/`CombineGVCFs`-produced inputs were not exercised for this rule;
  all fixtures here are single-input, single-sample except where noted.
* GATK's own culling policy differs in bookkeeping detail from `cull()` (a
  `PriorityQueue` keyed by end vs. a min-heap with the same key); the coverage
  predicate is identical, and culling cannot change an answer (a culled entry has
  `end <= pos` or sits on another contig), which is argued, not measured.
* Pre-existing, untouched items from earlier rounds: the
  `isSpanningDeletionOnly` row for star-only loci (`covered-star-only-record`),
  `FILTER=LowQual` coverage, `QUAL=Infinity` rendering, `-all-sites` alias, the
  `genotypeLikelihoods == null` arm of PREFER_PLS, reverse trimming of merged
  records, and the cross-sample same-start-record shape.
