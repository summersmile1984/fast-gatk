# Round: GenotypeGVCFs' own `FILTER=LowQual` (GenotypingEngine's call-threshold filter)

Scope: the residual the two previous rounds measured and deliberately left ungated
(`.diag/round-filter-rgq.md` §7 "GATK's own `LowQual` (residual, reported)" and
`.diag/round-header-filter.md` §7).  On the weak-locus fixture pinned GATK 4.6.2.0 writes

```
chr1	2	.	A	.	23.14	LowQual	DP=20;MLEAC=.;MLEAF=.	GT	./.
```

while native wrote `FILTER=.` in the same column.  Characterise first, then blast
radius, then gate, then fix.

## 0. Verdict

* **GATK's rule (measured, pinned 4.6.2.0):** `GenotypingEngine.calculateGenotypes()`
  recomputes a phred-scaled site confidence from the AF result

  ```java
  final double log10Confidence =
              !outputAlternativeAlleles.siteIsMonomorphic || configuration.annotateAllSitesWithPLs
                      ? AFresult.log10ProbOnlyRefAlleleExists() + 0.0 : AFresult.log10ProbVariantPresent() + 0.0;
  final double phredScaledConfidence = (-10.0 * log10Confidence) + 0.0;
  ```
  (`gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/genotyper/GenotypingEngine.java:158-163`),
  rebuilds the record from scratch and applies its **own** filter to it:

  ```java
  final VariantContextBuilder builder = new VariantContextBuilder(callSourceString(),
          vc.getContig(), vc.getStart(), vc.getEnd(), outputAlleles);   // :181
  builder.log10PError(log10Confidence);                                  // :183
  if ( ! passesCallThreshold(phredScaledConfidence) ) {                  // :184
      builder.filter(GATKVCFConstants.LOW_QUAL_FILTER_NAME);             // :185
  }
  ```
  with `passesCallThreshold(conf) = conf >= configuration.genotypeArgs.standardConfidenceForCalling`
  (`:430-432`) and `LOW_QUAL_FILTER_NAME = "LowQual"` (`GATKVCFConstants.java:179`).
  The header declaration is unconditional (`GenotypeGVCFsEngine.java:416`, text
  `GATKVCFHeaderLines.java:89`).
* **The value tested is the raw recomputed double, not the published QUAL token.**
  Measured by bracketing: on the weak-locus fixture GATK applies `LowQual` at
  `--standard-min-confidence-threshold-for-calling 23.1358` but not at `23.135`, so
  the tested confidence lies inside `[23.135, 23.1358)` while the QUAL column reads
  `23.14` (htsjdk rounds QUAL to two decimals).  An implementation that compares the
  rounded token is therefore wrong, and the gate's `near-cutoff-*` pair fails it.
  The tested confidence is also the **pre-posterior** one: `:184-186` runs before the
  optional `--use-posteriors-to-calculate-qual` update at `:192-199`.
* **The filter is site-level, and the per-allele intersection rule does not apply
  here.**  `builder.filter()` takes no allele; GenotypeGVCFs never calls
  `AlleleFilterUtils.addAlleleAndSiteFilters()` (whose intersection at
  `AlleleFilterUtils.java:115-120` belongs to VariantFiltration/Mutect2), and
  `--invalidate-previous-filters` is not a GenotypeGVCFs option: measured, GATK exits
  1 with `invalidate-previous-filters is not a recognized option` and native exits 2
  with `unknown option`.
* **Reach: dense mode only.**  `passesEmitThreshold()` (`:425-427`) requires
  `passesCallThreshold()` before a best-guess-reference locus is emitted at all, so in
  every non-dense fixture measured a below-threshold locus produces **no record** and
  the FILTER question never arises.  The divergence is confined to
  `--include-non-variant-sites`.
* **Fix kept in the tree** (confined, 1 production file): the raw confidence is
  recorded by the two engine stages that already derive QUAL, the LowQual id is
  declared in the output `bcf_hdr_t`, and the filter is applied in
  `apply_gatk_annotation_compatibility()` *after* the inherited-filter clear (which
  would otherwise wipe it).  107 insertions / 2 deletions (≈25 lines of code, the
  rest comments) in `fastgatk-native/src/genotype_gvcf_tool.cpp`.
* **New gate** `fastgatk-native/scripts/verify_genotype_gvcf_lowqual_gatk_oracle.py`:
  **exit 1 before** the fix (7 FILTER-only violations), **exit 0 after** (14 cases,
  0 violations).  Written and run BEFORE the fix, as required.
* Step 5: (a) exit 0; (b1) exit 0, (b2) 17/17 exit 0; (c) 22/22 exit 0;
  (d) **303/303 on both backends**, `REG_EXIT=0`.  Change **kept in the tree**.

## 1. STEP 1 — measured rule, case by case

Probes: `.diag/lowqual_probe.py` (log/JSON `.diag/lowqual-probe.log`) and
`.diag/lowqual_probe2.py` (`.diag/lowqual-probe2.log`).  Both build the 100 bp `chr1`
reference used by the sibling oracles, write one gVCF, index it with pinned
`IndexFeatureFile`, and run

```
GATK:   java -jar gatk-package-4.6.2.0-local.jar GenotypeGVCFs -R ref.fa -V in.g.vcf [args] -O out.vcf
native: fastgatk-genotype-gvcf -R ref.fa -V in.g.vcf [args] --gatk-compatible-annotations -O out.vcf
```

Fixtures (all single-sample, `REF=A` at `chr1:2`):

| name | body | why it was chosen |
| --- | --- | --- |
| weak locus | `chr1 2 . A G,<NON_REF> . PASS DP=20  GT:DP:AD:PL 0/1:20:19,1:0,0,40` | below any sane cutoff; its ALT is pruned, so the dense row is REF-only *and* monomorphic |
| strong monomorphic | `… PL 0/1:20:0,20:0,100,100` | the same materialization with a high confidence (control) |
| star-only | `AA→A,<NON_REF>` at 2 + `A→*,<NON_REF>` at 3 | an infinite confidence (`QUAL=Infinity`) |
| reference block | `chr1 2 . A <NON_REF> . PASS DP=20  GT:DP:GQ:PL 0/0:20:99:0,100,100` | a passthrough row the engine never regenotypes |
| two-ALT | `A→G,T,<NON_REF>` with G strong and T weak | per-allele vs site-level |
| G-plausible | `A→*,G,<NON_REF>` (the sibling oracle's fixture) | a surviving concrete ALT |

### 1.1 Literal rows (`QUAL` / `FILTER`, both sides)

| case (args) | GATK | native |
| --- | --- | --- |
| weak, `-all-sites`, cutoff 30 | `23.14` **`LowQual`** | `23.14` `.` |
| weak, `-all-sites`, cutoff 10 | `23.14` `.` | `23.14` `.` |
| weak, `-all-sites`, cutoff 23.131 / 23.135 | `23.14` `.` | `23.14` `.` |
| weak, `-all-sites`, cutoff **23.1358** / 23.136 / 23.14 / 23.15 / 24 / 100 | `23.14` **`LowQual`** | `23.14` `.` |
| weak, `-all-sites`, `--use-posteriors-to-calculate-qual` (GP-bearing header) | `23.14` **`LowQual`** | `23.14` `.` |
| weak, **no** `-all-sites`, cutoff 30 or 100 | *no records* | *no records* |
| strong monomorphic, `-all-sites`, cutoff 30 | `127.78` `.` | identical |
| star-only, `-all-sites`, cutoff 30 | `2: 92.60 .` and `3: Infinity .` | identical |
| star-only, `-all-sites`, cutoff **1e9** | `2: 0` **`LowQual`**, `3: Infinity` `.` | `2: 0` `.`, `3: Infinity` `.` |
| reference block, `-all-sites`, cutoff 30 | `QUAL='.' FILTER='.'` | identical |
| two-ALT, `-all-sites`, cutoff 30 | `25.59` **`LowQual`** | `25.59` `.` |
| G-plausible, `-all-sites`, cutoff 100 | `0` **`LowQual`** | `0` `.` |
| weak, `-all-sites`, `--invalidate-previous-filters` | GATK exit **1**, `is not a recognized option` | native exit **2**, `unknown option` |
| weak, `-all-sites`, cutoff **-5** | GATK exit **3** (invalid argument) | native exit **2** |

Every `LowQual` above is on a **monomorphic, REF-only** dense row, i.e. the filter is
applied to the record that `materialize_gatk_monomorphic_ref_call()` builds, not to a
variant row.  The one dense case that keeps a concrete ALT (`G-plausible`, cutoff 100)
prunes to REF-only as well, so the surviving-ALT variant of the question was answered
by reading the code rather than by measurement: `:155` computes the output-allele
subset *before* `:181`, and `:184-186` filters the site, so the decision is
independent of which alleles survived (this is what
`surviving-alt-below-threshold` pins).

### 1.2 What each question resolves to

* **Which confidence?** `phredScaledConfidence` = `-10 * log10Confidence`
  (`:163`), with `log10Confidence = log10ProbOnlyRefAlleleExists` for a variant site
  and `log10ProbVariantPresent` for a monomorphic one (`:158-163`).  For a
  monomorphic locus that is *not* the published QUAL formula used elsewhere: it is
  `-10*log10(1 - p_no_variant)`, which is exactly what native already computes at
  `genotype_gvcf_tool.cpp:2878-2880` (`gatk_log10_one_minus_pow10`), so the two sides
  agree bit-for-bit at the bracket: native reproduces both halves of the
  `[23.135, 23.1358)` measurement.
* **How does QUAL relate to it?**  QUAL is the *same* double
  (`builder.log10PError(log10Confidence)`, `:183`) rendered by htsjdk with two
  decimals and a stripped `.00`; the weak-locus row prints `23.14` for a confidence
  that is strictly below `23.14`.  `--use-posteriors-to-calculate-qual` can replace
  QUAL afterwards (`:192-199`) **without** changing the filter decision.
* **Dense / REF-only / monomorphic rows?**  Yes — they are the only place the rule is
  observable in practice, and they get `LowQual` whenever the confidence is below the
  cutoff.  A *reference-block* row (never regenotyped, QUAL missing) is **not**
  filtered.
* **Intersection rule / `--invalidate-previous-filters`?**  Not part of
  GenotypeGVCFs: the tool has no such argument (measured rejection, §1.1) and never
  calls `addAlleleAndSiteFilters()`.  `LowQual` participates at the **site** level
  only.

## 2. STEP 2 — blast radius, measured before touching anything

### 2.1 Which native outputs change

Only records (i) written through the GATK-compatibility writer, (ii) for which the
genotyping engine computed a confidence, and (iii) whose confidence is below the
cutoff.  Condition (ii) excludes reference-block passthrough rows; condition (iii) is
unreachable outside `--include-non-variant-sites` because `passesEmitThreshold()`
(`:425-427`) already requires the cutoff before such a locus is emitted (§1.1, row
"weak, no `-all-sites`" shows zero records on both sides).

* **Census of expected rows.**  `.diag/lowqual_row_census.py` scans all 335
  `fastgatk-native/scripts/*.py` for literal VCF rows: 405 rows, of which **9** have
  `QUAL < 30`.  Eight belong to other tools (`verify_select_variants.py`,
  `verify_variant_filtration.py`, `verify_variants_to_table.py`).  Exactly **one** is a
  GenotypeGVCFs output expectation — `verify_genotype_gvcf_spandel_gatk_oracle.py:875`,
  whose FILTER column already reads `LowQual` (it is that oracle's REPORTED ONLY
  case).  So no registered expectation has a low-`QUAL` genotype row carrying `.`.
* **Rows that flip on the existing fixtures: 1** — the spandel oracle's reported-only
  row (`chr1 2 . A . 23.14 .` → `… LowQual …`).  Every other gated fixture keeps its
  records byte-identical: the 18 FILTER-sensitive tests of §2.2 and the whole 303-test
  suite pass after the change (§5), which is the empirical bound on "no green row
  flipped".  Re-run standalone, that case now reports `"violations": []` and an
  unchanged header (`only_in_gatk: []`, `only_in_native: []`, `order_differs: false`,
  `.diag/lowqual-spandel-reported-case.log`).

### 2.2 FILTER-sensitive registered tests (named, counted)

20 registered tests exercise the genotype binary.  **18 of them are FILTER-sensitive**:
**17 compare whole GATK-derived records or whole files** (the FILTER column is
compared implicitly), and **1 compares the whole header including the `##FILTER`
group**:

| # | ctest name | comparison |
| --- | --- | --- |
| 1 | `fastgatk-genotype-gvcf-spandel-gatk-oracle` | every gated case's row vs GATK, byte for byte |
| 2 | `fastgatk-genotype-gvcf-malformed-gatk-oracle` | `GOOD_ROW` full row after the malformed-input contract |
| 3 | `fastgatk-genotype-gvcf-contract` (`verify_genotype_gvcf.py`) | 28 row/field assertions incl. the dense star row |
| 4 | `fastgatk-genotype-gvcf-gatk-oracle` (`verify_gatk_genotype_gvcf.py`) | record text + semantic header set (incl. `##FILTER`) |
| 5 | `fastgatk-genotype-gvcf-legacy-qual-gatk-oracle` | `native_rows == gatk_rows`, header equality |
| 6 | `fastgatk-genotype-gvcf-multisample-gatk-oracle` | `native_rows == gatk_rows`, header equality |
| 7 | `fastgatk-genotype-gvcf-multiallelic-gatk-oracle` | `native_rows == gatk_rows` |
| 8 | `fastgatk-genotype-gvcf-inbreeding-gatk-oracle` | `native_rows == gatk_rows` |
| 9 | `fastgatk-genotype-gvcf-include-non-variant-gatk-oracle` | `body(native) == gatk_rows`, 101 dense rows ×2 writers |
| 10 | `fastgatk-genotype-gvcf-assignment-gatk-oracle` | `gatk_rows != native_rows` check |
| 11 | `fastgatk-genotype-gvcf-multisample-reference-confidence-gatk-oracle` | 101 dense rows, aggregate + stream |
| 12 | `fastgatk-genotype-gvcf-max-alternate-alleles-gatk-oracle` | single-row field assertions |
| 13 | `fastgatk-genotype-gvcf-gp-input-gatk-oracle` | `native_rows == gatk_rows` |
| 14 | `fastgatk-genotype-gvcf-exclude-intervals-gatk-oracle` | `expected == aggregate == stream` |
| 15 | `fastgatk-genotype-gvcf-starts-in-intervals-gatk-oracle` | `observed == expected` |
| 16 | `fastgatk-genotype-gvcf-spanning-deletion-gatk-oracle` | `actual == expected` |
| 17 | `fastgatk-hc-dense-gvcf-genotype-gatk-oracle` | whole joint VCF text equality (rows + header) |
| 18 | `fastgatk-genotype-gvcf-header-order-gatk-oracle` | every header line in order (header only) |

Not sensitive: `fastgatk-hc-genotype-priors-contract` and
`fastgatk-hc-spanning-prior-genotype-gq-oracle` (HaplotypeCaller side).
`grep -n '\[6\]'` over all genotype scripts returns **0** explicit FILTER-column
assertions: the sensitivity is entirely implicit, via whole-row comparison — which is
why the divergence survived this long and why the list above, not a FILTER grep, is the
right blast-radius measure.

**Assessment: confined → proceed.**  The change is one file, one tool, one output
profile; the sensitive list is long but *explainable* (all 18 are whole-row
comparisons against pinned GATK, i.e. exactly the tests that should see a FILTER
change), and the census says no registered expectation pins the old `.` on a
below-threshold row.  Step 5 then runs all of them.

## 3. STEP 3 — the gate, written and run BEFORE the fix

`fastgatk-native/scripts/verify_genotype_gvcf_lowqual_gatk_oracle.py` (house
conventions: `main() -> int`, `--native`/`$FASTGATK_GENOTYPE_BINARY`, pinned
GATK/JDK17 from `third_party/`, `TemporaryDirectory` scratch, single-line JSON status
payload, `--expect-divergence`; it prints `oracle_not_verified` when the oracle is
absent and honours `FASTGATK_REQUIRE_GATK_ORACLE`).  It is **not registered** by me —
`fastgatk-native/CMakeLists.txt` must not be edited this round; the orchestrator
registers it.

14 cases, each one GATK run and one native run on the same input and arguments.  Every
case pins **GATK's own rows as literals** and then compares native to GATK, so the gate
fails if either side moves; the `##FILTER` header group is compared as well:

`below-threshold-dense-monomorphic`, `above-threshold-dense-monomorphic`,
`near-cutoff-just-above` (cutoff 23.135), `near-cutoff-just-below` (23.1358),
`threshold-lowered-below-confidence` (10), `threshold-raised-above-confidence` (100),
`surviving-alt-below-threshold`, `two-alt-locus-below-threshold`,
`infinite-confidence-never-filtered` (cutoff 1e9), `reference-block-row-not-filtered`,
`posterior-qual-option-keeps-call-confidence`,
`default-mode-emits-no-low-confidence-row`, `default-mode-strong-variant-unfiltered`,
`invalidate-previous-filters-unsupported`.

**BEFORE the fix: exit 1**, `"status": "divergence"`, **7 violations**, every one the
FILTER column of a dense row (`.diag/lowqual-gate-before.log`) — no GATK-side
"expected rows moved" violation, so every pinned literal reproduced:

```
[below-threshold-dense-monomorphic]   GATK 23.14 LowQual | native 23.14 .
[near-cutoff-just-below]              GATK 23.14 LowQual | native 23.14 .
[threshold-raised-above-confidence]   GATK 23.14 LowQual | native 23.14 .
[surviving-alt-below-threshold]       GATK 0     LowQual | native 0     .
[two-alt-locus-below-threshold]       GATK 25.59 LowQual | native 25.59 .
[infinite-confidence-never-filtered]  GATK 0     LowQual | native 0     .
[posterior-qual-option-keeps-call-confidence] GATK 23.14 LowQual | native 23.14 .
```

**AFTER the fix: exit 0**, `"status": "pass"`, 0 violations
(`.diag/lowqual-gate-after.log`).

The gate also caught an incomplete first attempt, which is the reason it is worth
having: applying the filter without declaring the id in the aggregate writer's
`bcf_hdr_t` made native exit 2 with
`OUTPUT_CONTRACT_FAILURE: cannot write LowQual FILTER: the id is not declared in the
output header` on all seven cases — the fail-loud path, not a silently wrong token.

## 4. STEP 4 — the fix

`git diff --stat` = **1 file changed, 107 insertions(+), 2 deletions(-)**
(`fastgatk-native/src/genotype_gvcf_tool.cpp`; the new gate script is untracked).
Six code sites, all inside the `--gatk-compatible-annotations` boundary:

1. **`struct Record`** — two host-side fields next to the existing cohort state:
   `double call_confidence` and `bool call_confidence_available`, with the `file:line`
   derivation of `phredScaledConfidence` and the reason the **raw, pre-posterior**
   double is kept (the `23.14`-token / `[23.135, 23.1358)` measurement) in the
   comment.  They travel with the record through the existing
   `GenotypeDecoded → GenotypeComputed → GenotypeEncoded` moves; no serialization
   format is affected (the stage queues are in-process, `genotype_record_bytes()` is
   only a capacity hint).
2. **`update_cohort_af_annotations()`** (variant branch) — records
   `call_confidence = -10.0 * cohort.log10_p_no_variant`, the same quantity it already
   publishes as QUAL one line above (`= AFresult.log10ProbOnlyRefAlleleExists`, and
   `-10 * log10ProbOnlyRefAlleleExists` is exactly `:163` for a non-monomorphic site).
3. **`materialize_gatk_monomorphic_ref_call()`** — records `-10.0 * complement` for
   the finite branch and `+Infinity` for the `complement == -Infinity` branch
   (`QUAL=Infinity` passes every finite cutoff, which is why the star-only row is never
   filtered).  The NaN branch keeps the previous leave-unchanged behaviour and is
   **not** recorded, so that degenerate shape stays exactly where it was measured.
4. **`ensure_gatk_low_qual_filter_declaration()`** (new, next to
   `kGatkLowQualFilterLine`) — appends GATK's unconditional `##FILTER` line to the
   output **`bcf_hdr_t`** when the compat profile is on and the id is absent, called
   from both header-setup sites (the streaming one and the aggregate one, which builds
   its header inline rather than through `add_genotype_output_header_fields`).  The
   header *text* is unaffected: `gatk_compatible_header_text()` collects the group from
   whatever the header holds, de-duplicates against this exact line, sorts it and
   re-inserts it where the input's filter lines were — and HTSlib resolves the id we
   write against this header, which is why the line has to be there at all.
5. **`apply_gatk_annotation_compatibility()`** — the application point, *after* the
   inherited-filter clear (the clear at `:4258-4261` would otherwise wipe a newly
   applied LowQual, as the previous round's finding warned):

   ```cpp
   if (record.call_confidence_available &&
       !(record.call_confidence >= options.standard_confidence_for_calling)) {
       int low_qual = bcf_hdr_id2int(output_header, BCF_DT_ID, "LowQual");
       if (low_qual < 0) throw std::runtime_error(/* fail loud */);
       if (bcf_update_filter(output_header, record.value, &low_qual, 1) != 0)
           throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write LowQual FILTER");
   }
   ```

   The `!(a >= b)` spelling is GATK's `!passesCallThreshold(...)` verbatim
   (`:430-432`), so NaN behaves as in Java.

Not restructured: the genotyping engine and its Kokkos kernels (only the existing
`cohort.qual` / `cohort_log10_p_no_variant` outputs are consumed), the kernel API/ABI,
the AF/PL math, the output-allele subset, the inherited-filter clear, every other tool,
`fastgatk-native/CMakeLists.txt` and every root `*.md`.  The non-compatibility native
diagnostic profile is untouched by construction — both new calls are behind
`options.gatk_annotation_compatibility`.

### 4.1 Collateral, explicitly

* `Record` grows by 16 bytes; the two copy sites (`Record split = original;`) are plain
  struct copies, so the fields follow.
* The stub `if (std::isfinite(complement))` became a braced `if/else if` (the 2
  deletions) — no logic change in the existing branches.
* Failed first attempt (kept in the log as evidence): the aggregate writer does not use
  `add_genotype_output_header_fields()`, so the id was missing there; the fail-loud
  check caught it before any wrong output was produced.

## 5. STEP 5 — gate results (mandatory order)

| Step | Command | Result |
| --- | --- | --- |
| a | `python3 fastgatk-native/scripts/verify_genotype_gvcf_lowqual_gatk_oracle.py` | **exit 0**, `"status": "pass"`, 14 cases, 0 violations (**exit 1 / 7 violations before** the fix) |
| b1 | `python3 fastgatk-native/scripts/verify_genotype_gvcf.py` | `B1_EXIT=0`, `{"status": "pass", "output_records": 1}` (`.diag/lowqual-verify-gvcf.log`) |
| b2 | `ctest --test-dir fastgatk-native/build -R 'genotype-gvcf' -V` | `B2_EXIT=0`, **17/17 passed**, 713.68 s (`.diag/lowqual-ctest-genotype-gvcf.log`) |
| c | the 20-name strict-gate filter from the task | `C_EXIT=0`, **22/22 matched tests passed**, 1610.29 s, incl. `fastgatk-genotype-gvcf-spandel-gatk-oracle` (391.03 s) and `fastgatk-genotype-gvcf-header-order-gatk-oracle` (`.diag/lowqual-strict-gates.log`) |
| d | `fastgatk-native/scripts/run_regression.sh --label lowqual-filter` | see below |

(d) both trees were rebuilt (`cmake --build fastgatk-native/build|build-serial --target
fastgatk-genotype-gvcf`, both exit 0, logs
`.diag/lowqual-build-{omp,serial}.log`) before the suite.  Result:

```
REG_EXIT=0
| 后端   | 构建目录                              | 结果 | 通过/总数 | 耗时        |
| omp    | OpenMP (fastgatk-native/build)        | 通过 | 303/303   | 1528.42 sec |
| serial | Serial (fastgatk-native/build-serial) | 通过 | 303/303   | 1399.13 sec |
```

Evidence block `.diag/regression/20260912-191608/{omp.log,serial.log,evidence.md}`
(+ `.diag/lowqual-regression.out` for `REG_EXIT`), run window 19:16:08 → 19:41:36,
label `lowqual-filter`, git `c1cc415` + the two uncommitted files of §7, and the
runner's staleness check reported no warning (`grep -c 警告` = 0).

Additionally, the new gate was re-run against the **serial** binary:
**exit 0**, `"status": "pass"`, 0 violations (`.diag/lowqual-gate-after-serial.log`).

**Nothing was edited after the suite run.**  `md5sum -c .diag/lowqual-binaries.md5`
re-confirms

```
93a5e453786aea89d01ca96672168597  fastgatk-native/build/fastgatk-genotype-gvcf        (omp)
08c1c2aa81a5525548711e6e942e09fa  fastgatk-native/build-serial/fastgatk-genotype-gvcf (serial)
db4f3412dbd1fd753648dfb4086cbd50  fastgatk-native/src/genotype_gvcf_tool.cpp
8fa6d9747581cafc8671ff6608685c1c  fastgatk-native/scripts/verify_genotype_gvcf_lowqual_gatk_oracle.py
```

and the hashes were recorded *before* the suite window.  The only file touched
afterwards is this report, which is not a test input (no source, test script, build
input or fixture changed), so no rebuild and no re-run were needed.

## 6. What is deliberately left, and what remains unproven

* **Legacy AF calculator (`--use-new-qual-calculator false`), dense mode — measured,
  NOT fixed.**  Native's cohort stage returns early without the new calculator, so no
  confidence is recorded and no LowQual is applied.  Measured on the weak locus with
  `-all-sites`:

  ```
  GATK   chr1 2 . A . 23.14 LowQual DP=20;MLEAC=.;MLEAF=. GT ./.
  native chr1 2 . A G .     .       AC=0;AF=0.00;AN=2;DP=20;ExcessHet=0.0000 GT:AD:DP:GQ:PL 0/0:19,1:20:0:0,0,40
  ```

  The row *shape* already diverges there (GATK materializes REF-only, native keeps the
  ALT), so the missing filter is a consequence of a larger, pre-existing legacy-path
  gap and cannot be repaired by this round's change.  The registered legacy-qual gate
  runs the non-dense mode, where GATK emits no below-threshold record at all.
* **Non-compatibility diagnostic profile.**  Both new calls are behind
  `gatk_annotation_compatibility`, so that profile still writes `.` where GATK would
  write `LowQual`.  It is explicitly not a parity surface (§7.1 of the filter-round
  report) and is left as it was.
* **`--force-output-intervals`** (GATK can force a below-threshold site into the
  output, `GenotypeGVCFs.java:324-326`) is not implemented by native at all, so the
  rule is untested on that path.
* **The rounded-token question is pinned, not proved exhaustively.**  The bracket
  `[23.135, 23.1358)` is one measurement on one fixture; it establishes that GATK
  tests a value strictly below the printed token, but I did not instrument the JVM to
  print the double.  Native reproduces both halves of the bracket, i.e. its recorded
  confidence is in the same interval — evidence that the two engines agree far beyond
  the two decimals they publish, not a proof of bit equality.
* **Untested shapes.**  `--dbsnp`, GenomicsDB inputs, `.bcf`, multi-shard inputs,
  `--keep-combined-raw-annotations`, `--force-output-intervals`,
  `USE_POSTERIORS_ANNOTATION`/other assignment modes combined with a below-threshold
  dense row, triploid/polyploid dense rows, and a source `##FILTER=<ID=LowQual,…>`
  declared with a different `Description` (which would make native's header carry two
  `LowQual` declarations — pre-existing, from the header round, not touched here).
  *Speculation:* the rule is a scalar comparison at the writer boundary, so those
  shapes should follow, but that is not measured.
* **Whether "no green row flipped" is exhaustive** rests on the 18 FILTER-sensitive
  tests plus the 303-test suite, not on a proof: a test could build a below-threshold
  fixture dynamically at runtime.  The suite result is evidence, not proof.
* **The NaN-complement branch** (a posterior with `log10Confidence > 0`, i.e.
  `p(variant) > 1`) is unreachable for a real posterior and left as it was; GATK would
  hand NaN to both `log10PError` and `passesCallThreshold`, so it would print `NaN`
  with `LowQual`.  Not measured, not changed.

## 7. Tree state

* Change **kept in the tree** (no commit, no branch, no branch switch).
* Modified: `fastgatk-native/src/genotype_gvcf_tool.cpp` (+107/-2).  Added:
  `fastgatk-native/scripts/verify_genotype_gvcf_lowqual_gatk_oracle.py` (untracked; the
  orchestrator registers it).  Nothing else in the tree — no existing test was edited,
  and no stale assertion was found to correct.
* Not touched: `fastgatk-native/CMakeLists.txt`, every root `*.md`, `third_party/`,
  Mutect2, HaplotypeCaller, every other tool, the Kokkos kernels and their ABI.
* Scratch/evidence under `.diag/` only (git-ignored): probes
  `.diag/lowqual_probe.py`, `.diag/lowqual_probe2.py`, `.diag/lowqual_row_census.py`
  (+ logs/JSON `.diag/lowqual-probe*.log`), gate logs
  `.diag/lowqual-gate-{before,after}.log`, `.diag/lowqual-spandel-reported-case.log`,
  build logs `.diag/lowqual-build-{omp,serial}.log`, step-5 logs
  `.diag/lowqual-{verify-gvcf,ctest-genotype-gvcf,strict-gates}.log`,
  `.diag/lowqual-regression.out`, binary hashes `.diag/lowqual-binaries.md5`.
