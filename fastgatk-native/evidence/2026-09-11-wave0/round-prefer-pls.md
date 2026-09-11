# Round: orphan spanning-deletion `*` with a surviving concrete ALT (PREFER_PLS subset fallback)

Scope: the case the previous round deliberately left **reported only** in
`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`
(`surviving-concrete-alt`), i.e. a gVCF record whose ALT is `*,G,<NON_REF>` where
`G` best supports the site: GenotypeGVCFs must prune the orphan `*` and still
publish a call for the surviving concrete ALT.

## 0. Verdict

* The previous round's residual description ("a **PREFER_PLS genotype
  re-derivation** difference after the orphan `*` is pruned") was **half right
  and half wrong**, and the wrong half was the actionable half:
  * right: the code path is
    `GATKVariantContextUtils.makeGenotypeCall()`'s PREFER_PLS branch;
  * wrong: the genotype is **not** re-derived from the projected PLs.  The PL
    row projected onto the kept alleles is *constant* (`[100,100,100]` for this
    fixture), so `isInformative()` is false and GATK takes PREFER_PLS's **second
    arm**: it ignores the PLs and projects the **source genotype** with
    `bestMatchToOriginalGT()`.  The genotype it publishes (`0/1`) is the input
    call, not a likelihood argmax.  The previous round's note was explicitly
    labelled a hypothesis; it is now measured and refuted.
* Consequence: the fix is not "re-derive the call differently" but "carry the
  merged **source** genotype through the output-allele subset and project it by
  allele identity".  Native had destroyed that call at the union stage, which
  made the naive version of the fix order-sensitive; see §3.
* Fix applied and kept in the tree: `fastgatk-native/src/genotype_gvcf_tool.cpp`
  (+~200 lines, one translation unit) plus promotion of the oracle case and
  three new gated cases.  No kernel, no `CMakeLists.txt`, no Mutect2, no other
  tool, no registered test script edited.
* Oracle: **exit 1 before** the fix (`status: divergence`), **exit 0 after**.
* Step 5: (a) 0, (b) 0 + 15/15, (c) 0 + 19/19, (d) see §6.

## 1. STEP 1 — measured GATK truth (literal rows)

### 1.1 Fixture and command

The fixture is the one already carried by the oracle
(`G_PLAUSIBLE_RECORD`, `verify_genotype_gvcf_spandel_gatk_oracle.py:97-100`):

```
chr1	2	.	A	*,G,<NON_REF>	.	PASS	DP=20	GT:DP:AD:PL	0/2:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100
```

Reference: a generated 100 bp `chr1` of `A` with hand-written `.fai`/`.dict`;
input is a plain `.vcf` indexed by `gatk IndexFeatureFile` (both tools get
byte-identical inputs).  Command:

```
third_party/jdk17/bin/java -Xmx1g -jar \
  third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
  GenotypeGVCFs -R <chr1.fa> -V <fixture>.g.vcf -O out.vcf \
  --create-output-variant-index false
```

### 1.2 Literal row (GATK 4.6.2.0)

```
chr1	2	.	A	G	82.26	.	AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11	GT:AD:DP:PL	0/1:0,20:20:0,0,0
```

Before the fix native published, for the same input:

```
chr1	2	.	A	G	82.26	.	AC=0;AF=0.00;AN=0;DP=20;MLEAC=1;MLEAF=0.00	GT:AD:DP:GQ:PL	./.:0,20:20:0:0,0,0
```

Six columns differ: `INFO` (missing `ExcessHet`/`QD`, wrong `AC`/`AF`/`AN`/
`MLEAF`) and `FORMAT`/sample (`GT`, plus an extra `GQ`).  Note that the **PL
vector and QUAL already agreed** (`0,0,0` / `82.26`) before the fix.

### 1.3 Supporting measurements (all with pinned GATK 4.6.2.0)

Allele order for `*,G,<NON_REF>` is `[A, *, G, <NON_REF>]`, so PL index 4 is the
`(*,G)` cell.  `-all-sites`/`--include-non-variant-sites` are two spellings of
one `@Argument` (`GenotypeGVCFs.java:116-117,126-127`).

| # | fixture (PL of the best cell, GT field) | default | `--include-non-variant-sites` |
| --- | --- | --- | --- |
| E1 | `*,G,<NON_REF>`, PL[idx4]=0, GT `0/2` | the §1.2 row | same row (site is polymorphic) |
| E2 | `*,G,<NON_REF>`, PL[idx4]=0 **and** PL[0]=0 | **0 rows** (`G` fails the AF threshold) | `A . 127.78 … GT ./.` |
| E5 | `*,G,<NON_REF>`, PL[idx4]=0, background 40 (not 100) | **0 rows** (`log10 P(G absent)` misses `-3`) | `A . 127.78 … GT ./.` |
| E6 | E1 with background 200 | same row, QUAL `182.26`, PL `0,0,0` | – |
| C1 | `G,*,<NON_REF>` with the best cell `(G,*)`, GT `0/1` | the §1.2 row (byte-identical) | – |
| Q1 | E1 with GT field `2/0` | `… GT:AD:DP:PL 1/0:…` | – |
| Q2 | E1 with GT field `2|0` | `… GT:AD:DP:PL 1|0:…` | – |
| Q3 | E1 with GT field `1/2` (`*`,`G`) | `… GT:AD:DP:PL 0/1:…` | – |
| Q4 | E1 with GT field `./.` | **0 rows** | `A . 127.78 … GT ./.` |
| C2 | `AA` / `*,A,<NON_REF>` (deletion owns the best cell, PL idx4=0) | `AA A 82.19 … GT:AD:DP:PL 0/1:…` | – |

`Q1`/`Q2` are the decisive ones: GATK preserves the **copy order and the phase
bit of the source GT** (`2/0 → 1/0`, `2|0 → 1|0`), so the published genotype is
a projection of the source call, not any canonicalised/argmax genotype.  `Q3`
shows the same rule from the other side (`(1,2) = (*,G) → (ref,G) = 0/1`).

## 2. The GATK rule, with `file:line`

1. `GenotypeGVCFsEngine.createMinimalArgs()` forces
   `args.genotypeArgs.genotypeAssignmentMethod = GenotypeAssignmentMethod.PREFER_PLS`
   (`GenotypeGVCFsEngine.java:378`).
2. `GenotypingEngine.calculateGenotypes()` builds the genotypes of the *published*
   allele set with
   `AlleleSubsettingUtils.subsetAlleles(vc.getGenotypes(), defaultPloidy, vc.getAlleles(), outputAlleles, gpc, configuration.genotypeArgs.genotypeAssignmentMethod)`
   (`GenotypingEngine.java:190`) — i.e. the subsetting runs **after** the output
   ALT subset is decided by `calculateOutputAlleleSubset` (`:296-327`).
3. `AlleleSubsettingUtils.subsetAlleles()` projects the PL row onto the kept
   alleles with `subsettedPLIndices()` (`:426-447`), rescales it
   (`scaleLogSpaceArrayForNumericalStability`, `:92`), writes it with
   `gb.PL(newLikelihoods)` (`:131`) and then calls
   `GATKVariantContextUtils.makeGenotypeCall(...)` (`:133`).  The projected GQ is
   kept only when the **source** genotype had GQ (`:127-128`).
4. `GATKVariantContextUtils.makeGenotypeCall()`'s PREFER_PLS branch is a
   **two-way split** (`GATKVariantContextUtils.java:331-352`):
   * `if ((genotypeLikelihoods == null || !isInformative(genotypeLikelihoods)) && assignmentMethod == PREFER_PLS)`
     → `gb.alleles(bestMatchToOriginalGT(allelesToUse, originalGT.getAlleles()))`
     (`:333-338`).  **No GQ is assigned on this arm** (the `gb.log10PError(gq)` at
     `:350-352` is inside the other arm).
   * otherwise the argmax of the projected row decides the call, including the
     "uninformative hom-ref → no-call" arm
     (`maxLikelihoodIndex == 0 && gq > SUM_GL_THRESH_NOCALL`, `:351-352`).
   * `isInformative(gls)` is `MathUtils.sum(gls) < SUM_GL_THRESH_NOCALL` with
     `SUM_GL_THRESH_NOCALL = -0.1` (`:54-58`).
   * `bestMatchToOriginalGT(allelesToUse, originalGT)` maps the **source allele
     list element by element**: an allele present in `allelesToUse` keeps its
     identity, any other allele — a pruned ALT, the symbolic `*`, `<NON_REF>` —
     becomes `allelesToUse.get(0)`, i.e. the reference; a no-call stays a no-call
     (`:397-403`).
5. **The projected row that is actually stored and published is min-shifted**:
   `GenotypeBuilder.PL(double[])` → `GenotypeLikelihoods.fromLog10Likelihoods().getAsPLs()`
   → `GLsToPLs()` computes `round(-10 * (log10L[i] - max(log10L)))`, i.e. the
   stored PL row is `pl_i - min(pl_i)` (htsjdk 4.2.0, verified by disassembling
   `htsjdk/variant/variantcontext/GenotypeBuilder.PL(double[])` and
   `GenotypeLikelihoods.GLsToPLs`).  For this fixture the stored row is `0,0,0`,
   which is what both tools publish.
   Consequently, on integer PLs, `isInformative` is exactly
   `sum(pl_i - min_pl) > 1` — the test used by the fix (§3).  This is
   *shift-invariant*, so it is also correct for a row native has not normalised.
6. `AC`/`AF`/`AN` are **not** recomputed from the pre-subset call: they come from
   the `StandardAnnotation` pass that runs after the genotypes are final —
   `ChromosomeCounts.annotate()` → `VariantContextUtils.calculateChromosomeCounts`
   over the published genotypes (`ChromosomeCounts.java:43-53`, applied from
   `GenotypeGVCFsEngine.regenotypeVC()`, `GenotypeGVCFsEngine.java:187-190`).
   `MLEAC`/`MLEAF` come from `GenotypingEngine.composeCallAttributes()`
   (`:434-441`) with `calculateMLEAlleleFrequencies()` dividing by the AN of the
   published genotypes (`:452-455`); the allele-frequency posterior itself
   (`AFresult`) is still computed on the **unsubsetted** record
   (`GenotypingEngine.java:154`), which is why `MLEAC=1` survives even though
   the emitting ALT list is built later.
   All six columns of §1.2 difference therefore follow from the single GT
   difference: `AN=2, AC=1, AF=0.500, MLEAF=1/2, QD=82.26/20=4.11`,
   `ExcessHet=0.0000` for the heterozygous call, and no `GQ` because the
   fallback arm never assigns one and the source FORMAT carried none.
7. **When a genotype becomes a no-call instead** is unchanged from the previous
   round: that is the *all-ALT-collapsed* case
   (`outputAlleles.size() == 1` → `GATKVariantContextUtils.subsetToRefOnly`,
   `GenotypingEngine.java:190` plus `cleanupGenotypeAnnotations`,
   `GenotypeGVCFsEngine.java:191-194,479-491`), which the previous round already
   fixed and which this round did not touch.  Within PREFER_PLS itself a no-call
   can also arise on the *informative* arm only
   (`GATKVariantContextUtils.java:351-352`); the fix now scopes native's existing
   hom-ref/GQ-0 no-call rule to that arm (§3).

## 3. STEP 3 — the fix

### 3.1 What was wrong in native

`apply_gatk_output_allele_subset()` remapped the record onto the output ALT list
and then let `derive_gt_gq_from_pl()` assign GT/GQ from the projected PL row.
For an orphan `*` that row is constant, so native derived hom-ref with GQ 0, and
its pre-existing "uninformative hom-ref → no-call" rule turned the sample into
`./.`.  GATK's second arm was missing entirely.

A first attempt projected `record.gt` as it stood just before that assignment —
and **that was wrong**, which the new gated cases expose: the union stage already
runs `derive_gt_gq_from_pl()` (`genotype_gvcf_tool.cpp`, the aggregate and
streaming union loops), so `record.gt` at that point is the PL argmax genotype
over the union alleles (`(*,G)` for both `*,G,<NON_REF>` and `G,*,<NON_REF>`),
not the source call.  Projecting it happened to give the right row for the
promoted fixture but produced `GT=1/0` where GATK publishes `0/1` for
`G,*,<NON_REF>`, i.e. the copy order came from the PL argmax instead of the
source genotype.  The `...-star-not-first`, `...-reversed-source-gt` and
`...-phased-source-gt` cases exist to pin exactly this.

### 3.2 The change (all in `fastgatk-native/src/genotype_gvcf_tool.cpp`)

1. `Record` gains `source_gt` and `source_alleles`: the **merged source call in
   the union allele space**, captured where GATK's
   `ReferenceConfidenceVariantContextMerger` would still have it.
2. Both union loops (aggregate and stream-by-locus) snapshot
   `record.source_gt = record.gt; record.source_alleles = record.alleles;`
   immediately after `remap_record_to_allele_union()` — before
   `remove_non_ref_allele()` and before the loop's `derive_gt_gq_from_pl()`
   consume it.  Capturing *before* `<NON_REF>` removal keeps the allele list
   identical to GATK's merged VC (including `<NON_REF>` and `*`), so the later
   projection by allele names reproduces `bestMatchToOriginalGT` exactly.
3. `merge_sample_fields()` merges `source_gt` with the same sample projection as
   GT (Host-only; no FORMAT tag) and keeps the first group member's
   `source_alleles`, which every member shares by construction.
4. `genotype_record_bytes()` accounts for both new vectors.
5. New helper `gatk_prefer_pls_row_is_uninformative(pl, begin, end)`: the exact
   integer form of `isInformative` (`sum(pl_i - min_pl) <= 1`).  An entirely
   absent row returns false, i.e. the `genotypeLikelihoods == null` arm of the
   Java condition is deliberately **not** implemented (native keeps its existing
   no-call there; see §8).
6. In `apply_gatk_output_allele_subset()`, after the Number=G remap and the
   orphan-span PL normalisation:
   * mark each sample whose projected row is uninformative;
   * for those samples (skipping any with a pre-existing explicit no-call), set GT
     to the **source** call projected by allele name — surviving allele keeps its
     index, anything else (pruned ALT, `*`, `<NON_REF>`) becomes the reference,
     no-call stays no-call, and the source copy order and phase bit are kept;
   * on that arm no GQ is assigned, so a source without FORMAT/GQ leaves the
     sample's GQ missing, and if **no** sample carries a GQ the FORMAT key is
     dropped entirely (which is what makes the published row `GT:AD:DP:PL`
     instead of `GT:AD:DP:GQ:PL`);
   * the existing "hom-ref with GQ 0 → no-call" rule is now scoped to the
     informative arm only (`GATKVariantContextUtils.java:351-352` requires
     `isInformative`), which is also why the fallback can publish a call.

Not restructured: the genotyping engine, the AF/PL Kokkos kernels and their ABIs
(`remap_genotype_pl_kokkos` still requires `target_allele_count >= 2`, so the
previous round's Host-side `subsetToRefOnly` boundary is untouched), the
`max_alternate_alleles` path and every other tool.  The change is confined to one
translation unit.  Full diff: `.diag/prefer-pls-fix.diff`.

### 3.3 Measured effect

| case | GATK | native before | native after |
| --- | --- | --- | --- |
| `surviving-concrete-alt` (E1) | `A G 82.26 . AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11 GT:AD:DP:PL 0/1:0,20:20:0,0,0` | `AC=0;AF=0.00;AN=0;…;GT:AD:DP:GQ:PL ./.:0,20:20:0:0,0,0` | byte-identical to GATK |
| `…-star-not-first` (C1) | byte-identical to E1 | `./.` as above | byte-identical to GATK |
| `…-reversed-source-gt` (Q1) | `… GT:AD:DP:PL 1/0:…` | `./.` | `1/0` (matches) |
| `…-phased-source-gt` (Q2) | `… GT:AD:DP:PL 1|0:…` | `./.` | `1|0` (matches) |

## 4. STEP 2 — the strict gate

`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`:

* `surviving-concrete-alt` **promoted from REPORTED ONLY to gated** (`gated:
  True`, `expect=[GATK_G_ROW]`), with its `why` replaced by the §2 citations.
* Three cases added, all gated, each with its measured GATK row in `expect`:
  * `surviving-concrete-alt-star-not-first` (`G,*,<NON_REF>`; pruned allele not
    first) — isolates the positional projection;
  * `surviving-concrete-alt-reversed-source-gt` (`GT 2/0`) — isolates "project the
    source call" from any argmax/canonicalised genotype;
  * `surviving-concrete-alt-phased-source-gt` (`GT 2|0`) — pins phase survival.
* Only one case stays REPORTED ONLY: `surviving-deletion-alt`
  (`AA` / `*,A,<NON_REF>`), added this round because it exposes a **different,
  unfixed** root cause (§8).  The three collapse cases that were already gated
  remain gated and unchanged.
* The module docstring gained an "orphan-`*`-with-surviving-ALT cases" section
  recording the rule so the gate is self-describing.

Exit status: **1 before the fix** (`status: divergence`; literal rows in
`.diag/prefer-pls-gate-before.log`), **0 after** (`.diag/prefer-pls-gate-after.log`,
`"status": "pass"`, 8 cases — 7 gated, 1 reported only — 0 violations).
The 7 gated cases are `default`, `include-non-variant-sites`,
`stand-call-conf-50-include-non-variant-sites`, `surviving-concrete-alt`,
`surviving-concrete-alt-star-not-first`,
`surviving-concrete-alt-reversed-source-gt`,
`surviving-concrete-alt-phased-source-gt`.

The exact "before" violations were:

```
[surviving-concrete-alt] VIOLATION: row 0 is not byte-identical:
  GATK='… AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11  GT:AD:DP:PL  0/1:0,20:20:0,0,0'
  NATIVE='… AC=0;AF=0.00;AN=0;DP=20;MLEAC=1;MLEAF=0.00                  GT:AD:DP:GQ:PL  ./.:0,20:20:0:0,0,0'
[surviving-concrete-alt-star-not-first] same two rows
```

## 5. STEP 4 — stale assertions

**No registered assertion needed correcting.**  Every candidate was checked
against the new contract and each one either passes unchanged or exercises a
different fixture:

| File | Checked | Outcome |
| --- | --- | --- |
| `verify_genotype_gvcf.py:369-435` | the same locus, corrected in the previous round (`assert star_records == []` at :400 and the literal dense-mode row at :416-417) | its fixture uses PL `0,0,100,100,100,100,100,100,100,100` with GT `0/1`, whose projection onto `[A,G]` is `[0,100,100]` — **informative** (`sum(pl-min) = 200 > 1`) — so both tools prune every ALT and the REF-only collapse applies; assertion unchanged and passing |
| `verify_genotype_gvcf.py:544-579,584-625` | no-call assertions | they run `--genotype-assignment-method BEST_MATCH_TO_ORIGINAL` / `SET_TO_NO_CALL_NO_ANNOTATIONS` (and no `--gatk-compatible-annotations`), so `apply_gatk_output_allele_subset()` returns before any of this; unchanged |
| `verify_genotype_gvcf_spanning_deletion_gatk_oracle.py:96-100` | `orphan` case asserting `0/1` for the sample and `./.` for the other | different fixture (`C,*,<NON_REF>`, two samples); unchanged, still passing |
| `verify_genotype_gvcf_max_alternate_alleles_gatk_oracle.py:95` | `./.` for the sample whose ALT is dropped at `--max-alternate-alleles 6` | max-ALT path, asserted on **GATK's own row** and then compared with native; unchanged, still passing |
| `verify_genotype_gvcf_{assignment,exclude_intervals,gp_input,include_non_variant,multisample_reference_confidence,starts_in_intervals}_gatk_oracle.py` | grep for `*`-in-ALT and no-call expectations | no dependence on the projected-constant-PL path; all pass unchanged (15/15 in §6) |

Lines touched in this round: none in any registered test script — the only test
changes are the oracle promotion/additions listed in §4.

## 6. STEP 5 — gate results (mandatory order)

| Step | Command | Result |
| --- | --- | --- |
| a | `python3 fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` | **exit 0** (`status: pass`, 8 cases, 0 violations); **exit 1 before the fix** |
| b1 | `python3 fastgatk-native/scripts/verify_genotype_gvcf.py` | **exit 0** — `{"status": "pass", "output_records": 1}` |
| b2 | `ctest --test-dir fastgatk-native/build -R 'genotype-gvcf' -V` | **exit 0 — 15/15 passed** (282.66 s) |
| c | the 17-name strict-gate filter from the task | **exit 0 — 19/19 matched tests passed** (1204.23 s) |
| d | `run_regression.sh --label prefer-pls-subset` | see §6.1 |

Logs: `.diag/prefer-pls-gate-before.log`, `.diag/prefer-pls-gate-after.log`,
`.diag/prefer-pls-ctest-genotype-gvcf.log`, `.diag/prefer-pls-strict-gates.log`,
`.diag/prefer-pls-regression.log`.

### 6.1 Full double-backend regression

`fastgatk-native/scripts/run_regression.sh --label prefer-pls-subset`
(omp + serial, `FASTGATK_REQUIRE_GATK_ORACLE=1` by default).

```
| 后端   | 构建目录                               | 结果 | 通过/总数 | 耗时  |
| omp    | OpenMP (fastgatk-native/build)         | 通过 | 300/300   | 1294.55 sec |
| serial | Serial (fastgatk-native/build-serial)  | 通过 | 300/300   | 1341.82 sec |
```

`REG_EXIT=0`; evidence block `.diag/regression/20260912-010049/{evidence.md,omp.log,serial.log,omp.status,serial.status}`.
Both trees were rebuilt **after** the last source edit and before the run
(`genotype_gvcf_tool.cpp` mtime `1789143718` < `build/fastgatk-genotype-gvcf`
`1789143731` < `build-serial/fastgatk-genotype-gvcf` `1789144849`).

Caveat, stated plainly: one `file:line` citation inside a comment (line 2849,
`:367-375` → `:333-338, definition at :397-403`) was corrected **after** the
suite had finished.  The edit is comment-only and line-count neutral; both trees
were rebuilt afterwards and the resulting binaries are **byte-identical** to the
ones the suite ran on
(`md5sum -c .diag/prefer-pls-binaries-before-comment-fix.md5` → both `OK`,
`9bfff57f3fbf87d13085012033b4351b` omp /
`7cd35836b01a2d7f7e8b08796d4acb73` serial), and the oracle and
`verify_genotype_gvcf.py` were re-run on them (exit 0 both).  No statement, type
or control flow was touched after the run.

## 7. Tree state

* Change **kept in the tree** (no commit, no branch).
* Modified: `fastgatk-native/src/genotype_gvcf_tool.cpp`,
  `fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`
  (`git diff --stat`: 2 files, +308/-8; the diff is saved as
  `.diag/prefer-pls-fix.diff`).
* Not touched: `fastgatk-native/CMakeLists.txt`, root `*.md`, Mutect2 and every
  other tool; no registered test script needed an edit.

## 8. What remains unproven / open

* **A `./.` source call now implies "GATK writes no record".**  Measured (Q4):
  with `GT=./.` and the same PLs GATK emits **0 rows** while native emits the
  call for `G`.  Reason: after the fallback projects the source call it is a
  no-call, `GATKVariantContextUtils.isProperlyPolymorphic()` is false, and
  `regenotypeVC()` returns `null` (`GenotypeGVCFsEngine.java:188-190`).  Native
  has no equivalent gate here.  This round did **not** change that behaviour (it
  emitted a row with `./.` before too), it is not gated, and it is a different
  root cause; a fixture is cheap if a later round wants it.
* **`surviving-deletion-alt` (reported only, still divergent).**  With
  `AA` / `*,A,<NON_REF>` and the deletion owning the best cell, GATK prunes the
  orphan `*` (`AA A 82.19 … GT:AD:DP:PL 0/1:0,20:20:0,0,0`) but native keeps it
  (`AA *,A … GT:AD:DP:GQ:PL 1/2:0,0,20:20:99:…`).  Native's
  `orphan_spanning_deletion` flag is decided by the caller's
  `group_has_supported_spanning_deletion()`/`deletion_spans` ownership check,
  which is a different mechanism from GATK's `isVcCoveredByDeletion`
  (`GenotypingEngine.java:316`) and is out of this round's scope.  Added to the
  oracle as a non-gated case so the divergence is visible.
* **The `genotypeLikelihoods == null` arm** of PREFER_PLS is deliberately not
  implemented: a record whose projected PL row is entirely absent keeps native's
  existing no-call instead of projecting the source call.  Unmeasured (no fixture
  with a missing PL row on a pruned-allele record was built).
* **`source_has_gq` is still a per-record flag.**  The fallback clears a sample's
  GQ when the record's sources carried no GQ, mirroring
  `AlleleSubsettingUtils.java:127-128`; for a *multi-shard* input whose shards
  differ in FORMAT this is an approximation.  Unmeasured.
* **Multi-shard `source_gt` merging** is written by symmetry with the GT merge
  but has only been exercised with single-input-file fixtures; the fallback
  silently does nothing when the stored call and the GT vector disagree in shape.
  *Speculation*: this is conservative, but a multi-shard orphan-`*` oracle would
  be needed to call it verified.
* **QD/ExcessHet/PL/MLEAF** now agree for these fixtures, but they were only
  measured here; other INFO keys are still carried through untouched by the
  subsetting path (unchanged from the previous round's note).
* The previous round's remaining items are untouched and still open:
  `FILTER=LowQual` is not implemented in native GenotypeGVCFs, `QUAL=Infinity` is
  rendered `inf`/`0`, and `-all-sites` is not accepted as a short alias.
