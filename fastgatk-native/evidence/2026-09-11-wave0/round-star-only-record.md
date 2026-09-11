# Round: the spanning-deletion-only record (`covered-star-only-record`)

Scope: the REPORTED-ONLY case `covered-star-only-record` in
`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` — a locus
whose only surviving ALT is the symbolic spanning deletion `*`, which GATK drops
(`GenotypeGVCFs.java:327-328`, `GATKVariantContextUtils.isSpanningDeletionOnly`)
while native emitted it.

## 0. Verdict

* GATK's rule is **two** rules, and the task statement quoted only the second
  one. Both are measured here:
  1. `GenotypingEngine.calculateGenotypes()` returns `null` whenever the
     **output allele set is exactly `[SPAN_DEL]`** and the traversal is not
     `EMIT_ALL_ACTIVE_SITES` (`GenotypingEngine.java:173-175`). This is the rule
     that would drop the record if the `*` survived the allele subset; in the
     measured default profile the drop is over-determined, because the `*` is
     pruned by the AF threshold first and `:167-169` then returns `null` for the
     empty ALT set (§1.3).
  2. `GenotypeGVCFs.apply()` writes the regenotyped record only when
     `forceOutput || !GATKVariantContextUtils.isSpanningDeletionOnly(...)`
     (`GenotypeGVCFs.java:326-328`), a second independent gate on the same
     condition. In dense mode `forceOutput` is true for **every** locus, so
     rule 2 never fires there — density disarms both gates at once.
* **Dense mode does change the outcome**: with `--include-non-variant-sites`
  GATK writes the locus as a **REF-only no-call** row
  (`chr1 3 . A . Infinity . DP=20;MLEAC=.;MLEAF=. GT ./.`), because the same
  `calculateOutputAlleleSubset()` prunes the `*` for the ordinary
  standard-confidence reason and the empty ALT set is materialized by
  `subsetToRefOnly()` + `cleanupGenotypeAnnotations()`. It is *not* written with
  `ALT='*'` in any configuration I could measure.
* The measured native defect is **not** only "native writes a `*` record". The
  surviving predicate was `(spanning_deletion && !orphan) || (!spanning_deletion
  && plausible)` — i.e. an **owned `*` was exempt from the AF threshold**
  entirely, where GATK applies `isPlausible` to `*` exactly as to a concrete ALT
  (`GenotypingEngine.java:312-316`). A second, previously unreported divergence
  falls out of that: `covered-star-implausible-plus-concrete-alt` (an implausible
  `*` beside a plausible concrete ALT) also disagreed in the *default* mode.
* Fix applied and kept in the tree: **one predicate** in one translation unit
  (`fastgatk-native/src/genotype_gvcf_tool.cpp`, `apply_gatk_output_allele_subset`),
  `if (owned && plausible)`. No kernel, no ABI change, no `CMakeLists.txt`, no
  Mutect2, no other tool, no restructuring.
* Oracle: **exit 1 before** the fix (3 violations, literal rows in §3),
  **exit 0 after** (15 gated + 2 reported-only cases, 0 violations).
  `covered-star-only-record` is now **GATED**, with three new gated isolating
  cases and two new REPORTED-ONLY cases (§3).
* One divergence in this area is **still open and deliberately not fixed**: the
  input-records-vs-emitted-alleles structural difference flagged by
  `.diag/round-star-ownership.md` §8. My fix does *not* remove it; it *unmasks*
  it, and I reproduced it in the default output mode (§7). Fixing it needs the
  ordered emitted-deletions state, which the previous round already ruled out of
  scope; it is carried as a REPORTED-ONLY case.

## 1. STEP 1 — measured GATK truth

### 1.1 Fixture and command

The fixture is the oracle's own `STAR_ONLY_COVERED_RECORD`
(`verify_genotype_gvcf_spandel_gatk_oracle.py:228-234`); the reference is the
generated 100 bp `chr1` of `A` with hand-written `.fai`/`.dict` that the oracle
already builds (`write_reference`, `:496-505`); the plain `.vcf` is indexed with
`gatk IndexFeatureFile` for both tools:

```
chr1	2	.	AA	A,<NON_REF>	.	PASS	DP=20	GT:DP:AD:PL	0/1:20:0,20,0:100,0,100,100,100,100
chr1	3	.	A	*,<NON_REF>	.	PASS	DP=20	GT:DP:AD:PL	0/1:20:0,20,0:100,100,100,0,100,100
```

```
third_party/jdk17/bin/java -Xmx1g -jar \
  third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
  GenotypeGVCFs -R <chr1.fa> -V <fixture>.g.vcf -O out.vcf \
  --create-output-variant-index false [OPTIONS]
```

### 1.2 Literal rows per configuration (measured, pinned GATK 4.6.2.0)

| fixture | OPTIONS | GATK 4.6.2.0 data rows |
| --- | --- | --- |
| `2 AA A` + `3 A *` (**the reported case**) | *(default)* | **1 row** — `chr1 2 . AA A 92.60 . AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63 GT:AD:DP:GQ:PL 0/1:0,20:20:99:100,0,100` — **no pos-3 row** |
| same | `--include-non-variant-sites` | the row above **plus** `chr1 3 . A . Infinity . DP=20;MLEAC=.;MLEAF=. GT ./.` |
| same, star-only record `FILTER=LowQual`, `DP=7` | *(default)* | identical: only the pos-2 row |
| same, `FILTER=LowQual`, `DP=7` | `--include-non-variant-sites` | pos-2 row plus `chr1 3 . A . Infinity . DP=7;MLEAC=.;MLEAF=. GT ./.` |
| `2 AA A` + `3 A *,G` with the best cell `(*,G)` | *(default)* | **2 rows** — pos-2 row plus `chr1 3 . A *,G 82.26 . AC=1,1;AF=0.500,0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1,1;MLEAF=0.500,0.500;QD=4.11 GT:AD:DP:GQ:PL 1/2:0,0,20:20:99:100,100,100,100,0,100` (the record is **KEPT**) |
| `2 AA A` + `3 A *,G` with the best cell `(A,G)` | *(default)* | 2 rows — pos-2 row plus `chr1 3 . A G 92.63 . AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63 GT:AD:DP:GQ:PL 0/1:0,20:20:99:100,0,100` (the `*` is pruned **allele-wise**, the record survives) |
| no upstream deletion, `3 A *` | *(default)* | **0 rows** |
| no upstream deletion, `3 A *` | `--include-non-variant-sites` | `chr1 3 . A . Infinity . DP=20;MLEAC=.;MLEAF=. GT ./.` |
| `2 AA A` (*implausible*) + `3 A *,G` | *(default)* | 1 row — `chr1 3 . A G 82.26 . AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11 GT:AD:DP:PL 0/1:0,20:20:0,0,0` (the `*` is spurious; see §7) |

Native rows for the same runs are in `.diag/star-only-probe-both.log` (measured
**before** the fix) and `.diag/star-only-probe-after.log`.

Four literal rows, tab-unfolded for readability (the real bytes are separated by
tabs). Where the record is dropped, the literal evidence is its absence:

```
# default, reported case — the ONLY row GATK writes
chr1	2	.	AA	A	92.60	.	AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63	GT:AD:DP:GQ:PL	0/1:0,20:20:99:100,0,100

# native BEFORE the fix (extra record; this is the divergence)
chr1	3	.	A	*	0	.	AC=0;AF=0.00;AN=2;DP=20;ExcessHet=0.0000;MLEAC=0;MLEAF=0.00	GT:AD:DP:GQ:PL	0/0:0,20:20:0:100,100,100

# dense mode, reported case — the row GATK materializes instead
chr1	3	.	A	.	Infinity	.	DP=20;MLEAC=.;MLEAF=.	GT	./.

# native AFTER the fix, dense mode (shape matches; QUAL differs — §7)
chr1	3	.	A	.	0	.	DP=20;MLEAC=.;MLEAF=.	GT	./.
```

Note what GATK does with the record's INFO/FORMAT **when it is kept in some
configuration**: it is kept only in dense mode, and then it is kept *without the
`*`* — the `*` is pruned, `subsetToRefOnly()` supplies the genotypes
(`GenotypingEngine.java:190`) and `cleanupGenotypeAnnotations(result, true, false)`
replaces them with `./.`, dropping every FORMAT field but `GT`
(`GenotypeGVCFsEngine.java:191-194`, `:479-491`). `INFO/DP` survives from the
merged record and the Number=A `MLEAC`/`MLEAF` vectors become missing. The
record's own `FILTER=LowQual` and `AD` never reach the output (only the *locus'*
`INFO/DP=7` does). So there is no configuration in which GATK publishes a
`*`-only record at all.

### 1.3 Is the `*` plausible or spurious in the reported case?

Decisive measurement: under `--include-non-variant-sites` GATK writes
`ALT='.'`, **not** `ALT='*'`. Had the `*` survived the output-allele subset,
`GenotypingEngine.java:173-175` would be skipped (`emitAllActiveSites()` is true)
and the rebuilt VariantContext would carry `ALT='*'`. It does not, so
`calculateOutputAlleleSubset()` pruned the `*` — and since the `*` *is* covered
by the emitted upstream deletion (proved by the `covered-star-plus-concrete-alt`
control, where the same covering deletion keeps a `*`), the pruning is the
`passesThreshold()` arm, not the `isSpuriousSpanningDeletion` arm
(`GenotypingEngine.java:312-316`). The default-mode drop is therefore
over-determined: it comes from `:167-169` *and* would come from `:173-175`. Both
are cited in the gate's `why` text.

## 2. The GATK rule, with `file:line`

1. `GenotypingEngine.calculateOutputAlleleSubset()` (`:296-328`) decides the ALT
   subset. For every non-reference allele:
   ```java
   final boolean isPlausible = afCalculationResult.passesThreshold(allele, configuration.genotypeArgs.standardConfidenceForCalling);   // :312
   final boolean isSpuriousSpanningDeletion = GATKVCFConstants.isSpanningDeletion(allele) && !isVcCoveredByDeletion(vc);              // :314
   final boolean toOutput = (isPlausible || forceKeepAllele(allele) || isNonRefWhichIsLoneAltAllele || forcedAlleles.contains(allele))
                            && !isSpuriousSpanningDeletion;                                                                             // :316
   ```
   There is **no exemption for `*`**: a spanning deletion must be *both*
   plausible and owned. (`MinimalGenotypingEngine.forceKeepAllele()` returns
   `configuration.annotateAllSitesWithPLs`,
   `MinimalGenotypingEngine.java:62-64` — false for GenotypeGVCFs.)
2. `GenotypingEngine.calculateGenotypes()` drops the record outright when the
   surviving ALT list is exactly the spanning deletion (`:173-175`):
   ```java
   // return a null call if we aren't forcing site emission and the only alt allele is a spanning deletion
   if (! emitAllActiveSites() && outputAlternativeAlleles.alleles.size() == 1 && Allele.SPAN_DEL.equals(outputAlternativeAlleles.alleles.get(0))) {
       return null;
   }
   ```
   `emitAllActiveSites()` is `configuration.outputMode == OutputMode.EMIT_ALL_ACTIVE_SITES`
   (`:421-423`), and `GenotypeGVCFsEngine.createMinimalArgs(forceOutput)` sets
   `args.outputMode = forceOutput ? EMIT_ALL_ACTIVE_SITES : EMIT_VARIANTS_ONLY`
   (`GenotypeGVCFsEngine.java:381-383`), with `forceOutput = includeNonVariants ||
   inForceOutputIntervals` (`GenotypeGVCFs.java:322-323`).
3. `GenotypeGVCFs.apply()` is the second, independent gate
   (`GenotypeGVCFs.java:320-330`):
   ```java
   final boolean forceOutput = includeNonVariants || inForceOutputIntervals;                    // :322
   final VariantContext regenotypedVC = gvcfEngine.callRegion(..., forceOutput);                // :323
   if (regenotypedVC != null) {
       if ((forceOutput || !GATKVariantContextUtils.isSpanningDeletionOnly(regenotypedVC))) {   // :326-327
           vcfWriter.add(regenotypedVC);
       }
   }
   ```
   `isSpanningDeletionOnly()` is `getAlternateAlleles().size() == 1 &&
   isSpanningDeletion(getAlternateAllele(0))`
   (`GATKVariantContextUtils.java:2089-2091`).
4. **Does dense mode change it?** Yes, but not by relaxing rule 2: rule 2 needs
   `!emitAllActiveSites()`, and rule 3 needs `!forceOutput`. `includeNonVariants`
   makes `forceOutput` true for every locus, so *both* gates are disarmed at once.
   What is then published is decided one level up: the `*` is pruned by rule 1
   (it is implausible here), leaving an empty ALT set, and
   `calculateGenotypes()` builds the REF-only call through
   `GATKVariantContextUtils.subsetToRefOnly(vc, defaultPloidy)`
   (`GenotypingEngine.java:190`), which `regenotypeVC`'s monomorphic branch
   (`GenotypeGVCFsEngine.java:360-364`) turns into the no-call row. Nothing about
   density re-enables a `*`-only record.
5. The record-level rule is therefore *not* "drop every record carrying a `*`":
   `covered-star-plus-concrete-alt` (two surviving ALTs) is kept and genotyped
   `1/2`, and `covered-star-implausible-plus-concrete-alt` keeps the record while
   dropping only the `*` allele.

### 2.1 Native before the fix, and why it was wrong

`fastgatk-native/src/genotype_gvcf_tool.cpp:2911-2936` (before the fix; the same
region is `:2911-2937` after):

```cpp
// GATK treats '*' as a structural non-variant only when a concrete
// deletion owns the span.  An orphan '*' is removed regardless of
// its AF posterior; ordinary ALTs remain subject to the posterior
// confidence threshold.
if ((spanning_deletion && !record.orphan_spanning_deletion) ||
    (!spanning_deletion && plausible))
    output_alleles.push_back(allele);
else
    ++pruned;
```

An **owned** `*` was admitted without ever consulting `plausible`, which is the
`GATKTF`-inverse of `GenotypingEngine.java:316`: GATK's `!isSpuriousSpanningDeletion`
is an *additional* condition on top of `isPlausible`, not a replacement for it.
The comment's claim ("the AF threshold applies to ordinary concrete ALTs") is the
precise defect.

### 2.2 Interaction with the previous round's "input records, not emitted alleles"

The previous round found (`round-star-ownership.md` §8, *speculation* there) that
native's `deletion_spans` is filled from every **input** record with a concrete
deletion allele (`genotype_gvcf_tool.cpp:5173` + `:5304` in the aggregate path,
`:6199-6203` in the stream-by-locus path, predicate `record_has_concrete_deletion`
at `:1208`), while GATK's `upstreamDeletionsLoc` only ever holds deletions it
actually **emitted** (`GenotypingEngine.java:178-179` → `:343-357`).

Interaction, measured (§7): before this fix, a record whose only ALT was an owned
`*` was being written by the star-only bug anyway, so the input-vs-emitted
difference was **masked** in the default profile. It is not masked any more once
the `*` is required to be plausible *and* the record is dropped when `*` alone
survives — but the difference itself is untouched: for a case where GATK drops
the upstream deletion record and the downstream `*` is *plausible*, native still
calls it owned. That case is now
`unemitted-upstream-deletion-star-plus-concrete-alt` (REPORTED ONLY). My fix does
**not** address it and does not need to; it needs the ordered emitted-deletions
restructuring the previous round ruled out of scope. I did not attempt it.

## 3. STEP 2 — the strict gate (written before the fix)

`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`:

* `covered-star-only-record` **promoted from REPORTED ONLY to GATED**
  (`gated: True`, `expect=[GATK_DEL_UPSTREAM_LOCUS_ROW]`), `why` replaced by the
  §1/§2 citations.
* Three new **GATED** isolating cases, each with measured GATK rows in `expect`:
  * `covered-star-only-record-non-pass` — the star-only record is
    `FILTER=LowQual`, `DP=7`; GATK's output is bit-identical to the promoted
    case, proving the rule is a property of the output allele set, not of the
    record's FILTER/INFO/FORMAT content.
  * `covered-star-plus-concrete-alt` — **positive control**: `*,G` with the
    `(*,G)` cell best, so the `*` is plausible *and* covered; GATK keeps the
    record and calls `1/2`. Without this, "drop every record carrying a `*`"
    would look correct.
  * `covered-star-implausible-plus-concrete-alt` — the same locus with the best
    cell at `(A,G)`: the `*` is covered but implausible, GATK prunes the
    **allele** and publishes `ALT='G'` with `0/1`. This is the case that forces
    the `&& plausible` half of the fix; it was **not** covered by the reported
    divergence and it failed before the fix.
* Two new **REPORTED-ONLY** cases:
  * `covered-star-only-record-dense` — the dense-mode statement of the same
    locus, with the measured GATK REF-only row (and the residual QUAL
    divergence, §7) recorded in `expect`;
  * `unemitted-upstream-deletion-star-plus-concrete-alt` — the residual
    input-vs-emitted divergence of §2.2/§7.
* The module docstring gained a "The spanning-deletion-only record" section and a
  paragraph on the reported-only residual divergence, so the gate is
  self-describing. No existing case's expectation changed.

### 3.1 Exit status and literal rows — BEFORE the fix

**exit 1**, `"status": "divergence"`, **3 violations**
(`.diag/star-only-gate-before.log`):

```
[covered-star-only-record] record count differs: GATK=1 native=2
    NATIVE='chr1\t3\t.\tA\t*\t0\t.\tAC=0;AF=0.00;AN=2;DP=20;ExcessHet=0.0000;MLEAC=0;MLEAF=0.00\tGT:AD:DP:GQ:PL\t0/0:0,20:20:0:100,100,100'
[covered-star-only-record-non-pass] record count differs: GATK=1 native=2
    NATIVE='chr1\t3\t.\tA\t*\t0\tRGQ\tAC=0;AF=0.00;AN=2;DP=7;ExcessHet=0.0000;MLEAC=0;MLEAF=0.00\tGT:AD:DP:GQ:PL\t0/0:0,7:7:0:100,100,100'
[covered-star-implausible-plus-concrete-alt] row 1 is not byte-identical:
  GATK  ='chr1\t3\t.\tA\tG\t92.63\t.\tAC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63\tGT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100'
  NATIVE='chr1\t3\t.\tA\t*,G\t92.63\t.\tAC=0,1;AF=0.00,0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=0,1;MLEAF=0.00,0.500;QD=4.63\tGT:AD:DP:GQ:PL\t0/2:0,0,20:20:99:100,100,100,0,100,100'
```

The `covered-star-plus-concrete-alt` positive control and all eleven previously
gated cases were green **before** the fix, so the three violations isolate exactly
the allele-level qualification of `*`.

**AFTER the fix: exit 0**, `"status": "pass"`, 0 violations
(`.diag/star-only-gate-after.log`). Case census — **15 gated, 2 reported only**:

| case | gated | GATK rows | native rows |
| --- | --- | --- | --- |
| `default` … `surviving-deletion-alt-only-alt-implausible` (11 pre-existing) | yes | — | byte-identical |
| **`covered-star-only-record`** | **yes** | 1 | 1 |
| **`covered-star-only-record-non-pass`** | **yes** | 1 | 1 |
| **`covered-star-plus-concrete-alt`** | **yes** | 2 | 2 |
| **`covered-star-implausible-plus-concrete-alt`** | **yes** | 2 | 2 |
| `covered-star-only-record-dense` | no | 2 | 2 (QUAL differs) |
| `unemitted-upstream-deletion-star-plus-concrete-alt` | no | 1 | 1 (ALT differs) |

## 4. STEP 3 — the fix

One predicate in `fastgatk-native/src/genotype_gvcf_tool.cpp`
(`apply_gatk_output_allele_subset`), plus its comment; full patch
`.diag/star-only-fix.diff`:

```diff
-        // A spanning deletion is retained as a structural non-variant allele
-        // when it is present in the merged context.  The upstream deletion
-        // ownership check happens in the caller's reference-confidence
-        // merge, while the AF threshold applies to ordinary concrete ALTs.
         const bool spanning_deletion = allele == "*";
         const auto absent = record.cohort_log10_p_allele_absent[index];
         const bool plausible = std::isfinite(absent) && absent + 1.0e-10 < threshold;
-        // GATK treats '*' as a structural non-variant only when a concrete
-        // deletion owns the span.  An orphan '*' is removed regardless of
-        // its AF posterior; ordinary ALTs remain subject to the posterior
-        // confidence threshold.
-        if ((spanning_deletion && !record.orphan_spanning_deletion) ||
-            (!spanning_deletion && plausible))
+        // ... GenotypingEngine.java:312-316 quoted in full ...
+        const bool owned = !spanning_deletion || !record.orphan_spanning_deletion;
+        if (owned && plausible)
             output_alleles.push_back(allele);
         else
             ++pruned;
```

Why this is sufficient for both gates without new code: once the `*` is pruned,
`output_alleles.size() == 1` and the existing branch at
`genotype_gvcf_tool.cpp:2941-2952` already implements the whole of §2's rule 2/4 —
default profile destroys the record (nothing is written), dense profile calls
`materialize_gatk_monomorphic_ref_call()`, which is the native counterpart of
`subsetToRefOnly()` + `cleanupGenotypeAnnotations(true)` that the earlier
`-all-sites` round implemented.

Why `owned` and not simply `plausible`: `record.orphan_spanning_deletion` is set
per **record** — for every record in a group containing an unsupported `*`
(`:5582-5583`, `:6242-6244` after this edit's line shift) — so folding it into
the concrete-ALT arm would also drop a legitimate concrete ALT that shares the
record with an orphan `*`. That would regress `surviving-deletion-alt` (fixture
`AA *,A,<NON_REF>`, where GATK keeps `A`). The predicate therefore applies
`owned` only to the spanning deletion.

Reachability is unchanged: the whole function short-circuits on
`!options.gatk_annotation_compatibility` (`:2900`), so only the
`--gatk-compatible-annotations` profile is affected; the native diagnostic
profile — which `verify_genotype_gvcf.py` also exercises — is structurally
untouched.

Not restructured: the genotyping engine, the AF/PL Kokkos kernels and their ABIs,
the union/merge stages, `max_alternate_alleles`, dense-mode materialization, the
deletion-span construction, and every other tool. **In particular I did not touch
native's input-derived `deletion_spans`** (§2.2) — that is the restructuring the
task said to stop and report instead of attempting, and it is not needed for this
case.

## 5. STEP 4 — stale assertions

**No registered test pinned the old behaviour, so no test line needed
correcting.** The only test-script change in this round is the oracle itself
(§3). Verified rather than assumed:

| File / check | Why it could have been stale | Outcome |
| --- | --- | --- |
| `verify_genotype_gvcf.py:380-400` (star fixture `A *,G,<NON_REF>`, `assert star_records == []` at `:400`) | it is the only registered `*` fixture without a deletion allele | no `*` ownership state is involved at all (no record carries a concrete deletion, so no span exists either way); unchanged and passing (`B1_EXIT=0`, `{"status": "pass", "output_records": 1}`) |
| `verify_genotype_gvcf_spanning_deletion_gatk_oracle.py:47-59` (`T C,*,<NON_REF>`, the one registered fixture sitting on the ownership boundary, `"C,*" in ALT` asserted) | with `*` now also needing `passesThreshold()`, an over-eager `&& plausible` could have pruned a `*` GATK keeps | unchanged and passing — ctest `fastgatk-genotype-gvcf-spanning-deletion-gatk-oracle`, 16.20 s; GATK and native rows byte-identical |
| `verify_genotype_gvcf_assignment_gatk_oracle.py:144-145` (`assert output_alleles_pruned == 1`, `output_allele_pruning_calls > 0`) | the fix changes what gets counted by the pruning telemetry | `grep` on its fixture shows ALT `C,G,<NON_REF>` — **no `*`** — so the counted allele is a concrete one and the predicate is unchanged for it; passing (`fastgatk-genotype-gvcf-assignment-gatk-oracle`, 12.34 s) |
| the other `verify_genotype_gvcf_*_gatk_oracle.py` siblings (`multisample`, `multiallelic`, `include-non-variant`, `max-alternate-alleles`, `gp-input`, `exclude-intervals`, `starts-in-intervals`, `legacy-qual`, `multisample-reference-confidence`, `contract`, `gatk-oracle`) | any of them could carry a `*` ALT fixture | `grep` over every `verify_*genotype_gvcf*.py` for `*` in an ALT column: the only hits are the spandel/spanning-deletion scripts above and Python `*args`; all siblings pass in §6 step (b2) |
| `benchmark_genotype_gvcf.py` | reads telemetry | only reports `orphan_spanning_deletion_loci`; no asserted value |
| C++ / ctest sources | a unit test could pin `apply_gatk_output_allele_subset` | `grep -rn "output_alleles_pruned\|orphan_spanning_deletion\|spanning_deletion_supported"` over `*.cpp/*.hpp/*.py/*.sh/*.txt/*.cmake` → only `genotype_gvcf_tool.cpp` and the two scripts above; no unit test pins the predicate |

## 6. STEP 5 — gate results (mandatory order)

| Step | Command | Result |
| --- | --- | --- |
| a | `python3 fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` | **exit 0** (15 gated, 0 violations); **exit 1 before the fix** (3 violations, §3.1) |
| b1 | `python3 fastgatk-native/scripts/verify_genotype_gvcf.py` | `B1_EXIT=0` — `{"status": "pass", "output_records": 1}` |
| b2 | `ctest --test-dir fastgatk-native/build -R 'genotype-gvcf' -V` | `B2_EXIT=0` — **15/15 passed**, incl. `fastgatk-genotype-gvcf-spandel-gatk-oracle` (149.94 s) and `fastgatk-genotype-gvcf-spanning-deletion-gatk-oracle` (16.20 s) |
| c | the 17-name strict-gate filter from the task | see §6.1 |
| d | `fastgatk-native/scripts/run_regression.sh --label star-only-record` | see §6.1 |

Binaries (rebuilt in both trees **after** the last source edit and before every
gate above; the oracle script was re-frozen before the re-runs, see §6.1):
`md5sum .diag/star-only-binaries.md5` →
`265cdbfa5c22dc35e4f2ea81b3fa2bca` (omp),
`f0239c2a5fbf8e196e3e2cd8e23d6ce7` (serial).

Logs: `.diag/star-only-gate-before.log`, `.diag/star-only-gate-after.log`,
`.diag/star-only-verify-gvcf.log`, `.diag/star-only-ctest-genotype-gvcf.log`,
`.diag/star-only-strict-gates.log`, `.diag/star-only-probe-both.log`,
`.diag/star-only-probe-after.log`, `.diag/star-only-regression.log`.

### 6.1 Strict gates and the double-backend regression

* **(c)** the 17-name filter from the task, run **after** the final oracle
  revision was frozen: `C_EXIT=0`, **19/19 matched tests passed**
  (`.diag/star-only-strict-gates.log`, total 1185.16 s). The filter matched
  `fastgatk-hc-ploidy-window-invariance-gatk-oracle` and a second `mutect2-recheck`
  test in addition to the named 17. The `spandel-gatk-oracle` entry
  (`fastgatk-genotype-gvcf-spandel-gatk-oracle`, 162.17 s) is the promoted gate
  itself and ran the final revision.
* **(b2) re-run on that same final revision** (the oracle script gained two
  REPORTED-ONLY cases after the first (b2) run, so the whole family was re-run to
  keep every gate on one revision): `B2_EXIT=0`, **15/15 passed**
  (`.diag/star-only-ctest-genotype-gvcf.log`, total 318.20 s, including the
  spandel oracle and the `verify_genotype_gvcf_spanning_deletion_gatk_oracle.py`
  sibling).
* **(d)** `fastgatk-native/scripts/run_regression.sh --label star-only-record`
  (omp + serial, `FASTGATK_REQUIRE_GATK_ORACLE=1` by default, ctest parallelism
  8), `REG_EXIT=0`:

  ```
  | 后端   | 构建目录                               | 结果 | 通过/总数 | 耗时  |
  | omp    | OpenMP (fastgatk-native/build)         | 通过 | 300/300   | 1314.78 sec |
  | serial | Serial (fastgatk-native/build-serial)  | 通过 | 300/300   | 1294.90 sec |
  ```

  Evidence block `.diag/regression/20260912-041129/{omp.log,serial.log}`; the
  runner's summary names the same two trees and reports the staleness check as
  clean. Suite log `.diag/star-only-regression.log`.

**Nothing was edited after the suite run.** `md5sum -c
.diag/star-only-binaries.md5` re-confirms
`265cdbfa5c22dc35e4f2ea81b3fa2bca` (omp) and `f0239c2a5fbf8e196e3e2cd8e23d6ce7`
(serial) after the suite. Timestamps, all **before** the suite window
04:11:29 → 04:33:24 (`.diag/regression/20260912-041129/omp.log` mtime):
`genotype_gvcf_tool.cpp` 03:24:50 < omp binary 03:25:12 < serial binary 03:25:33,
and the oracle `verify_genotype_gvcf_spandel_gatk_oracle.py` 03:48:14. The only
file touched afterwards is this report, which is not a test input.

## 7. What remains unproven / open

* **The input-records-vs-emitted-alleles divergence is real and now unmasked**
  (measured this round, default output mode, not fixed, not gated):

  ```
  fixture  2 AA A,<NON_REF>   GT:DP:AD:PL 0/0:20:20,0,0:0,100,100,100,100,100   # deletion implausible
           3 A *,G,<NON_REF> GT:DP:AD:PL 0/2:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100
  GATK   chr1 3 . A G   82.26 . AC=1;AF=0.500;AN=2;DP=20;...    GT:AD:DP:PL 0/1:0,20:20:0,0,0
  NATIVE chr1 3 . A *,G 82.26 . AC=1,1;AF=0.500,0.500;AN=2;...  GT:AD:DP:GQ:PL 1/2:0,0,20:20:99:100,100,100,100,0,100
  ```
  GATK drops the upstream record (`GenotypingEngine.java:167-169`), so
  `recordDeletions()` never sees it (`:178-179`) and the downstream `*` is
  spurious (`:314`) even though it is plausible; native's span list comes from
  the **input** record and still owns it. Measured **before** as well as after
  my fix (the old predicate kept any non-orphan `*`; the new one keeps a
  plausible non-orphan `*`), so this is **pre-existing**, not a regression.
  Identity of the reported-only case
  `unemitted-upstream-deletion-star-plus-concrete-alt`. Fixing it needs the
  ordered per-locus emitted-deletions state (a real restructuring) — the same
  conclusion as `round-star-ownership.md` §8, now with a default-mode reproducer
  instead of a masked one.
* **Dense mode: native's QUAL for the materialized monomorphic row.** For the
  promoted fixture native writes `chr1 3 . A . 0 . …` where GATK writes
  `Infinity`. **Pre-existing and unrelated to the allele rule**, proven by the
  pre-fix probe: `uncovered-star-only-dense` (no star involved) already showed
  GATK `Infinity` vs native `0`, and `implausible-upstream-plausible-star-dense`
  showed GATK `Infinity` vs native `159.55` at the same time. The ALT/GT/FORMAT
  shape now matches. *Speculation*: it is the
  `log10ProbOnlyRefAlleleExists` vs `log10ProbVariantPresent` branch of
  `GenotypingEngine.java:158-163` for a site whose ALT set is empty at
  `--include-non-variant-sites`, i.e. the same family as the `127.78` row that
  does match; not investigated further, and not gated.
* **Multi-sample and multi-shard** forms of the star-only rule are unmeasured;
  every case here is single-sample, single-input. *Speculation*: the rule is
  per-locus, so it should hold, but that is untested.
* **`--force-output-intervals`** moves `forceOutput` without
  `--include-non-variant-sites` (`GenotypeGVCFs.java:322`), so rule 2 is
  disarmed while rule 1 is not (`forceOutputGenotypingEngine` is selected by
  `outputNonVariants`, `GenotypeGVCFsEngine.java:247-248`). Untested; native has
  no counterpart flag in this path.
* The pre-existing, still-open items from earlier rounds are untouched:
  `FILTER=LowQual` is not implemented in native GenotypeGVCFs, `QUAL=Infinity` is
  rendered `inf`/`0`, `-all-sites` is not accepted as a short alias, the
  `genotypeLikelihoods == null` arm of PREFER_PLS is not implemented, and a
  `./.` source call may still produce a record where GATK writes none.

## 8. Tree state

* Change **kept in the tree** (no commit, no branch, no branch switch).
* Modified — exactly two files, `git status --short`:

  ```
   M fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py
   M fastgatk-native/src/genotype_gvcf_tool.cpp
  ```

  `git diff --stat` = **2 files changed, 221 insertions(+), 21 deletions(-)**
  (oracle +213/-… , tool +29/-…). Scratch lives under `.diag/` only, which is
  git-ignored (`.gitignore:42`), so it does not appear in `git status`:
  `round-star-only-record.md` (this file), `star_only_probe.py`,
  `star-only-fix.diff`, `star-only-binaries.md5`,
  `star-only-{gate-before,gate-after,verify-gvcf,ctest-genotype-gvcf,strict-gates,probe-both,probe-after,regression}.log`,
  `regression/20260912-041129/`.
* Not touched: `fastgatk-native/CMakeLists.txt`, every root `*.md`, Mutect2,
  every other tool, and every registered test script other than the oracle.
* Both build trees were rebuilt after the last **source** edit and before every
  gate in §6; the oracle's last edit (a docstring section plus two
  REPORTED-ONLY cases) preceded steps (a), (b2), (c) and (d) as re-run in §6.1.
* Binaries are byte-identical to the ones the suite ran:
  `md5sum -c .diag/star-only-binaries.md5` → both `OK`
  (`265cdbfa5c22dc35e4f2ea81b3fa2bca` omp,
  `f0239c2a5fbf8e196e3e2cd8e23d6ce7` serial).
