# Round: spanning-deletion ownership — the locus' own deletion must not own its own `*`

Scope: the REPORTED-ONLY case `surviving-deletion-alt` in
`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`
(`AA` / `*,A,<NON_REF>` with the concrete deletion owning the best genotype).

## 0. Verdict

* The hand-off hypothesis — *"native's `orphan_spanning_deletion` uses its own
  ownership check rather than GATK's `isVcCoveredByDeletion`"* — is
  **CONFIRMED in substance, REFUTED in wording**. Native does implement a private
  predicate (`spanning_deletion_supported_at`), but the predicate is *the same
  test* as GATK's apart from **one boundary comparison**: native uses an
  inclusive deletion start (`span.begin <= record.pos`) where
  `GenotypingEngine.isVcCoveredByDeletion()` requires a **strictly** earlier
  start (`loc.getStart() < vc.getStart()`, `GenotypingEngine.java:369`). It is
  not a different mechanism, not a different interval convention (the end
  semantics agree exactly), and not a record-ordering bug at the call site.
* Measured cause: native's ownership set is a pre-computed list of every input
  record carrying a concrete deletion allele, so the locus' **own** deletion
  allele is in that set, and with an inclusive start it covers its own position
  (`span.begin == record.pos`). GATK cannot reproduce that: it collects
  deletions in `recordDeletions()` **after** the output-allele subset has already
  been computed (`GenotypingEngine.java:155` then `:178-179`), and its `<` test
  would exclude an equal start anyway.
* Fix applied and kept in the tree: **one line** (`<=` → `<`) plus a 9-line
  comment, in one translation unit
  (`fastgatk-native/src/genotype_gvcf_tool.cpp`, `spanning_deletion_supported_at`).
  No kernel, no `CMakeLists.txt`, no Mutect2, no other tool; the oracle is the
  only test script touched, and only because step 2 required promoting the case.
* Oracle: **exit 1 before** the fix (3 violations, literal rows in §4),
  **exit 0 after** (11 gated + 1 reported-only case, 0 violations).
* `surviving-deletion-alt` is now **GATED**, together with three new isolating
  gated cases and one new REPORTED-ONLY case (§4).

## 1. STEP 1 — measured GATK truth (literal rows + command)

### 1.1 Fixture and command

Byte-for-byte the fixture the oracle already carried
(`G_DELETION_ALT_RECORD`, now `verify_genotype_gvcf_spandel_gatk_oracle.py:161-165`):

```
chr1	2	.	AA	*,A,<NON_REF>	.	PASS	DP=20	GT:DP:AD:PL	0/2:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100
```

Reference: a generated 100 bp `chr1` of `A` with hand-written `.fai`/`.dict`;
input is a plain `.vcf` indexed with `gatk IndexFeatureFile`. Allele order is
`[AA, *, A, <NON_REF>]`, so PL index 4 is the `(*,*)` cell (value 0, the best
cell) and the source call `0/2` is `(AA, A)`. Command (identical on both sides
apart from `--gatk-compatible-annotations` for native):

```
third_party/jdk17/bin/java -Xmx1g -jar \
  third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
  GenotypeGVCFs -R <chr1.fa> -V <fixture>.g.vcf -O out.vcf \
  --create-output-variant-index false
```

### 1.2 Literal rows (pinned GATK 4.6.2.0, measured)

| # | fixture | GATK 4.6.2.0 row(s) | native BEFORE | native AFTER |
| --- | --- | --- | --- | --- |
| R | `2 AA *,A,<NON_REF>` (the reported case) | `chr1 2 . AA A 82.19 . AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11 GT:AD:DP:PL 0/1:0,20:20:0,0,0` | `chr1 2 . AA *,A 82.19 . AC=1,1;AF=0.500,0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1,1;MLEAF=0.500,0.500;QD=4.11 GT:AD:DP:GQ:PL 1/2:0,0,20:20:99:100,100,100,100,0,100` | byte-identical to GATK |
| N1 | `2 AA A,<NON_REF>` (spans 2-3) + `4 AA *,A,<NON_REF>` | `chr1 2 . AA A 92.60 . AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63 GT:AD:DP:GQ:PL 0/1:0,20:20:99:100,0,100` and `chr1 4 . AA A 82.19 . AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11 GT:AD:DP:PL 0/1:0,20:20:0,0,0` | row 1 identical; row 2 = `AA *,A … GT:AD:DP:GQ:PL 1/2:…` | byte-identical to GATK |
| N2 | `2 AAA A,<NON_REF>` (spans 2-4) + `4 AA *,A,<NON_REF>` | `chr1 2 . AAA A 92.60 . … GT:AD:DP:GQ:PL 0/1:0,20:20:99:100,0,100` and `chr1 4 . AA *,A 82.19 . AC=1,1;AF=0.500,0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1,1;MLEAF=0.500,0.500;QD=4.11 GT:AD:DP:GQ:PL 1/2:0,0,20:20:99:100,100,100,100,0,100` | byte-identical (already correct) | byte-identical to GATK |
| S | `2 AA *,A,<NON_REF>` with the only other ALT implausible (`PL 0,100,100,100,100,100,100,100,100,100`, `GT 0/1`) | **0 rows** | `chr1 2 . AA * 0 . AC=0;AF=0.00;AN=2;DP=20;… GT:AD:DP:GQ:PL 0/0:10,0:20:99:0,100,100` | **0 rows** (matches GATK) |
| C | `2 AA A,<NON_REF>` + `3 A *,<NON_REF>` (star **covered**, star-only record) | `chr1 2 . AA A 92.60 . … ` only — **no pos-3 row** | pos-2 row identical + `chr1 3 . A * 0 . AC=0;AF=0.00;AN=2;DP=20;… GT:AD:DP:GQ:PL 0/0:0,20:20:0:100,100,100` | unchanged: still writes the pos-3 row (reported only, §8) |
| D | `2 AA A,<NON_REF>` with the deletion implausible (`PL 0,100,100,100,100,100`) + `3 A *,<NON_REF>` | **0 rows** | `chr1 3 . A * 0 . …` | unchanged: still writes the pos-3 row (reported only, §8) |
| O | `2 AAA A,<NON_REF>` + `3 AA A,<NON_REF>` + `4 A *,<NON_REF>` | pos 2 and pos 3 only | pos 2, pos 3 identical + `chr1 4 . A * 0 . …` | unchanged (same star-only-record rule) |
| X | `2 AAA AA,<NON_REF>` (spans 2-3) + `4 A *,<NON_REF>` | `chr1 2 . AA A 92.60 …` only (reverse-trimmed) | `chr1 2 . AAA AA 92.60 …` (not trimmed — unrelated pre-existing divergence, §8) + pos-4 `*` row | unchanged |

Rows R, N1 and S are the divergence; N2 is the positive control that shows the
fix did **not** over-prune. The task statement's premise for R is reproduced
exactly (ALT `A`, `GT 0/1`, `AC=1;AF=0.500;AN=2`, `GT:AD:DP:PL`, `PL 0,0,0`); no
part of it needed correcting.

Evidence: `.diag/star-ownership-probe-default.log` (before the fix, all rows) and
`.diag/star-ownership-probe-after.log` (after), produced by
`.diag/star_ownership_probe.py` (scratch, `tempfile.TemporaryDirectory`).

## 2. The GATK rule, with `file:line`

1. `GenotypingEngine.calculateGenotypes()` decides the output ALT subset at
   `:155` (`calculateOutputAlleleSubset(AFresult, vc, forcedAlleles)`), computes
   the emitted allele list at `:178` (`outputAlleles(vc.getReference())`) and
   **only then** calls `recordDeletions(vc, outputAlleles)` at `:179`.
2. `calculateOutputAlleleSubset()` (`:296-328`) marks an ALT spurious when
   `isSpuriousSpanningDeletion = GATKVCFConstants.isSpanningDeletion(allele) &&
   !isVcCoveredByDeletion(vc)` (`:314`) and drops it from the output
   (`:316`); `siteIsMonomorphic` follows (`:318`).
3. `isVcCoveredByDeletion()` (`:365-371`):
   ```java
   return !upstreamDeletionsLoc.isEmpty() && upstreamDeletionsLoc.stream()
           .anyMatch(loc -> loc.getContig().equals(vc.getContig())
                         && loc.getStart() < vc.getStart()
                         && vc.getStart() <= loc.getEnd());
   ```
   Two things matter: the start test is **strict** (`<`), and the set only
   contains deletions recorded by *earlier* `recordDeletions()` calls. The state
   itself is a per-engine `PriorityQueue<Locatable> upstreamDeletionsLoc`
   (`:52`), i.e. **accumulated across the whole ordered walk**, lazily culled at
   the head of the queue against the current locus (`:344-346`); the only reset
   helper, `clearUpstreamDeletionsLoc()` (`:330-332`), is never called from
   production code (grep over `gatk-source/src/main/java/`).
4. `recordDeletions()` (`:343-357`) adds, for each **emitted** ALT,
   `new SimpleInterval(contig, vc.getStart(), vc.getStart() + deletionSize)` with
   `deletionSize = vc.getReference().length() - allele.length()` (`:349-355`) —
   an **inclusive** end, and only from a locus that was actually emitted (the
   `null` returns at `:167-169` and `:172-175` sit *before* `:179`).
5. Therefore: **a deletion allele of the very record being genotyped can never
   own that record's own `*`.** Its interval is not in the set yet, and its start
   equals the locus start, which `:369` rejects.
6. Consequences for fixture R, all measured: `*` is spurious → dropped; the
   surviving `A` is published; the genotype is the source call projected onto the
   kept alleles by `GATKVariantContextUtils.subsetAlleles()` →
   `makeGenotypeCall()`'s PREFER_PLS fallback (`0/2 → 0/1`, copy order and phase
   preserved, no `GQ` assigned, PL row min-shifted to `0,0,0`) — the mechanism the
   previous round already implemented and which needed no change here; `AC/AF/AN`
   come from the post-subset `StandardAnnotation` pass. This is why the published
   row is `AA A 82.19 … GT:AD:DP:PL 0/1:0,20:20:0,0,0`.
7. Native's counterpart (`fastgatk-native/src/genotype_gvcf_tool.cpp:1222-1230`
   before the fix):
   ```cpp
   bool spanning_deletion_supported_at(const Record& record,
                                       const std::vector<VariantSpan>& deletion_spans) {
       if (record.rid < 0) return false;
       for (const auto& span : deletion_spans) {
           if (span.rid == record.rid && span.begin <= record.pos && record.pos < span.end)
               return true;
       }
       return false;
   }
   ```
   `deletion_spans` is filled by a full pre-pass over the inputs from every
   record for which `record_has_concrete_deletion()` holds
   (`:6190-6194` aggregate/stream probe; `:5164` declaration, `:5295` fill),
   `record.pos` is 0-based and `record_span_end()` = `pos + strlen(REF)`
   (`:1188-1200`), i.e. a 0-based **exclusive** end that is numerically identical
   to GATK's 1-based inclusive end. So the end semantics agree; only the start
   boundary differed, and the *current* locus' own deletion is present in the set.
8. Difference summary (measured, not assumed):
   * interval: identical (`[pos, pos+len(REF))` ≡ GATK `[start, start+deletionSize]`).
   * deletion record ordering: GATK = ordered "emitted at an earlier locus"
     state; native = unordered global pre-pass **including the current locus**.
   * predicate: GATK `start' < pos`; native (before) `start' <= pos`.
   * **residual, still divergent** (unchanged by this round): GATK records only
     deletions it actually **emitted**, native counts every input record with a
     concrete deletion allele even when GATK dropped that record or pruned that
     allele (fixtures C and D, §1.2).

## 3. STEP 3 — the fix

`fastgatk-native/src/genotype_gvcf_tool.cpp`, one line of behaviour plus a
comment (`git diff`, full patch in `.diag/star-ownership-fix.diff`):

```diff
+// GATK only honours a spanning deletion that some *previously emitted* deletion
+// owns: GenotypingEngine.isVcCoveredByDeletion() requires
+// ``loc.getStart() < vc.getStart() && vc.getStart() <= loc.getEnd()``
+// (GenotypingEngine.java:365-371) and the deletion is recorded only after the
+// current locus has been subsetted (GenotypingEngine.java:178-179 calls
+// recordDeletions() *after* calculateOutputAlleleSubset()).  A deletion allele of
+// the very locus being genotyped therefore never owns that locus' own '*', which
+// is why the boundary test below is strict on the start (``span.begin <
+// record.pos``) rather than inclusive.
 bool spanning_deletion_supported_at(const Record& record,
                                     const std::vector<VariantSpan>& deletion_spans) {
     if (record.rid < 0) return false;
     for (const auto& span : deletion_spans) {
-        if (span.rid == record.rid && span.begin <= record.pos && record.pos < span.end)
+        if (span.rid == record.rid && span.begin < record.pos && record.pos < span.end)
             return true;
     }
     return false;
 }
```

Why `<` and not "skip the record itself": the two are equivalent for fixture R,
but GATK's rule is a property of the **locus**, not of the record —
`recordDeletions()` runs after the whole merged locus has been subsetted, so a
deletion carried by a *sibling* record at the same start also cannot own the
locus (`GenotypeGVCFsEngine.java:349-361` merges same-start records in dense
mode). The strict start expresses that directly, keeps the change to one
comparison, and leaves `deletion_spans` (an untouched, purely pre-computed
structure) alone.

Not restructured: the genotyping engine, the AF/PL Kokkos kernels and their ABIs,
the union/merge stages, `max_alternate_alleles`, dense-mode materialization of
REF-only calls, and every other tool. The change is confined to one predicate in
one translation unit. Native's `orphan_spanning_deletion` flag and the
`orphan_spanning_deletion_loci` telemetry keep their existing meaning; the
reported fixture now correctly counts as one orphan locus.

Reachability (read, not assumed): both call sites short-circuit on
`!options.gatk_annotation_compatibility` before consulting the predicate
(`:5557-5559` aggregate/stream union, `:6211-6213` stream-by-locus union), so the
predicate is only reached under `--gatk-compatible-annotations`. The
native-diagnostic profile that `verify_genotype_gvcf.py` also exercises is
therefore structurally untouched by this change; the reported fixture's
`orphan_spanning_deletion_loci` telemetry goes 0 → 1 only in the compat profile.

## 4. STEP 2 — the strict gate

`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`
(+172/-12 vs. the previous round's file):

* `surviving-deletion-alt` **promoted from REPORTED ONLY to GATED**
  (`gated: True`, `expect=[GATK_DEL_ROW]`), with its `why` replaced by the §2
  citations.
* Three new **GATED** isolating cases, each with its measured GATK rows in
  `expect`:
  * `surviving-deletion-alt-noncovering-upstream` — deletion emitted at 2 spans
    2-3, `*` at 4: prunes, proving the record's own deletion does not own it even
    with an upstream deletion present elsewhere;
  * `surviving-deletion-alt-covering-upstream` — **positive control**: the
    upstream deletion now spans 2-4 and *does* cover locus 4, so GATK keeps
    `*,A` and publishes `GT 1/2`. Without this case an over-broad "prune every
    `*`" fix would look correct;
  * `surviving-deletion-alt-only-alt-implausible` — pruning `*` leaves a single
    implausible ALT, the site turns monomorphic and **both** tools must write no
    record.
* One new **REPORTED ONLY** case, `covered-star-only-record`, documenting a
  *different* divergence found while measuring (§8).
* `HEADER` gained `##FORMAT=<ID=GQ,…>`: it was missing, and GATK aborts with
  `Key GQ found in VariantContext field FORMAT … isn't defined in the VCFHeader`
  (htsjdk `VCFEncoder.fieldIsMissingFromHeaderError`) on any fixture whose
  deletion record carries an informative PL row. The three pre-existing case
  expectations were re-measured after that change and are unchanged.
* The module docstring gained a "The deletion-ownership cases" section recording
  the rule and the native defect, so the gate is self-describing.

### 4.1 Exit status and literal "before" rows

**Before the fix: exit 1**, `"status": "divergence"`, 3 violations, log
`.diag/star-ownership-gate-before.log`:

```
[surviving-deletion-alt] row 0 is not byte-identical:
  GATK  ='chr1\t2\t.\tAA\tA\t82.19\t.\tAC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11\tGT:AD:DP:PL\t0/1:0,20:20:0,0,0'
  NATIVE='chr1\t2\t.\tAA\t*,A\t82.19\t.\tAC=1,1;AF=0.500,0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1,1;MLEAF=0.500,0.500;QD=4.11\tGT:AD:DP:GQ:PL\t1/2:0,0,20:20:99:100,100,100,100,0,100'
[surviving-deletion-alt-noncovering-upstream] row 1 is not byte-identical:
  GATK  ='chr1\t4\t.\tAA\tA\t82.19\t.\tAC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11\tGT:AD:DP:PL\t0/1:0,20:20:0,0,0'
  NATIVE='chr1\t4\t.\tAA\t*,A\t82.19\t.\tAC=1,1;AF=0.500,0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1,1;MLEAF=0.500,0.500;QD=4.11\tGT:AD:DP:GQ:PL\t1/2:0,0,20:20:99:100,100,100,100,0,100'
[surviving-deletion-alt-only-alt-implausible] record count differs: GATK=0 native=1
  NATIVE='chr1\t2\t.\tAA\t*\t0\t.\tAC=0;AF=0.00;AN=2;DP=20;ExcessHet=0.0000;MLEAC=0;MLEAF=0.00\tGT:AD:DP:GQ:PL\t0/0:10,0:20:99:0,100,100'
```

The pre-existing gated cases (`default`, `include-non-variant-sites`,
`stand-call-conf-50-…`, the four `surviving-concrete-alt…`) and the new positive
control were green **before** the fix, so the three violations isolate exactly the
ownership boundary.

**After the fix: exit 0**, `"status": "pass"`, 0 violations, log
`.diag/star-ownership-gate-after.log`. Case census (11 gated, 1 reported only):

| case | gated | GATK rows | native rows |
| --- | --- | --- | --- |
| `default` | yes | 0 | 0 |
| `include-non-variant-sites` | yes | 1 | 1 |
| `stand-call-conf-50-include-non-variant-sites` | yes | 1 | 1 |
| `surviving-concrete-alt` | yes | 1 | 1 |
| `surviving-concrete-alt-star-not-first` | yes | 1 | 1 |
| `surviving-concrete-alt-reversed-source-gt` | yes | 1 | 1 |
| `surviving-concrete-alt-phased-source-gt` | yes | 1 | 1 |
| **`surviving-deletion-alt`** | **yes** | 1 | 1 (byte-identical) |
| **`surviving-deletion-alt-noncovering-upstream`** | **yes** | 2 | 2 |
| **`surviving-deletion-alt-covering-upstream`** | **yes** | 2 | 2 |
| **`surviving-deletion-alt-only-alt-implausible`** | **yes** | 0 | 0 |
| `covered-star-only-record` | no | 1 | 2 |

## 5. STEP 4 — stale assertions

**No registered test asserted the old behaviour, so no test line needed
correcting.** The only test-script change in this round is the oracle itself
(§4). Verified, not assumed:

| File | Checked | Outcome |
| --- | --- | --- |
| `verify_genotype_gvcf.py:369-379` (star comment), `:390-400` (`assert star_records == []`), `:407-417` (dense literal row), `:418-438` (`--gp-qual` run) | its fixture is `A *,G,<NON_REF>` with **no deletion allele**, so `record_has_concrete_deletion()` is false and there is no deletion span either way; the comment already describes GATK's "no *emitted* deletion covers the locus" rule correctly | unchanged and passing (`{"status": "pass", "output_records": 1}`) |
| `verify_genotype_gvcf_spanning_deletion_gatk_oracle.py:53-60` (`supported`, asserted at `:103` as `"C,*" in ALT`) | its deletion is at 17:69066 (0-based 69065) spanning 69066-69067 and the `*` is at 69067 (0-based 69066) — `69065 < 69066`, so the strict test still supports it | unchanged and passing; it is the one registered fixture that sits on the new boundary, so its GATK-side assertion (`"C,*" in expected[1]`) was re-run — `fastgatk-genotype-gvcf-spanning-deletion-gatk-oracle` passed in b2 (16.53 s) with identical GATK and native rows |
| `…_max_alternate_alleles_…`, `…_multisample_reference_confidence_…`, `…_include_non_variant_…`, `…_assignment_…`, `…_exclude_intervals_…`, `…_gp_input_…`, `…_starts_in_intervals_…`, `…_multiallelic_…`, the `verify_gatk_genotype_gvcf*.py` family | grep for `*` in an ALT field: the only matches are Python `*args` splats and unrelated fixtures; no fixture carries both a `*` and a deletion allele | no change needed |
| `fastgatk-native/scripts/benchmark_genotype_gvcf.py:199` | reads `telemetry.orphan_spanning_deletion_loci` only for reporting | no asserted value, unaffected |
| C++/ctest sources | `grep -rn orphan_spanning_deletion\|spanning_deletion_supported` over `*.cpp/*.hpp/*.py/*.sh/*.txt/*.cmake` | only `genotype_gvcf_tool.cpp` and the two scripts above reference this logic; no unit test pins it |

## 6. STEP 5 — gate results (mandatory order)

| Step | Command | Result |
| --- | --- | --- |
| a | `python3 fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` | **exit 0** (11 gated, 0 violations); **exit 1 before the fix** (3 violations) |
| b1 | `python3 fastgatk-native/scripts/verify_genotype_gvcf.py` | `B1_EXIT=0` — `{"status": "pass", "output_records": 1}` |
| b2 | `ctest --test-dir fastgatk-native/build -R 'genotype-gvcf' -V` | `B2_EXIT=0` — **15/15 passed**; includes `fastgatk-genotype-gvcf-spandel-gatk-oracle` (117.28 s) and the sibling `fastgatk-genotype-gvcf-spanning-deletion-gatk-oracle` (16.53 s) |
| c | the 17-name strict-gate filter from the task | `C_EXIT=0` — **19/19 matched tests passed** |
| d | `run_regression.sh --label star-ownership` | see §6.1 |

Both trees were rebuilt **after** the last source edit and before every gate
above: `genotype_gvcf_tool.cpp` mtime `1789149727` <
`build/fastgatk-genotype-gvcf` `1789149738` <
`build-serial/fastgatk-genotype-gvcf` `1789149874`; the oracle script
(`1789149599`) was frozen before all of them. No source or test-script edit was
made after that point, and `md5sum -c .diag/star-ownership-binaries.md5` confirms
both binaries still match the ones the gates ran on
(`04cd234f3f56514f627288254bd060b7` omp, `6ce22c679e93bb3ad1455214330a1dd1`
serial).

Logs: `.diag/star-ownership-gate-before.log`,
`.diag/star-ownership-gate-after.log`, `.diag/star-ownership-verify-gvcf.log`,
`.diag/star-ownership-ctest-genotype-gvcf.log`,
`.diag/star-ownership-strict-gates.log`, `.diag/star-ownership-regression.log`.

Step (c) matched 19 tests (the task's 17 patterns also match
`fastgatk-hc-ploidy-window-invariance-gatk-oracle` and a second `mutect2-recheck`
test), all passed: window-invariance, ploidy-window-invariance, alleles-overlap,
span-del-qual, gvcf-symbolic-prior, arbitrary-ploidy-span-del-prior,
polyploid-gvcf-span-del-prior, spanning-prior-genotype-gq, af-zero-format,
multialt-owner-annotation, mutect2-recheck ×2, gvcf-indel-end, culprit,
select-variants-refonly, asfilterstatus, variant-filtration-flag-only,
droplowqual, spandel.

### 6.1 Full double-backend regression

`fastgatk-native/scripts/run_regression.sh --label star-ownership`
(omp + serial, `FASTGATK_REQUIRE_GATK_ORACLE=1` by default, ctest parallelism 8):

```
| 后端   | 构建目录                               | 结果 | 通过/总数 | 耗时  |
| omp    | OpenMP (fastgatk-native/build)         | 通过 | 300/300   | 1292.59 sec |
| serial | Serial (fastgatk-native/build-serial)  | 通过 | 300/300   | 1285.31 sec |
```

`REG_EXIT=0`; the runner's staleness check reported no warning for either build
directory. Evidence block `.diag/regression/20260912-022834/{evidence.md,omp.log,serial.log,omp.status,serial.status}`
(suite log `.diag/star-ownership-regression.log`).

**No edit of any kind — source, script or comment — was made after this run.**
Both source files were already frozen before the gates in §6 and their mtimes
(`verify_genotype_gvcf_spandel_gatk_oracle.py` `1789149599`,
`genotype_gvcf_tool.cpp` `1789149727`) still precede the binaries
(`1789149738` omp, `1789149874` serial); `md5sum -c .diag/star-ownership-binaries.md5`
re-confirms `04cd234f3f56514f627288254bd060b7` (omp) and
`6ce22c679e93bb3ad1455214330a1dd1` (serial) after the suite, i.e. the suite ran on
exactly the binaries described everywhere above.

## 7. Tree state

* Change **kept in the tree** (no commit, no branch).
* Modified: `fastgatk-native/src/genotype_gvcf_tool.cpp` (+10/-1: one comparison
  plus a comment), `fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`
  (+172/-12). Scratch added under `.diag/` only.
* Not touched: `fastgatk-native/CMakeLists.txt`, root `*.md`, Mutect2, every
  other tool, and every registered test script other than the oracle.
* Both build trees rebuilt after the last source edit, and every gate in §6 ran
  on those binaries; nothing was edited afterwards.
* `git status --short` (working tree, no commit, no branch): exactly 2 modified
  files and nothing else — `.diag/` is git-ignored (`.gitignore:42`), so the
  scratch (`round-star-ownership.md`, `star_ownership_probe.py`, the logs,
  `star-ownership-fix.diff`, `star-ownership-binaries.md5`,
  `regression/20260912-022834/`) does not appear.

## 8. What remains unproven / open

Measured this round but **not fixed and not gated**:

* **A star-only record is still written by native where GATK drops it.** With a
  genuinely covering upstream deletion (`2 AA A` + `3 A *,`), ownership is
  correct in both tools, but `GenotypeGVCFs.apply()` discards the record because
  `GATKVariantContextUtils.isSpanningDeletionOnly()` is true
  (`GenotypeGVCFs.java:327-328`, `GATKVariantContextUtils.java:2089-2091`),
  whereas native writes `chr1 3 . A * 0 . AC=0;AF=0.00;AN=2;…`. Carried as the
  new REPORTED-ONLY case `covered-star-only-record`. This is a record-emission
  rule, a different root cause from ownership, and it was already divergent
  before this round (fixtures C, D, O, X in §1.2).
* **Native still counts deletions GATK never emitted.** Fixture D: GATK drops the
  upstream deletion record (its allele fails the threshold, `GenotypeGVCFsEngine.java:157-158`
  and `:170-171` return `null` before `recordDeletions()` can run), so the
  downstream `*` is an orphan and GATK writes nothing — native still treats it as
  owned. Native's `deletion_spans` is built from the *inputs*, not from *emitted*
  output. *Speculation*: fixing this properly needs an ordered per-locus
  "emitted deletions" state threaded through the walk (the GATK model), which is
  a real restructuring and was correctly out of scope for a one-line round; it is
  currently masked in the default profile by the star-only-record rule above.
* **Dense mode (`--include-non-variant-sites`) ownership** was not re-measured
  for the deletion fixtures; the four new gated cases all use default options.
  *Speculation*: the same predicate governs it, but that is untested.
* **Incidental, unrelated divergence measured while probing** (fixture X):
  for `2 AAA AA,<NON_REF>` GATK publishes the reverse-trimmed `chr1 2 . AA A`
  (`GATKVariantContextUtils.reverseTrimAlleles`,
  `GenotypeGVCFsEngine.java:167`) while native publishes `chr1 2 . AAA AA`.
  Identical before and after this round's change, not gated, not fixed here; it
  is an allele-trimming question with no connection to ownership.
* **Two input records starting at the same locus** (the dense-mode merge path,
  `GenotypeGVCFsEngine.java:349-361`) is unmeasured; the `<` choice is argued from
  the source, not demonstrated.
* **Multi-sample / multi-shard** fixtures for this rule are unmeasured; all cases
  here are single-sample, single-input.
* The pre-existing, still-open items from the previous two rounds are untouched:
  `FILTER=LowQual` is not implemented in native GenotypeGVCFs, `QUAL=Infinity` is
  rendered `inf`/`0`, `-all-sites` is not accepted as a short alias, the
  `genotypeLikelihoods == null` arm of PREFER_PLS is not implemented, and a
  `./.` source call may still produce a record where GATK writes none.
