# Round: dense-mode record materialization — spanning-locus synthesis + reverse trim

Scope: the two behaviours the previous round left *measured but ungated*
(`.diag/round-emitted-ownership.md` §7). The hand-off framed both as dense-mode
(`--include-non-variant-sites`) behaviours; measurement shows the first one is
dense-only and the second one is not:

1. a row for every reference position covered by an input record's span — for a
   deletion allele that row carries the symbolic `*`;
2. the reverse trim of the merged record (`AA *` → `A *` in the previous
   round's measurement).

## 0. Verdict

* **(1) is STRUCTURAL.** GATK visits *every reference locus* that any input
  record spans (`VariantLocusWalker.java:150-176`), not only the loci where a
  record starts. Native's locus list is built from record starts (plus the
  per-position expansion of `<NON_REF>` blocks that dense mode already does), so
  the missing rows are missing *loci*, not missing fields: new records have to
  be produced at positions where native currently produces none, and the pass
  has to sit in **both** traversal implementations. Per the round's rule
  ("if it requires restructuring the merge or traversal, STOP"), this round
  **stops at step 4** and reports the plan (§4).
* **(2) is NOT structural — it is a contained, mode-independent allele rewrite**
  — but it is also *not* dense-only, which the hand-off framing got wrong:
  measured in **default** mode (`AAAA/AACA` → `AAA/AAC`). Isolated from (1) it
  could be implemented in ~50-80 lines, but it cannot make this round's gate
  green on its own, so nothing was implemented (§4.2).
* **Newly measured (REPORTED ONLY, not this round's scope):** (a) the
  materialization is *broader than `*`* — GATK also publishes a REF-only no-call
  row for positions covered by a spanning record with no deletion allele
  (§1.5); (b) native publishes a `*`-only row in **default** mode where GATK
  refuses to build the call (`GenotypingEngine.java:173-175`, §1.5); (c) the
  synthesized/`*`-only rows carry `QD=-0.00` where native writes `0.00`; (d) a
  sample with no record at a locus renders as `./.` in GATK and `./.:.:.:.:.` in
  native.
* **Production change: NONE.** Constraint: no `CMakeLists.txt`, no root `.md`,
  no `third_party/`. The tree gains exactly one new file: the gate (§3).
* Gate: `fastgatk-native/scripts/verify_genotype_gvcf_dense_materialize_gatk_oracle.py`
  — 11 cases (9 gated + 2 REPORTED ONLY); **exit 1 with 7 violations** in strict
  mode, exit 0 under `--expect-divergence`. Left failing **by design** (§5a).
* Step 5: (b1) exit 0, (b2) 18/18, (c) 23/23, (d) skipped — no production code
  changed (§5).

## 1. STEP 1 — characterisation

All fixtures are 100 bp single-contig poly-A `chr1` (the registered gates'
reference), single sample `STAR`, dense run =
`--include-non-variant-sites`. Literal rows are the output data rows; probes and
logs are `.diag/dense_materialize_probe{,2,3,4,5,6}.py` and
`.diag/dense-materialize-probe{,2,3,4,5,6}.log`.

### 1.1 WHEN GATK synthesises the row: the traversal, not the writer

| gatk-source `file:line` | what it does |
| --- | --- |
| `engine/VariantLocusWalker.java:153-176` | `traverse()` walks shards with `ShardedIntervalIterator(... getDrivingVariantCacheLookAheadBases())`, and for a shard containing any variant iterates **every locus** via `getLocusStream()` (`IntervalLocusIterator`, `:188-190`), calling `apply(locus, overlappingVariants, ...)` whenever `drivingVariants.query(locus)` is non-empty |
| `engine/VariantWalkerBase.java:39` | `DEFAULT_DRIVING_VARIANTS_LOOKAHEAD_BASES = 100_000` — `GenotypeGVCFs` does **not** override `getDrivingVariantCacheLookAheadBases()` (`:146-148`), so a record starting up to 100 kb upstream is still returned by `query(locus)` |
| `walkers/GenotypeGVCFs.java:320-330` | `apply(loc, variants, ...)`; `forceOutput = includeNonVariants \|\| inForceOutputIntervals`; the record is written when `forceOutput \|\| !isSpanningDeletionOnly(vc)` |
| `walkers/ReferenceConfidenceVariantContextMerger.java:150-151` | `isSpanningEvent = loc.getStart() != vc.getStart()`; a spanning record's alleles go through `replaceWithNoCallsAndDels(vc, doSomaticMerge)` instead of `remapAlleles` |
| `:222-244` | `replaceWithNoCallsAndDels()`: ref → `NO_CALL`; `NON_REF` kept; **every ALT with `allele.length() < vc.getReference().length()` → `SPAN_DEL`**; everything else → `NO_CALL` |
| `:297-303` | `filterAllelesForFinalSet()` drops `NON_REF`, the reference, symbolic-vs-symbolic and `NO_CALL` — leaving exactly `*` |
| `:325-345` | `collectTargetAlleles()` keeps `*` when a spanning deletion was seen and either a non-spanning event was seen or `removeNonRefSymbolicAllele` is false |
| `walkers/GenotypeGVCFsEngine.java:128,136` | `callRegion()` → `getVariantSubsetToProcess()` → `merger.merge(variantsToProcess, loc, ref.getBase(), true, false)` (NON_REF removed) |
| `:349-367` | **DENSE-ONLY precedence**: with `includeNonVariants`, if any variant starts exactly at the locus only those are processed ("starting record wins"), otherwise all overlapping (spanning) variants are |
| `:371-384` | `createMinimalArgs(forceOutput)` sets `OutputMode.EMIT_ALL_ACTIVE_SITES`, which is what lets a locus whose only ALT is `*` survive `GenotypingEngine.java:173-175` |
| `:191-194` | the monomorphic dense path (`cleanupGenotypeAnnotations(result, true, false)` + annotation engine) |

So the synthesis is done **during merging of spanning events** (the
`ReferenceConfidenceVariantContextMerger` path above), is driven by the
**locus-level traversal**, and is **published only in dense mode** — confirmed
by measurement (§1.2, `default-mode-no-synthesis`).

### 1.2 The literal rows (pinned GATK 4.6.2.0 vs native `build/`)

Case A — `chr1 2 . AAA A,<NON_REF>` alone, dense. GATK 3 rows, native 1:

```
GATK   chr1 2 . AAA A 92.60 . AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63   GT:AD:DP:GQ:PL 0/1:0,20:20:99:100,0,100
GATK   chr1 3 . A   * 0     LowQual AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=-0.00 GT:AD:DP:GQ:PL 0/1:0,20:20:99:100,0,100
GATK   chr1 4 . A   * 0     LowQual AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=-0.00 GT:AD:DP:GQ:PL 0/1:0,20:20:99:100,0,100
NATIVE chr1 2 . AAA A 92.60 . (row 2 only)
```

The synthesized row's **full column set** is therefore an ordinary regenotyped
record: `QUAL=0`, `FILTER=LowQual`, `INFO=AC;AF;AN;DP;ExcessHet;MLEAC;MLEAF;QD`
(no MQ/SOR/rank-sums — no reads in the pipeline), `FORMAT=GT:AD:DP:GQ:PL` with
the *spanning record's own projected call* (GT `0/1`, AD `0,20`, PL
`100,0,100`). `QD` is the literal **`-0.00`** (GATK's QUAL is `-0.0`).

Case C — the previous round's shape `2 AAA A` + `4 AA *,<NON_REF>` +
`5 A *,G,<NON_REF>`, dense. GATK 4 rows, native 3:

```
GATK   chr1 4 . A  * 0 LowQual AC=2;AF=1.00;AN=2;DP=20;ExcessHet=0.0000;MLEAC=2;MLEAF=1.00;QD=0.00 GT:AD:DP:GQ:PL 1/1:0,20:20:99:100,100,0
NATIVE chr1 4 . AA * 0 LowQual AC=2;AF=1.00;AN=2;DP=20;ExcessHet=0.0000;MLEAC=2;MLEAF=1.00;QD=0.00 GT:AD:DP:GQ:PL 1/1:0,20:20:99:100,100,0
```

i.e. the only difference at locus 4 (whose own record wins the dense
priority rule) is the trimmed `REF`, and locus 3 is missing entirely. Native
renders the `*`-only locus *byte-identically* when the record exists in the
input (`2 AA A` + `3 A *,<NON_REF>`, dense: 2 rows on both sides, §3
controls) — so the missing piece really is record supply, not row rendering.

**The projection is index-preserving for the common shape** (probe 6,
`.diag/dense-materialize-probe6.log`): with the source record
`2 AAA A,<NON_REF>` and `GT:DP:AD:PL 0/1:20:10,10,0:100,0,100,100,100,100`
(REF AD non-zero), GATK's row at 3 is `GT:AD:DP:GQ:PL 0/1:10,10:20:99:100,0,100`
— i.e. source allele 0 (`AAA`) → target allele 0 (the ref base `A`) keeping its
AD and PL row, source allele 1 (`A`) → the `*`, `<NON_REF>` → `<NON_REF>`. A
homozygous `1/1` source stays `1/1`, and a `./.` source is called by PREFER_PLS.
So for a record with a single deletion ALT the projection is the **identity map
on allele indices**; only records with several ALTs (GATK maps every non-deletion
ALT to `NO_CALL` and drops it, `:222-244` + `:297-303`) need a deduplicating
re-index (`newToOldGenotypeMap`).

### 1.3 The reverse-trim rule

| gatk-source `file:line` | what it does |
| --- | --- |
| `walkers/GenotypeGVCFsEngine.java:160-169` | reached for every locus with `originalVC.isVariant() && DP > 0` and (`isProperlyPolymorphic(regenotypedVC) \|\| includeNonVariants`) — i.e. **both modes**; `:167` `reverseTrimAlleles(withAnnotations)` runs after `finalizeAnnotations` (`:165`) and **before** the site-annotation engine (`:189`) |
| `utils/variant/GATKVariantContextUtils.java:1443-1445` | `reverseTrimAlleles(vc)` = `trimAlleles(vc, trimForward=false, trimReverse=true)` |
| `:1458` | **guard**: return unchanged if `getNAlleles() <= 1` or any allele has `length()==1 && !equals(SPAN_DEL)` — this is why SNPs and HaplotypeCaller's minimally represented indels are never trimmed |
| `:1462-1467` | the trimmed group is the non-symbolic, non-`*` alleles; `shifts = AlignmentUtils.normalizeAlleles(sequences, ranges, 0, true)`; only `endTrim = shifts.getRight()` survives (forward trimming is off), `startTrim = -shifts.getLeft()` is used only by the emptiness guard |
| `:1469-1475` | `restoreOneBaseAtEnd = emptyAllele && startTrim == 0` (and the mirror for the start); `endBasesToClip = restoreOneBaseAtEnd ? endTrim - 1 : endTrim` — **an allele can never be emptied** |
| `:1489-1515` | `trimAlleles(vc, fwdTrimEnd, revTrim)`: `:1492-1493` no-op when `fwdTrimEnd == -1 && revTrim == 0`; `:1497-1499` `*` and symbolic alleles are copied untouched; `:1512` `builder.stop(start + alleles.get(0).length() - 1)` |
| `:1516-1528` | genotypes are remapped by identity (`AlleleMapper`), so GT indices, phase and PLs are unchanged |

Measured:

| fixture (mode) | GATK | native |
| --- | --- | --- |
| `2 AAAA AACA,<NON_REF>` (default) | `chr1 2 . AAA AAC 92.64 . …` | `chr1 2 . AAAA AACA 92.64 . …` |
| `2 ACGTACGT ACGT,<NON_REF>` (default, `ACGT` reference) | `chr1 1 . ACGTA A 92.60 . …` | `chr1 1 . ACGTACGT ACGT …` |
| `2 CGTACG CGT,<NON_REF>` (default) | `chr1 2 . CGTACG CGT …` (unchanged) | identical |
| `2 AAA A,<NON_REF>` (default) | unchanged (ALT is 1 base → guard `:1458`) | identical |
| `4 AA *,<NON_REF>` (dense) | `chr1 4 . A * …` | `chr1 4 . AA * …` |

* How many bases: the common trailing run over the non-symbolic, non-`*`
  alleles (`endTrim`), minus one when the clip would consume an allele and no
  leading bases are clipped (`:1469-1475`). One base in both `AAAA/AACA` and
  `AA/*`; three in the `ACGTACGT/ACGT` corner.
* Can it empty an allele: **no** — the `:1469-1475` guard restores one base;
  measured in the `ACGTACGT/ACGT` corner where `endTrim` is 4 and the ALT keeps
  the single base `A`.
* Which records: every emitted record with ≥2 alleles and no single-base
  non-`*` allele — in practice multi-base REF with all-multi-base ALTs (MNPs
  and `*`-records with a multi-base REF). Multi-sample merges whose ALT union
  contains any 1-base ALT are protected by the guard, which is why the
  registered corpus (HaplotypeCaller-minimal alleles) never sees it.
* Placement constraint: GATK applies the trim **after** `recordDeletions()`
  (`GenotypingEngine.java:179`, inside `calculateGenotypes` at
  `GenotypeGVCFsEngine.java:156`) and **before** `annotateContext` (`:189`). A
  native implementation must keep that order, otherwise the
  `*`-records' own `deletionSize` (previous round's case T) would change.

### 1.4 Do (1) and (2) interact?

Measured/derived answer: **only in one harmless direction.**

* Coverage decisions are made on the **input** records' spans inside the merger
  (`loc.getStart() != vc.getStart()`, `ReferenceConfidenceVariantContextMerger.java:150`),
  long before the trim at `:167`, so a trimmed record cannot create or delete a
  locus to visit. The trim cannot "produce records that then need synthesis".
* Synthesized rows are themselves subject to the trim, but their merged REF is
  one base, so the `:1458` guard short-circuits — measured: GATK's synthesized
  rows at 3 and 4 keep `REF=A`.
* The one thing a trim would break if placed wrongly is the previous round's
  `EmittedDeletions` interval (`AA *` has `deletionSize=1`, `A *` has `0`), which
  is why the placement constraint above matters.

### 1.5 Three further measured divergences (REPORTED ONLY, out of scope)

| # | fixture | GATK | native |
| --- | --- | --- | --- |
| a | dense, `2 AAAA AACC,<NON_REF>` (spanning substitution, no deletion) | 4 rows: the site + REF-only `A . . . DP=20 GT:AD ./. :0` at 3, 4, 5 | 1 row |
| b | **default**, `2 AAA A,<NON_REF>` + `3 A *,<NON_REF>` (owned `*`) | 1 row (locus 3 refused, `GenotypingEngine.java:173-175`) | 2 rows — native publishes `chr1 3 . A * 0 LowQual …` |
| c | dense, `2 AA A,<NON_REF>` + `3 A *,<NON_REF>` | `…QD=-0.00…` | `…QD=0.00…` |
| d | two samples, S1's deletion at 2 spans 3, S2 has `3 A G` | S1 column `./.` | `./.:.:.:.:.` |

(a) shows the gap is **wider than `*`**: the rule is "every covered position",
and a covered position whose spanning record contributes no `*` is published as a
REF-only no-call row (two shapes: `QUAL 192.21 … GT ./.` when a pruned `*`
survived the merge, and `QUAL . … GT:AD ./. :0` when the merge returned null).
(b) is a **second, previously unmeasured** call-site gap: native's
`apply_gatk_output_allele_subset()` refuses only the *empty*-ALT-set case
(`output_alleles.size() == 1`, `genotype_gvcf_tool.cpp:3053-3073`), not the
*exactly-`*`* case. The registered `covered-star-only-record` fixture reaches
the empty-ALT-set path, which is why the corpus is green.

### 1.6 Reachability on production input

`gatk-source/src/test/resources/.../haplotypecaller/expected.testGVCFMode.gatk4.g.vcf`
(GATK's own pinned HaplotypeCaller gVCF output, 1291 records, chr20) contains
**139 positions** that are covered by a spanning variant record's span, are not
a record start, and are not inside a `<NON_REF>` block — every one of them
produced by a record carrying a deletion allele (e.g.
`20:10068158 GTGTATATATATA > G,<NON_REF>` covering 10068159-10068170 behind the
`*` record at 10068160; `20:10004769 TAAAACTATGC > T,<NON_REF>`;
`20:10097436 CTTTTCTTT… > C,<NON_REF>`). In dense mode GATK materializes a `*`
row at each; native materializes none, because its locus set is exactly
{record starts} ∪ {positions inside dense-expanded blocks}
(`genotype_gvcf_tool.cpp:6611-6690` for the block expansion, `:6785-6800` for
splitting + grouping). Scan: `.diag/real_gvcf_coverage_scan.py`.
(The repo's other real gVCFs are git-lfs pointers, so no second corpus.)

## 2. STEP 2 — blast radius

**Structural or additive?** For (1): **structural**, in the sense the task
defines it — rows must appear at loci where native currently emits none, so the
record-supply/traversal changes in *both* traversal implementations (aggregate
`genotype_gvcf_tool.cpp:6370-7000`, streaming `:5230-6350`). It is *contained*,
though: no stage, no kernel, no ABI, no reordering, no change to
`apply_gatk_output_allele_subset()`, and the nearest existing analogue is already
in the file (the dense reference-block expansion at `:6611-6690` and
`emit_dense_stream_reference_piece()` at `:5320`). For (2): **additive**, a local
allele rewrite in the compute stage.

**Which registered outputs change?** None — measured, not assumed:

* `ctest -N -R 'genotype-gvcf'` → **18 tests** (spandel, malformed,
  header-order, lowqual, contract, gatk-oracle, legacy-qual, multisample,
  multiallelic, inbreeding, include-non-variant, assignment,
  multisample-reference-confidence, max-alternate-alleles, gp-input,
  exclude-intervals, starts-in-intervals, spanning-deletion).
* **Sensitive to (1)**: a test is sensitive only if a *dense* run's input has a
  deletion record covering a position where no record starts. Seven of the 18
  scripts ever pass `--include-non-variant-sites` (§contract/spandel/lowqual/
  header-order/include-non-variant/gatk-oracle/multisample-reference-confidence);
  `.diag/scan_dense_gap_cases.py` imports each of them and evaluates **every**
  case's fixture: **0 cases** contain a gap-shaped deletion (the spandel and
  lowqual dense fixtures all use `2 AA A`, whose single covered position 3 has
  its own starting record, and the three real-gVCF dense fixtures run on
  `17:69000-69100`, which `--min-alt-support 1` HaplotypeCaller finds indel-free
  — the 1001-row dense comparison there is byte-identical today,
  `.diag/dense-materialize-probe.log`/oracle baseline).
* **Sensitive to (2)**: a test is sensitive only if some *emitted* record has no
  single-base non-`*` allele and a shared trailing base.
  `.diag/scan_fixture_sensitivity.py` evaluates every inline VCF line of every
  script under `fastgatk-native/scripts/`: **0 TRIM hits** in the
  `genotype-gvcf` family (the only TRIM hits are in `verify_left_align.py`
  and `verify_reblock_gvcf.py`, i.e. other tools, structurally unreachable from
  `genotype_gvcf_tool.cpp`).
* No C++ unit test references the functions a fix would use (`grep` over
  `*.cpp/*.hpp/CMakeLists.txt`: `genotype_gvcf_tool.cpp` only).

**Cheap ones run before changing anything** (this round changed nothing, so
"before" and "after" coincide): the whole `genotype-gvcf` ctest subset and the
23-test strict-gate subset were run on the unchanged tree (§5).

**Consequence:** the registered corpus *cannot* validate a fix for (1) — it is
green precisely because it contains no gap-shaped fixture. Only a new gate
(§3) can.

## 3. STEP 3 — the gate

`fastgatk-native/scripts/verify_genotype_gvcf_dense_materialize_gatk_oracle.py`
(new file, 11 cases = 9 gated + 2 REPORTED ONLY; not registered in
`fastgatk-native/CMakeLists.txt` — constraint 1 forbids editing it and the
orchestrator registers gates).

Gated: 4 failing cases pinning the two behaviours + 5 controls that are green
today and must stay green:

| case | mode | today |
| --- | --- | --- |
| `deletion-alone-dense-materializes-covered-positions` | rows | **fail** (GATK 3 / native 1) |
| `deletion-plus-downstream-star-dense` | rows | **fail** (GATK 4 / native 2) |
| `merged-star-record-dense-reverse-trim` | rows | **fail** (GATK 4 / native 3, untrimmed REF) |
| `suffix-substitution-default-reverse-trim` | rows | **fail** (untrimmed, default mode) |
| `no-deletion-alt-no-star-rows` (control) | no-star | pass |
| `pruned-deletion-no-star-rows` (control) | no-star | pass |
| `starting-record-wins-dense` (control) | rows | pass |
| `default-mode-no-synthesis` (control) | rows | pass |
| `no-trim-needed-multi-base-alt` (control) | rows | pass |

REPORTED ONLY (GATK truth pinned, native recorded): `owned-star-only-locus-
dense-negative-zero-qual` (§1.5c) and `owned-star-only-locus-default-mode-
refused` (§1.5b).

**Before the fix: exit 1**, `"status": "divergence"`, exactly 7 violations
(`.diag/dense-materialize-oracle-before.log`):

```
[deletion-alone-dense-materializes-covered-positions] record count differs: GATK=3 native=1;
    missing native positions=['chr1:3', 'chr1:4']
[deletion-plus-downstream-star-dense] record count differs: GATK=4 native=2;
    missing native positions=['chr1:3', 'chr1:4']
[deletion-plus-downstream-star-dense] row 1 is not byte-identical:
    GATK  ='chr1\t3\t.\tA\t*\t0\tLowQual\tAC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=-0.00\tGT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100'
    NATIVE='chr1\t5\t.\tA\tG\t82.26\t.\t…\tGT:AD:DP:PL\t0/1:0,20:20:0,0,0'
[merged-star-record-dense-reverse-trim] record count differs: GATK=4 native=3; missing=['chr1:3']
[merged-star-record-dense-reverse-trim] row 1 / row 2 are not byte-identical
    GATK row 2 ='chr1\t4\t.\tA\t*\t0\tLowQual\t…\tGT:AD:DP:GQ:PL\t1/1:0,20:20:99:100,100,0'
    NATIVE row2='chr1\t4\t.\tAA\t*\t0\tLowQual\t…\tGT:AD:DP:GQ:PL\t1/1:0,20:20:99:100,100,0'
[suffix-substitution-default-reverse-trim] row 0 is not byte-identical:
    GATK  ='chr1\t2\t.\tAAA\tAAC\t92.64\t…'
    NATIVE='chr1\t2\t.\tAAAA\tAACA\t92.64\t…'
```

`--expect-divergence` exits 0 (diagnostic mode), as house convention requires.

## 4. STEP 4 — STOP (no production change) and the plan

### 4.1 (1) spanning-locus materialization — structural

Where native stands: the locus list is `std::sort`+`key` grouping of the decoded
records (`:6785-6800`), where dense mode has already expanded `<NON_REF>` blocks
one position at a time (`:6611-6690`); the streaming path has the analogue
(`emit_dense_stream_reference_piece`, `:5320-5340`, fed by
`split_stream_reference_block`, `:5232`). Nothing in either path knows that a
*non-block* record covers downstream positions.

Plan (≈150-250 lines, both traversals):

1. New pass `materialize_spanning_loci(output_header, reference_index, records,
   options)`, run **only** when `options.include_non_variant_sites`, inserted
   between `split_reference_blocks_at_variants()` (`:6785`) and the sort
   (`:6786`) in the aggregate path, and at the equivalent point of the streaming
   producer; same function reused by both.
   * collect the record starts (a `std::set` of (rid,pos)) and the spans, reusing
     `record_span_end()` (`:1213`) exactly as `split_reference_blocks_at_variants`
     does (`:3446-3452`);
   * for every position `p` of the union of spans with **no** record start and
     not already inside a dense-expanded block, append one synthetic `Record` at
     `p` whose alleles are
     `[refbase(p), "*", "<NON_REF>"]` if any covering record has an ALT shorter
     than its own REF, else `[refbase(p), "<NON_REF>"]` (measured §1.5a);
   * its sample data is GATK's reference-confidence projection
     (`AlleleSubsettingUtils.getIndexesOfRelevantAllelesForGVCF` +
     `newToOldGenotypeMap` + `generateAD`, called from
     `ReferenceConfidenceVariantContextMerger.mergeRefConfidenceGenotypes():575-612`).
     Measured (probe 6): for a record with a single deletion ALT the projection
     is the **identity map on allele indices** — source allele 0 (the long REF)
     keeps index 0 and its AD/PL row, the deletion ALT becomes the `*`, and
     `<NON_REF>` stays last — so GT/PL/AD can be copied verbatim; `Record::gt`,
     `pl`, `ad`, `dp` are already the Host vectors the existing
     `remap_record_to_allele_union()` (`:2375`) machinery consumes. Records with
     several ALTs (every non-deletion ALT becomes `NO_CALL` and is dropped) need
     a deduplicating re-index instead;
   * do **not** mark it `reference_block` (it must be genotyped, unlike a block).
2. Byte-exactness details that the gate pins and that must not be forgotten:
   * `QUAL = -0.0` for these loci → `QD` must render as `-0.00` (§1.5c). Native
     today renders `0.00` even for the *explicit*-record version of the same row,
     so this is a prerequisite, not a consequence, of (1).
   * the dense "starting record wins" precedence is already implicit in native's
     group-by-key, but the pass must skip positions that have a start, or the
     `starting-record-wins-dense` control breaks.
   * the synthesized loci, when emitted, must feed the previous round's
     `EmittedDeletions` state in traversal order exactly as GATK's
     `recordDeletions()` does; because they are inserted in sorted order before
     the pipeline, this falls out of the existing call at
     `apply_gatk_output_allele_subset()` (`:3084`).
   * the default-mode half of the star-only rule (§1.5b) should be fixed in the
     same change: `GenotypingEngine.java:173-175` refuses a locus whose only
     surviving ALT is `*` regardless of ownership, which native does not
     implement.
3. Verification: the new gate, plus the 18 `genotype-gvcf` tests and the 23-test
   strict subset, plus a double-backend `run_regression.sh` (the pass runs on
   every dense invocation, so the aggregate/streaming pair must be checked with
   `--stream-by-locus` too — the gate currently exercises the aggregate path
   only).

### 4.2 (2) reverse trim — confined, mode-independent, **not** implemented here

Plan (≈50-80 lines):

1. New helper `reverse_trim_record(bcf_hdr_t*, Record&)` implementing
   `GATKVariantContextUtils.trimAlleles(vc,false,true)`:
   guard `NAlleles<=1` or any non-`*` allele of length 1 → return (`:1458`); the
   trimmed group is the non-symbolic, non-`*` alleles; `endTrim` = common
   trailing run; `endBasesToClip = (emptyAllele && startTrim == 0) ? endTrim-1 :
   endTrim` (`:1465-1475`); `0` → return (`:1492-1493`); clip the tail of every
   non-symbolic, non-`*` allele, leave `*`/symbolic and all genotypes untouched
   (`:1497-1519`).
2. Call it in the compute stage **after** `apply_gatk_output_allele_subset()`
   and `apply_genotype_assignment()` (GATK: `:156` → `:165` → `:167` → `:189`)
   and **before** `update_site_annotations()`, in both the aggregate compute
   lambda and the streaming compute function. It must be a no-op for
   `record.finalized_monomorphic_ref` rows (their `NAlleles == 1` guard already
   guarantees that).
3. Why not now: it is independent of (1) but cannot make this gate green, it
   changes output in **both** modes (a wider claim than the hand-off's "dense
   mode"), and the round's brief is to implement only what is confined *as a
   whole*. Recommended as its own round, with its own gate cases (the two trim
   cases already exist here and can be split out).

## 5. STEP 5 — gate results, mandatory order

_No production code was changed in this round_, so (d) is skipped by its own
rule ("if you changed NO production code, say so and skip (d)").

| step | command | result |
| --- | --- | --- |
| a | `python3 fastgatk-native/scripts/verify_genotype_gvcf_dense_materialize_gatk_oracle.py` | **exit 1**, `"status": "divergence"`, 7 violations in 4 of the 9 gated cases, 5 controls green — **stopped at step 4, left failing by design**. `--expect-divergence` → **exit 0**. Logs `.diag/dense-materialize-oracle-before{,-diag}.log` |
| b1 | `python3 fastgatk-native/scripts/verify_genotype_gvcf.py` | **exit 0** — `{"status": "pass", "output_records": 1}` |
| b2 | `ctest --test-dir fastgatk-native/build -R 'genotype-gvcf' -V` | **exit 0** — **18/18 passed**, 988.2 s |
| c | the 21-pattern strict-gate filter | **exit 0** — **23/23 matched tests passed**, 1873.7 s |
| d | `run_regression.sh --label dense-materialize` | **skipped** — no production code changed (the rule allows this; §6) |

Combined log `.diag/dense-materialize-step5.log` (`B1_EXIT=0`, `B2_EXIT=0`,
`C_EXIT=0`). Nothing in the tree was edited between these runs except comment and
citation text inside the *new, unregistered* gate script, which was re-run
afterwards with the same result (exit 1, 7 violations) — no production file, no
registered test and no built binary changed at any point, so no rebuild and no
byte-identity proof were required.

## 6. Tree state

* **The tree CONTAINS no production change.** Exactly one file is added:
  `fastgatk-native/scripts/verify_genotype_gvcf_dense_materialize_gatk_oracle.py`.
  Scratch (`.diag/`) is git-ignored.
* Not touched: `fastgatk-native/CMakeLists.txt`, root `*.md`, `third_party/`,
  `genotype_gvcf_tool.cpp`, every other tool, every pre-existing test script.
* No git commit, no branch.

## 7. What remains unproven / open

* The gate exercises the **aggregate** traversal only; the streaming path
  (`--stream-by-locus`) has its own dense record supply
  (`emit_dense_stream_reference_piece`, `:5320`) which the plan must extend and
  which nothing currently gates.
* Reachability was measured against GATK's pinned HaplotypeCaller gVCF
  (139 gap positions in 1291 records, §1.6); the repo's other real gVCFs are
  git-lfs pointers, so a second, larger production corpus was not available, and
  GATK's dense output on that file could not be *run* (its hg38 chr20 reference
  is not in the repository) — the row count there is derived from the traversal
  rule, not measured.
* The projection's re-indexing is measured only for the single-deletion-ALT
  shape (identity map, probe 6); the multi-ALT dedup path (a record whose ALT
  list carries several shorter ALTs, e.g. `G,*,<NON_REF>`) is read from the
  merger's `filterAllelesForFinalSet()`/`newToOldGenotypeMap` and is *argued*, not
  measured.
* The reverse trim's `normalizeAlleles()` internals (htsjdk 4.2.0) were not read
  from source (only the pinned jar is available); the rule is derived from
  GATK's `trimAlleles()` wrapper plus four measurements, including the
  cannot-empty corner.
* The `-0.00` rendering (§1.5c) is a prerequisite for byte parity on the
  synthesized rows but is a different root cause and is left REPORTED ONLY.
* Multi-sample / cross-sample dense shapes were only spot-checked (the
  cross-sample spanning case produced no extra `*` on either side, §1.5d);
  `CombineGVCFs`/GenomicsDB-produced inputs were not exercised at all.
* Not re-measured this round: the previous round's streaming-probe results and
  the `.diag/emitted-ownership-*` residuals still stand as reported.
```
