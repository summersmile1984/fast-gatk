# Round 49 — adversarial falsification attempt against `apply_gatk_reverse_trim()`

**Subject of the test.** The just-landed Host-side helper
`apply_gatk_reverse_trim()` in `fastgatk-native/src/genotype_gvcf_tool.cpp:3028-3086`, called from
three return points inside `apply_gatk_output_allele_subset()` (`:3122`, `:3159`, `:3389`), which
mirrors `GenotypeGVCFsEngine.java:160-169` → `GATKVariantContextUtils.reverseTrimAlleles()`
(`GATKVariantContextUtils.java:1443-1521` + `AlignmentUtils.normalizeAlleles()`,
`AlignmentUtils.java:818-845`).

**Oracle.** Pinned GATK 4.6.2.0 (`third_party/jdk17/bin/java -Xmx1g -jar
third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar`, htsjdk 4.2.0 per GATK's own
startup banner) against `fastgatk-native/build/fastgatk-genotype-gvcf` (and
`build-serial/…` for the registered gate), same input, same semantic arguments, native always with
`--gatk-compatible-annotations`.

**Everything below is either MEASURED (a command was run and its output pasted) or REASONED (a
reading of source code). Each claim is labelled.** No production file was modified; no commit or add
was run; scratch lives in `.diag/round49/` (plus `wide/` and `work/` beneath it) because `/tmp` does
not survive between tool calls.

---

## 0. Verdict in one page

| Question the fix claims | Result |
| --- | --- |
| Clip length = common trailing run capped by the shortest candidate | **not falsified** — 104 shapes, 0 allele differences (measured) |
| Clip-one-fewer emptiness rule when the run would empty an allele | **not falsified** — `ACGTACGT/ACGT→ACGTA/A`, 12A/11A → `AA/A` (measured) |
| Guard at `:1458` on the *post-subset* allele list | **not falsified** — one-base ALT protects the record; pruned one-base ALT does not (measured) |
| Symbolic alleles and `*` copied untouched | **not falsified** — measured in `*` and `<NON_REF>` fixtures |
| POS never moves | **not falsified** — no native-only shifted position in any fixture (measured) + helper writes no `pos` (reasoned) |
| Placement constraint: trim runs **after** `EmittedDeletions::record()` | **not falsified** — purpose-built chain fixture matches GATK at every record start, aggregate *and* `--stream-by-locus` (measured) |
| FORMAT (GT/AD/DP/GQ/PL) untouched by the trim | **not falsified** — no fixture in which a FORMAT field differed while the alleles matched, apart from out-of-scope dense row shapes (measured) |
| **Native trims exactly where GATK trims** | **FALSIFIED** — see counterexamples 1 and 2 below |
| Native's post-trim record is a faithful GATK record | **FALSIFIED in one respect** — counterexample 3 (stale `INFO/END`) |

Counterexamples (details in §6): **C1** `INFO/DP=0` loci — native trims where GATK provably does not
(in dense mode both tools publish a row at the same POS with *different alleles*); **C2** the same
gate removes GATK's emitted-deletion bookkeeping, so the two tools disagree on `*` ownership;
**C3** the trim leaves native's pass-through `INFO/END` inconsistent with the record's own span.

All three are *reach/annotation* defects around a clip algorithm that itself matched GATK on every
one of my ~140 fixture-shaped probes.

---

## 1. What was read, and what was measured

Read (source-level, for hypothesis generation only):

* `fastgatk-native/src/genotype_gvcf_tool.cpp:2986-3086` (helper + its doc comment),
  `:3088-3400` (`apply_gatk_output_allele_subset`, including the three call sites), `:1262-1300`
  (`EmittedDeletions`), `:2796-2900` (`materialize_gatk_monomorphic_ref_call`), `:3401-3560`
  (`materialize_reference_only`), `:5540-5660` (`split_stream_reference_block`,
  `emit_lazy_stream_reference_piece`), `:6277-6350` and `:6979-7075` (the two compute stages).
* `gatk-source/.../GenotypeGVCFsEngine.java:100-224`, `:345-360`,
  `.../genotyper/GenotypingEngine.java:120-200`, `:295-370`,
  `.../ReferenceConfidenceVariantContextMerger.java:128-215`, `:245-285`, `:353-362`, `:432-440`,
  `.../utils/variant/GATKVariantContextUtils.java:1443-1521`, `:489-500`,
  `.../utils/read/AlignmentUtils.java:812-845`, `.../GenotypeGVCFs.java:230-340`,
  `.../engine/VariantLocusWalker.java:95-150`.

Measured with the pinned binaries (all commands in §11, raw logs under
`fastgatk-native/evidence/2026-09-13-round49/raw/`):

1. `GATKVariantContextUtils.reverseTrimAlleles()` called directly from Java against hand-built
   `VariantContext`s compiled with the pinned GATK fat jar (`raw/TrimProbe.java`, `raw/trimprobe.out`)
   — 16 allele shapes; plus `raw/P2.java` / `raw/p2.out` for attribute behaviour.
2. A **wide sweep**: 60 loci (stage A) and 44 loci (stage B) in one gVCF each, so one GATK startup
   covers every shape (`raw/wideA.log`, `raw/wideB.log`), compared against GATK *and* against an
   independent Python transcription of the helper.
3. ~35 targeted fixtures across the eight requested axes, four drivers (`raw/probe.py`, `probe2.py`,
   `probe3.py`, `probe4.py`), raw comparison output in `raw/all-rows.txt`, per-field attribution in
   `raw/attribution-table.txt`.
4. The registered gate, `--expect-divergence` and strict, on both backends
   (`raw/gate-divergence-parallel.log`, `raw/gate-strict-build.log`,
   `raw/gate-strict-build-serial.log`): **status pass, 0 violations, exit 0 on both backends**
   (independently reproduced, not taken on trust).

Canonical invocation form used for every fixture:

```
GATK:   third_party/jdk17/bin/java -Xmx1g -jar \
          third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar GenotypeGVCFs \
          -R <ref.fa> -V <in.g.vcf> [args] -O <out.vcf> --create-output-variant-index false
native: fastgatk-native/build/fastgatk-genotype-gvcf \
          -R <ref.fa> -V <in.g.vcf> [args] --gatk-compatible-annotations -O <out.vcf>
```

Every gVCF input is indexed with `gatk IndexFeatureFile` first (GATK needs it for
`--include-non-variant-sites`); the reference is a 400 bp single-contig FASTA with a hand-written
`.fai`/`.dict`.

---

## 2. Facts about the oracle that the fixture design depends on (measured)

These change what a "GATK-faithful" implementation even means, and several of them are what the
counterexamples exploit.

**(a) GATK's input contract.** A gVCF record must carry `<NON_REF>`:
`A USER ERROR has occurred: The list of input alleles must contain <NON_REF> as an allele but that is
not the case at position 12` (wide stage A first attempt, `raw/wideA.log`); native accepts such input
and trims it. All later fixtures carry `<NON_REF>`.

**(b) `<NON_REF>` never reaches the genotyper.** `GenotypeGVCFsEngine.callRegion()` calls
`merger.merge(variantsToProcess, loc, ref.getBase(), true, false)` and
`ReferenceConfidenceVariantContextMerger.merge(..., removeNonRefSymbolicAllele=true, ...)` drops the
`<NON_REF>` allele when the target allele list is built (reasoned: `ReferenceConfidenceVariantContextMerger.java:137,195`
+ `collectTargetAlleles`). Both native paths likewise call `remove_non_ref_allele()` before the compute
stage (reasoned: `genotype_gvcf_tool.cpp:6315`, `:6979`). **I tested the consequence anyway**:
`prune-refonly-dense` (hom-ref evidence so every concrete ALT is pruned) is byte-identical —
`chr1 2 . AAAA . 127.78 . DP=20;MLEAC=.;MLEAF=. GT ./.` on *both* sides — so native does not clip a
REF-only emitted list, and no phantom `<NON_REF>` deletion span is recorded (`g-longref-phantom`,
below).

**(c) In DEFAULT mode GATK never merges two input records.** `GenotypeGVCFs.onTraversalStart()`
calls `changeTraversalModeToByVariant()` whenever `includeNonVariants` is false
(reasoned: `GenotypeGVCFs.java:296-298`), and `VariantLocusWalker.traverse()` then does
`apply(variant, Collections.singletonList(variant), …)` (reasoned:
`VariantLocusWalker.java:126-140`; measured consequence in §6/§7).

**(d) In DENSE mode GATK refuses two records that start at the same locus.**
`GenotypeGVCFsEngine.getVariantSubsetToProcess()` throws `IllegalStateException("Variant input
contains more than one variant starting at location: chr1:2-5")` (measured, `d2-merged-two-het-dense`,
which GATK exited with a USER ERROR; reasoned:
`GenotypeGVCFsEngine.java:347-364`).

**(e) The `DP>0` precondition.** `regenotypeVC()` only regenotypes when
`originalVC.isVariant() && originalVC.getAttributeAsInt(DP_KEY,0) > 0`; otherwise
`result = originalVC` and `:167` never runs (reasoned: `GenotypeGVCFsEngine.java:152-175`). The merged
VC's `INFO/DP` is `calculateVCDepth()` = **the input record's `INFO/DP` when the key is present**,
else the sum of genotype depths (reasoned: `ReferenceConfidenceVariantContextMerger.java:353-362`,
`:382`). This is the mechanism behind C1/C2, and it is measured: `INFO/DP=1` trims
(`dp1-dense` EQUAL, `AAA/AAC` both sides), `INFO/DP=0` does not (`dp0-dense` DIFFER).

**(f) htsjdk 4.2.0 gives symbolic alleles `length() == 0`.** Measured:
`Allele.create("<NON_REF>", false)` → `display=<NON_REF> len=0 baseslen=0 isSymbolic=true`
(`raw/p2.out`), while `Allele.SPAN_DEL` is `*` with `len=1 isSymbolic=false`. So the `:1458` guard's
`a.length()==1` never fires for a symbolic allele, and `recordDeletions`'s
`refLen - allele.length()` would score a span of `refLen` for any emitted symbolic allele (reasoned).
Native's `EmittedDeletions::record()` uses `std::string::size()`, i.e. 9 for `"<NON_REF>"` (reasoned:
`:1291-1298`). I could not turn this into a divergence because of (b) — see `g-longref-phantom` in §7.

**(g) A trimmed record must not carry an `END` attribute inconsistent with its stop.**
`VariantContextBuilder.make()` throws
`Badly formed variant context at location chr1:100; getEnd() was 102 but this VariantContext contains
an END key with value 103` when `reverseTrimAlleles()` is handed a VC with `END` set (measured,
`raw/p2.out`, thrown from `trimAlleles` at `GATKVariantContextUtils.java:1521`). GATK's pipeline never
hits this because `removeStaleAttributesAfterMerge()` deletes `END` (reasoned:
`ReferenceConfidenceVariantContextMerger.java:437`) — measured consequence in C3.

---

## 3. Unit-level oracle: what the shipped Java does to an allele list

`GATKVariantContextUtils.reverseTrimAlleles()` on synthetic records (`raw/trimprobe.out`); start=100
for every probe, so `start=100 … end=100+len(REF)-1` after the trim in every trimmed case:

```
multibase-ref-nonref        in=AAAA/<NON_REF>          out=A/<NON_REF>   start=100 end=100  trimmed
multibase-ref-nonref-8      in=ACGTACGT/<NON_REF>      out=A/<NON_REF>   start=100 end=100  trimmed
multibase-ref-star          in=AA/*                     out=A/*           start=100 end=100  trimmed
multibase-ref-star-nonref   in=AA/*/<NON_REF>           out=A/*/<NON_REF> start=100 end=100  trimmed
ref-alt-shared3             in=AAAA/AA                  out=AAA/A         start=100 end=102  trimmed
ref-alt-shared1             in=AAAA/AACA                out=AAA/AAC       start=100 end=102  trimmed
ref-alt-noshared            in=AAAA/AACC                unchanged=true
guard-one-base-alt          in=AAA/A                   unchanged=true    (the :1458 guard)
guard-star-only-alt         in=AAA/*                     out=A/*           (only REF is a candidate)
empty-guard-acgtacgt        in=ACGTACGT/ACGT            out=ACGTA/A       (clip 4-1)
ref-lt-alt-insertion        in=AA/AAAA                   out=A/AAA
ref-2alt-mixed-len          in=AAAAA/AAA/CAAAA          out=AAA/A/CAA
ref-star-long               in=AAAAA/*                  out=A/*
symbolic-only               in=ACGT/<DEL>                out=A/<DEL>
symbolic-plus-nonref        in=ACGT/<DEL>/<NON_REF>      out=A/<DEL>/<NON_REF>
ref-alt-symbolic-shared     in=ACGTAC/ACGT/<NON_REF>     unchanged=true    (no common trailing base)
star-ref-multibase-3        in=AAA/*/<NON_REF>           out=A/*/<NON_REF>
```

Two things this settles by measurement rather than reading: the clip is a pure trailing rewrite with
`start` unchanged, and `end` recomputed from the new REF length; and an allele list whose only
non-symbolic, non-`*` candidate is the REF is clipped **to one base** (`AAAA/<NON_REF> → A/<NON_REF>`,
`ACGTACGT/<NON_REF> → A/<NON_REF>`), i.e. the REF itself can end up one base long.

---

## 4. Wide sweep: clip length, guard and emptiness rule over 104 shapes

One GATK run + one native run per stage; every locus is a separate position with 12 bp spacing
(reference = deterministic 1600 bp ACGT sequence, each `REF` taken from the reference, each record
carrying `<NON_REF>`).

* **Stage A** (`raw/wideA.log`): lengths 2-7 × {common trailing run of exactly 0,1,2,3; trailing
  deletion of 1,2,3 bases; one-base ALT; `ins-repeat`; `ins-diff`; `ins-shared`} = **60 loci,
  0 not EQUAL, 0 disagreements between GATK and the Python transcription of the helper.**
* **Stage B** (`raw/wideB.log`): the same shapes with a *second* concrete ALT (3 alleles + `<NON_REF>`,
  `GT=1/2` so both concrete ALTs are supported) = **44 loci, 0 allele differences, 0 model
  disagreements.** All 16 "DIFFER" flags are `INFO/QD` only (`QD=4.15` GATK vs `4.16` native on the
  loci with two concrete ALTs; verified field-by-field), an unrelated rounding divergence.

Representative rows (POS, input REF/ALT, GATK = native after the trim):

```
L2-run0   12 TC/TG,<NON_REF>          -> TC/TG       (no common trailing base)
L2-run1   24 TC/AC,<NON_REF>          -> T/A
L3-run2  120 TGC/AGC,<NON_REF>        -> T/A         (run 2 < shortest 3, so clip 2)
L4-del3  276 GTTG/G,<NON_REF>         -> GTTG/G      (one-base ALT, guard)
L4-run0   12 TCGT/TCGA,TGGT,<NON_REF> -> TCGT/TCGA,TGGT   (no trailing match: unchanged)
L4-run2   36 TCGT/TGGT,TCTT,<NON_REF> -> TCG/TGG,TCT       (clip 1 from REF and both ALTs)
L5-run3  180 GATCG/GCTCG,GAACG,<NON_REF> -> GAT/GCT,GAA
L7-ins-repeat 696 ATCGACG/ATCGACGG,<NON_REF> -> ATCGAC/ATCGACG   (emptiness rule: clip 5, not 6)
L8-ins-shared 528 GCACGATC/GCACGATCTC,GGACGATC,<NON_REF> -> GCACGA/GCACGATC,GGACGA
```

**Not falsified.** No shape produced a different clip length, a different emptiness decision, or a
different guard decision between GATK and native, and my independent transcription agreed with GATK
on all 104.

---

## 5. Axis-by-axis results

Verdicts are for the row sets as emitted; "alleles" says whether the `REF`/`ALT` columns agree at the
positions both tools publish **from a record start**. Dense-mode fixtures publish extra materialised
rows on the GATK side — the known, out-of-scope covered-locus materialisation gap — so those are
compared at record starts only (`raw/attribution-table.txt` has the full per-field table).

### (a) `INFO/END` together with a multi-base REF

| fixture | input | GATK | native | verdict |
| --- | --- | --- | --- | --- |
| `a-end-with-concrete-alt` | `chr1 2 . AAAA AACA,<NON_REF> . PASS DP=20;END=5` | `chr1 2 . AAA AAC 92.64 . AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63` | same **plus `;END=5`** | **DIFFER (INFO/END)** |
| `end-no-trim-control` | `AAAA AACC,<NON_REF>` + `END=5` (no trim) | `…AAAA AACC…` | same **plus `;END=5`** | DIFFER (INFO/END) |
| `end-trim-consistent` | `AAAA AACA,<NON_REF>` + `END=4` | `…AAA AAC…` | same **plus `;END=4`** | DIFFER (INFO/END) |
| `end-guard-no-trim` | `AAA A,<NON_REF>` + `END=4` | `…AAA A…` | same **plus `;END=4`** | DIFFER (INFO/END) |
| `a-end-with-concrete-alt-dense` | as row 1 + `--include-non-variant-sites` | 4 rows, `AAA AAC` at 2 then materialised 3,4,5 | 1 row, `AAA AAC;END=5` at 2 | DIFFER (END + materialisation) |
| `end-block-dense` | `ACGT <NON_REF>` `END=5`, dense | `chr1 2 . ACGT .` + 3 materialised rows | `chr1 2 . C .` | see §7 (attribution probe) |

The **alleles** are identical in every case: the trim itself behaves as GATK's. The `INFO/END`
difference is pre-existing (it is present in the two no-trim controls) — but see C3: the trim is what
makes native's `END` self-inconsistent. Also measured here: the multi-base-`REF` block with `END` does
**not** make GATK crash in dense mode (`end-block-dense` exit 0), confirming that the merger really
strips `END` before `:167` (fact (g)).

### (b) Multi-base REF with a symbolic-only ALT list, and with a concrete ALT

* `b-block-multibase-ref-default` (`chr1 2 . AAAA <NON_REF>`, `GT 0/0`): **both tools emit nothing →
  EQUAL.** Reasoned+measured mechanism: the merged VC is non-variant *and* not properly polymorphic,
  so GATK returns `null` (default mode) — and native drops it too.
* `b-block-multibase-ref-end` (block + `END`), `b-block-plus-symbolic-del`
  (`ACGT <NON_REF>,<DEL>`): **both emit nothing → EQUAL.** (Both of these, and the later
  `end-block-dense`/`x-attrib-block-ref-base` probes, use a reference that is an `ACGT` repeat, so the
  record's `REF=ACGT` at POS 2 disagrees with the reference's own bases `CGTA` there. GATK does not
  reject the mismatch; for the two fixtures that emit nothing the mismatch is inert, and for
  `x-attrib-block-ref-base` it is the point of the probe — see §7.1.)
* `b-block-multibase-ref-dense`: at record start 2 GATK emits `chr1 2 . AAAA .` (untrimmed, because
  the merged emitted list is the REF alone and `getNAlleles()<=1` returns it untouched) while native
  emits `chr1 2 . A .` → **DIFFER**, attributed in §7 to native's reference-block materialisation, not
  to the trim (proved with `x-attrib-block-ref-base`).
* Multi-base REF **with** a concrete ALT: covered by 104 wide loci + all the `(d)`/`(g)` fixtures; at
  every shared record start the alleles agree.

### (c) Spanning deletion `*` with a multi-base REF, with and without an upstream deletion

* `c-star-multibase-ref-no-upstream` (`chr1 4 . AA *,<NON_REF>` `1/1`): **both emit nothing →
  EQUAL** (a lone `*` is spurious → pruned → the record is not properly polymorphic).
* `c-star-multibase-ref-with-upstream` (concrete deletion `AAAA>AA` at 2, then `AA *` at 4):
  **DIFFER — GATK 1 row, native 2 rows**:

```
GATK   : chr1 2 . AAA A 92.60 . AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63  GT:AD:DP:GQ:PL 0/1:0,20:20:99:100,0,100
NATIVE : chr1 2 . AAA A 92.60 . (identical row)
         chr1 4 . A   * 0     LowQual AC=2;AF=1.00;AN=2;DP=20;ExcessHet=0.0000;MLEAC=2;MLEAF=1.00;QD=0.00 GT:AD:DP:GQ:PL 1/1:0,20:20:99:100,100,0
```

The extra native row at 4 is the trimmed `AA→A` `*` record. GATK drops that record in default mode
because `isProperlyPolymorphic([REF,*])` is false; the `*` ownership decision agrees with native
(it is *kept* by both), only the emission differs. The same input under `--stream-by-locus`
(`h-stream-star-locus-group`) gives the identical 1-vs-2 split, so the streaming path shares this
behaviour too. Attributed in §7 to the emission rule, not to the clip.

### (d) Two-sample / three-sample merged locus with a surviving one-base ALT

* `d-three-sample-one-base-alt-pruned` (3 samples; the one-base ALT's evidence is hom-ref so it is
  pruned before `:1458`): **EQUAL** — GATK `chr1 2 . AAA AAC …` with three sample columns, native
  identical, so the guard really is consulted on the post-subset list in native.
* `d-two-sample-one-base-alt-survives` and `d2-merged-two-het-default` (two same-POS records, each
  het with strong PLs): **DIFFER — GATK emits 2 rows, native 1:**

```
GATK   : chr1 2 . AAA  AAC 190.50 … GT:AD:DP:GQ:PL 0/1:…  0/1:…     <- record 1, trimmed
         chr1 2 . AAAA A   190.46 … GT:AD:DP:GQ:PL 0/1:…  0/1:…     <- record 2, guard: UNTRIMMED
NATIVE : chr1 2 . AAA  AAC 190.50 … (byte-identical to GATK's first row, both samples)
```

GATK's default traversal regenotypes each input record separately (fact (c)), so the one-base ALT of
record 2 survives *within its own record* and `:1458` protects it; native merges the locus, the
one-base ALT loses the merged subset, and the merged allele list is trimmed. The guard of
`apply_gatk_reverse_trim()` is therefore faithful to the list it is given, but native's list is a
merged list no GATK default-mode path ever forms.
* `d2-merged-two-het-dense` (same input, dense): **GATK exits with a USER ERROR**
  (`Variant input contains more than one variant starting at location: chr1:2-5`, fact (d)); native
  emits the merged trimmed row.

### (e) Polyploid samples

`--sample-ploidy` is **not accepted by native** (`fastgatk-genotype-gvcf: unknown option:
--sample-ploidy`, exit 2, measured in `raw/run1.log` / `raw/run2.log`), so the requested flag cannot be
given to both tools. I therefore carried the ploidy in the GT/PL fields instead of the flag, on both
sides (`p3-nop-*`, `p4-nop-*`, `p3-nop-guard`) plus a wider `Number=G` diploid case
(`e-wide-diploid-3alt`):

```
p3-nop-triploid-gt  chr1 2 . AAA AAC 94.02 . AC=1;AF=0.333;AN=3;DP=20;MLEAC=1;MLEAF=0.333;QD=4.70 GT:AD:DP:GQ:PL 0/0/1:0,20:20:99:100,0,100,100      EQUAL
p4-nop-tetraploid   chr1 2 . AAA AAC 94.92 . AC=1;AF=0.250;AN=4;DP=20;MLEAC=1;MLEAF=0.250;QD=4.75 GT:AD:DP:GQ:PL 0/0/0/1:0,20:20:99:100,0,100,100,100  EQUAL
p3-nop-guard        chr1 2 . AAAA A 93.98 . AC=1;AF=0.333;AN=3;DP=20;MLEAC=1;MLEAF=0.333;QD=4.70 GT:AD:DP:GQ:PL 0/0/1:0,20:20:99:100,0,100,100      EQUAL (guard, untrimmed)
e-wide-diploid-3alt chr1 2 . AAA AAC … GT:AD:DP:GQ:PL 0/1:0,20:20:99:100,0,100                EQUAL (3rd ALT pruned by both, PL remapped to width 3)
```

So the trim behaves identically with 3- and 4-copy genotypes and the wider `Number=G` vectors are
untouched on both sides. **Not falsified** — with the caveat that GATK's own `--sample-ploidy` code
path was not exercised on the native side, because native has no such option.

### (f) Dense mode

All dense fixtures were compared **only at positions that are record starts** (explicitly out of
scope: native materialises far fewer covered coordinates than GATK — e.g. `g-star-chain-dense` GATK
positions `[2,3,4,5,6,7,8,9,10]` vs native `[2,4,7]`; the trigger fixture's GATK rows at 3 and 5
exist and native has none). At record starts, alleles agreed everywhere except the cases reported
below. Native's materialised-row *shape* differs too (`GT:AD:DP:RGQ` vs `GT:AD`, missing `INFO/DP`),
which is also part of the same out-of-scope gap.

### (g) The placement constraint (highest-value target)

Design. GATK records each emitted allele's deletion span at `GenotypingEngine.java:178-179`
(`recordDeletions`) using `vc.getReference().length() - allele.length()` and
`SimpleInterval(start, start+deletionSize)`, i.e. **before** `GenotypeGVCFsEngine.java:167` trims.
For two concrete alleles the reverse trim removes the same number of trailing bases from both, so the
recorded size is invariant and the placement is unobservable (reasoned; consistent with
`g-placement-deletion-size`, `g-placement-long-ref-star-4` where the alleles agree). The placement
*is* observable when an emitted span is later interpreted as "covering" a downstream `*` whose REF
is clipped to one base (span would collapse to nothing). Fixture `g-star-chain-dense`
(`raw/run2.log`), dense mode:

```
input:
 chr1 2 . AAAAAA AA,<NON_REF>   DP=20  GT:DP:AD:PL 0/1:20:0,20,0:100,0,100,100,100,100   (deletion size 4 -> span 2..6)
 chr1 4 . AAAAA  *,<NON_REF>    DP=20  GT:DP:AD:PL 1/1:20:0,20:100,100,0
 chr1 7 . AAAA   *,<NON_REF>    DP=20  GT:DP:AD:PL 1/1:20:0,20:100,100,0

compared at [2,4,7]:
GATK   chr1 2 . AAAAA A  92.60 . AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63  GT:AD:DP:GQ:PL 0/1:0,20:20:99:100,0,100
GATK   chr1 4 . A     *  0     LowQual AC=2;AF=1.00;AN=2;DP=20;…;QD=0.00  GT:AD:DP:GQ:PL 1/1:0,20:20:99:100,100,0
GATK   chr1 7 . A     *  0     LowQual AC=2;AF=1.00;AN=2;DP=20;…;QD=0.00  GT:AD:DP:GQ:PL 1/1:0,20:20:99:100,100,0
NATIVE  (three rows byte-identical to the three above)
# VERDICT: EQUAL
```

The `*` at 7 existing at all requires the locus at 4 to have contributed a span reaching 7 — which is
only possible if its contribution used the **untrimmed** REF (`AAAAAA`, span `[4,9]`); had native
trimmed first, the merged REF there would be one base and the contribution would be span `[4,4]`, and
the `*` at 7 would have been pruned as an orphan (reasoned counterfactual — I may not modify the
production file to measure it).

Controls (measured, `raw/run4.log`):
* `g-star-nomid-control` (same fixture without the record at 4): the only recorded span is `[2,6]`,
  which does not reach 7 → **both** tools prune the `*` at 7 and emit `chr1 7 . AAAA .` (alleles
  equal; only QUAL differs, §7). So the coverage test is live and both tools agree on the negative
  case, which is what makes the positive case informative.
* `g-star-mid-only-control` (records at 2 and 4): **EQUAL**, `chr1 4 . A *`.
* `h-stream-placement-chain`: the same chain through `--stream-by-locus`: rows at `[2,4,7]` agree.

**Not falsified**, in the aggregate path or the streaming path.

### (h) The second traversal

`--help` lists `--stream-by-locus` ("bounded k-way joint-locus merge"); by inspection both compute
stages route through `apply_gatk_output_allele_subset()` (`genotype_gvcf_tool.cpp:6349` streaming,
`:7069` aggregate), which is where the helper is called. By execution:

| fixture (same input, same flags) | aggregate path | `--stream-by-locus` |
| --- | --- | --- |
| suffix substitution `AAAA AACA` (default mode) | **EQUAL**, `chr1 2 . AAA AAC …` both | `h-stream-suffix-substitution` **EQUAL**, same row |
| upstream deletion + `*` locus group (default) | `c-star-multibase-ref-with-upstream` DIFFER (native-only row 4) | `h-stream-star-locus-group`, identical 1-vs-2 split |
| dense single locus | `h-stream-dense` alleles EQUAL (GATK-only materialised 3,4,5) | same |
| `*` chain, placement (dense) | `g-star-chain-dense` EQUAL at `[2,4,7]` | `h-stream-placement-chain` EQUAL at `[2,4,7]` |
| `INFO/END` (default) | `a-end-with-concrete-alt` DIFFER (END) | `h-stream-end` DIFFER (END) |

("EQUAL"/"DIFFER" for the dense rows are at record-start positions only, as everywhere in this report;
the GATK-only positions are the known materialisation gap.)

So the fix is present in both paths (verified by execution, not only inspection) and the paths behave
the same on my fixtures.

---

## 6. Counterexamples

### C1 (strongest, trim-specific): `INFO/DP=0` — native trims where GATK provably does not

Input (one record, no exotic alleles):

```
chr1 2 . AAAA AACA,<NON_REF> . PASS DP=0 GT:DP:AD:PL 0/1:20:0,20,0:100,0,100,100,100,100
```

Dense mode (`dp0-dense`, `raw/run2.log`; identical result in `r-info-dp-zero-dense`, `raw/run1.log`):

```
GATK   : chr1 2 . AAAA AACA . . DP=0 GT:AD ./.:0,20          <- untrimmed (the `result = originalVC` arm)
NATIVE : chr1 2 . AAA  AAC 92.64 . AC=1;AF=0.500;AN=2;DP=0;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63 GT:AD:DP:GQ:PL 0/1:0,20:20:99:100,0,100
# VERDICT: DIFFER at position 2, a record-start position both tools publish
```

Default mode (`r-info-dp-zero-default`): **GATK emits 0 rows**, native emits the trimmed row
`chr1 2 . AAA AAC 92.64 . …;DP=0;… GT:… 0/1:…`.

Measured boundary: `INFO/DP=1` + dense (`dp1-dense`) → both `AAA AAC`, EQUAL; `INFO/DP` absent from
the input (`r-no-info-dp`) → the merged depth falls back to the sum of genotype depths (20), both tools
trim, and the two rows differ only in that native omits the `INFO/DP` key (§7.6). So the trim
divergence is exactly at `DP == 0`.

Mechanism (reasoned, from the sources cited in §2(e)): `GenotypeGVCFsEngine.regenotypeVC()` gates the
whole regenotyping block — including `:167` — on `getAttributeAsInt(DP,0) > 0`; the merged VC's
`INFO/DP` is the input record's own `INFO/DP` when present
(`ReferenceConfidenceVariantContextMerger.calculateVCDepth`). Native's
`apply_gatk_output_allele_subset()` has no such precondition and calls
`apply_gatk_reverse_trim()` unconditionally, so it trims a record GATK leaves alone (dense) or
discards (default).

**Answer to "does native ever trim where GATK does not?": yes.**

### C2: the same gate removes GATK's deletion bookkeeping — `*` ownership diverges

Input (dense, `dp0-star-ownership-dense`, `raw/run4.log`):

```
chr1 2 . AAAA AA,<NON_REF> . PASS DP=0  GT:DP:AD:PL 0/1:20:0,20,0:100,0,100,100,100,100
chr1 4 . AA   *,<NON_REF>  . PASS DP=20 GT:DP:AD:PL 1/1:20:0,20:100,100,0

compared at [2,4]:
GATK   : chr1 2 . AAAA AA .  .  DP=0                GT:AD          ./.:0,20
GATK   : chr1 4 . AA   .  163.06 . DP=20;MLEAC=.;MLEAF=. GT         ./.          <- '*' pruned as an orphan
NATIVE : chr1 2 . AAA  A  92.60 . AC=1;AF=0.500;AN=2;DP=0;…;QD=4.63  GT:AD:DP:GQ:PL 0/1:0,20:20:99:100,0,100
NATIVE : chr1 4 . A    *  0     LowQual AC=2;AF=1.00;AN=2;DP=20;…;QD=0.00 GT:AD:DP:GQ:PL 1/1:0,20:20:99:100,100,0
# VERDICT: DIFFER at both positions
```

Two independent divergences in one fixture, both traceable to the missing `DP>0` precondition: the
locus at 2 is trimmed by native and not by GATK, **and** because GATK never calls
`recordDeletions()` for that locus, GATK's `*` at 4 has no covering deletion and is pruned
(`AA .`), whereas native — which did record the deletion — keeps it and trims the REF (`A *`).
Note that the `*` ownership arithmetic itself (`loc.start < vc.start <= loc.end`) is identical in both
implementations (reasoned: `GenotypingEngine.java:365-372` vs `EmittedDeletions::covers`,
`genotype_gvcf_tool.cpp:1281-1285`), as every other `*` fixture shows.

### C3: the trim leaves native's pass-through `INFO/END` stale

`a-end-with-concrete-alt` (identical commands, identical alleles, one INFO field different):

```
input : chr1 2 . AAAA AACA,<NON_REF> . PASS DP=20;END=5  GT:DP:AD:PL 0/1:20:0,20,0:100,0,100,100,100,100
GATK  : chr1 2 . AAA AAC 92.64 . AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63
NATIVE: chr1 2 . AAA AAC 92.64 . AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63;END=5
```

Attribution, from two more measurements:
* The *presence* of `END` is pre-existing: `end-no-trim-control` (REF/ALT share no trailing base, so
  GATK returns the record untouched) shows native emitting `…;END=5` on a record whose span really is
  2..5, i.e. the difference exists with the trim disabled.
* The *staleness* is caused by the trim: after native's (GATK-correct) `AAAA→AAA` rewrite the record
  spans 2..4 while it still advertises `END=5` (rlen 3 vs advertised end 5). htsjdk's own
  `VariantContextBuilder` rejects that combination (`Badly formed variant context at location
  chr1:2; …` — measured at the unit level in `raw/p2.out`).
* Scope of the impact, measured: htsjdk's reader accepts native's file and re-emits the `END`
  (`SelectVariants` run, exit 0 on the row), and `gatk ValidateVariants -V <native output> -R <ref>`
  passes with no complaint. So this is a byte-level divergence plus a misleading coordinate, **not** a
  hard invalidity; I did not find a GATK tool that rejects the file.

---

## 7. Differences measured but attributed elsewhere (explicitly **not** trim findings)

1. **Dense-mode block REF comes from the FASTA base, not the trim.** `x-attrib-block-ref-base` is a
   deliberate attribution probe whose input REF disagrees with the reference on purpose: reference
   `ACGTACGT…`, record `chr1 2 . ACGT <NON_REF>`. GATK emits `chr1 2 . ACGT .`; native emits
   `chr1 2 . C .`, i.e. the FASTA base at POS 2 — whereas the reverse trim of `ACGT` would give `A`
   and "no rewrite" would give `ACGT`. Native's dense-block REF thus comes from
   `reference_base()` in the reference-block materialisation — aggregate path
   `split_reference_blocks_at_variants()` (`genotype_gvcf_tool.cpp:3636-3664`, the rewrite is
   `alleles = reference + ",<NON_REF>"`, `alleles.front() = reference`), streaming path
   `emit_lazy_stream_reference_piece()`/`split_stream_reference_block()` (`:5570`, `:5598`, `:5622`) —
   and never from `apply_gatk_reverse_trim()`. This also explains
   `b-block-multibase-ref-dense` (`AAAA .` vs `A .`), where both mechanisms coincide byte-for-byte.
   (The same probe doubles as the `END`-with-block case; note that its input REF is invalid on
   purpose — GATK does not reject it.)
2. **Default-mode traversal.** Two records at one POS: GATK emits one row per input record
   (fact (c)); native merges them into one locus. Dense mode: GATK refuses the input outright
   (fact (d)). Both measured (§5(d), §5(c)); the native merged row is byte-identical to GATK's
   equivalent row, so the trim is not what differs.
3. **Row-emission rule for `*` records in default mode.** GATK drops a record whose regenotyped
   allele list is `[REF,*]` (`isProperlyPolymorphic` false); native emits it (trimmed, `A *`). The
   `*` ownership decision itself agrees (measured in §5(c), §5(g)).
4. **Monomorphic REF-only QUAL in dense mode.** Both directions measured, with identical alleles:
   `star-only-nonref-dense` GATK `163.06` vs native `Infinity`; `g-longref-phantom` at 4 GATK
   `Infinity` vs native `159.55`; `g-star-nomid-control` at 7 GATK `163.06` vs native `Infinity`. In
   each of these the emitted list has one allele, so no trim ran inside either tool.
5. **`INFO/QD` at loci with two concrete ALTs**: 16/44 wide stage-B rows show `QD=4.15` (GATK) vs
   `4.16` (native) with identical alleles, QUAL, AC/AF/AN/MLEAC/MLEAF and FORMAT columns.
6. **`INFO/DP` omission**: `r-no-info-dp` (no `INFO/DP` in the input) — GATK writes `DP=20`, native
   omits the key; alleles identical.
7. **Dense-mode row shape / materialisation gap** (out of scope by instruction): GATK emits a row for
   every covered coordinate with `GT:AD` + `INFO/DP`, native emits only record starts with
   `GT:AD:DP:RGQ` and no `INFO/DP`. Position lists are printed for every dense fixture in
   `raw/run*.log`.
8. **`--sample-ploidy` is not implemented by native** (`unknown option`, exit 2). A capability gap,
   and the reason the polyploid axis was exercised through the GT/PL fields instead.

---

## 8. The five questions

**(1) Does native now trim in EVERY case where GATK trims, on your fixtures?**
**Yes on every fixture where the comparison is possible** — in each of the 104 wide-sweep loci and in
every targeted fixture in which GATK emitted a trimmed record at a position both tools published from
a record start, native emitted the same trimmed alleles (`AAA/AAC`, `ACGTA/A`, `AA/A`, `TCG/TGG,TCT`,
`GAT/…`, …). I could not construct a case where GATK trims and native does not; the observed
asymmetry is always the other way (C1/C2). Caveat: this is evidence, not a proof — native can only be
shown to trim *more*, never *less*, when its emitted allele list is a superset or equal.

**(2) Does native ever trim where GATK does NOT?**
**Yes — measured, three independent fixtures.** `dp0-dense` and `r-info-dp-zero-dense` (dense: GATK
`AAAA AACA` untrimmed at POS 2 vs native `AAA AAC`), `r-info-dp-zero-default` (GATK 0 rows vs native
1 trimmed row), and `dp0-star-ownership-dense` (both the trim at 2 and the `*` decision at 4 differ).
Mechanism: GATK gates the whole regenotyping-and-trim block on the merged `INFO/DP > 0`
(`GenotypeGVCFsEngine.java:152-175`), which native's helper call does not reproduce. Byte-level
consequence also measured in C3 (stale `END` after a trim GATK would not have made on an
`END`-bearing record — GATK strips `END` so the question is moot for GATK, but native's record is
self-inconsistent afterwards).

**(3) Is the clip length (including the clip-one-fewer emptiness rule) exactly GATK's on your
fixtures?**
**Yes.** 60 + 44 = 104 loci covering common trailing runs of 0-3, deletions of 1-3 bases, one-base
ALTs, insertions and two-concrete-ALT combinations: zero allele differences and zero disagreements
with an independent Python transcription of the helper (`raw/wideA.log`, `raw/wideB.log`). Emptiness
corners measured on both sides: `ACGTACGT/ACGT → ACGTA/A` (clip 3, not 4) and 12A/11A → `AA/A`
(clip 10, not 11) in `g-longref-default`/`g-longref-phantom`; a single-candidate list clips to one
base (`AAAA/<NON_REF> → A`). The unit-level oracle shows the same rule in the shipped Java.

**(4) Is POS ever moved by native's trim?**
**No.** No fixture produced a native-only position that would indicate a moved record start; every
divergence is either an extra/missing row or a difference in the same row. Every trimmed native row I
inspected kept GATK's POS with a shortened REF (`chr1 2 . AAA AAC`, `chr1 4 . A *`,
`chr1 2 . ACGTA A`-style cases in the gate). This is measured (position sets printed per fixture) and
consistent with the helper writing no `pos` and only calling `bcf_update_alleles_str()`
(`genotype_gvcf_tool.cpp:3083-3085`, reasoned).

**(5) Are the FORMAT fields (GT/AD/DP/GQ/PL) left untouched by the trim in native, as GATK's identity
allele remap implies?**
**Yes, compared to GATK, on every fixture where the alleles match.** The per-field attribution table
(`raw/attribution-table.txt`) lists `FORMAT`/`SAMPLE` differences for exactly one fixture,
`end-block-dense`, and only at positions 3, 4, 5 — GATK's materialised rows, whose FORMAT shape
(`GT:AD`) differs from native's record-start shape (`GT:AD:DP:RGQ`); at its record start (POS 2) the
FORMAT column is `GT:AD` vs `GT:AD:DP:RGQ`, again the out-of-scope dense shape gap, with no trim
involved. Concrete evidence for the trim itself: `p3-nop-*`/`p4-nop-*`/`e-wide-diploid-3alt` keep
identical GT (`0/0/1`, `0/0/0/1`), AD and PL vectors (widths 4, 5, 3) after the trim, and
`dp1-dense`, all wide loci, and the gate's fixtures have byte-identical FORMAT columns. Code-read
supports this: the helper only rewrites alleles (`record.alleles`, `bcf_update_alleles_str`) and never
touches GT/GQ/AD/DP/PL, and GATK's trim remaps genotypes by allele identity
(`updateGenotypesWithMappedAlleles`, `:1506-1512`).

---

## 9. Claims I could not falsify (negative results, with the command that shows them)

* **Clip algorithm** — `python3 wide.py --stage A` (60 loci) and `--stage B` (44 loci): `0 not EQUAL`
  / `0 allele differences`. Every attempt to make the run/cap/emptiness rule diverge failed.
* **Placement (`:178-179` before `:167`)** — `python3 probe2.py g-star-chain-dense` → EQUAL at
  `[2,4,7]`; negative control `python3 probe4.py g-star-nomid-control` → both tools prune the `*` at 7
  (alleles equal); streaming twin `probe2.py h-stream-placement-chain` → EQUAL.
* **Guard on the post-subset list** — `python3 probe2.py prune-refonly-dense` EQUAL (`AAAA .`), and
  `probe.py d-three-sample-one-base-alt-pruned` EQUAL.
* **Both traversal paths carry the fix** — `h-stream-suffix-substitution` EQUAL, `h-stream-dense`
  alleles EQUAL, `h-stream-placement-chain` EQUAL, plus the call-site inspection at `:6349` / `:7069`.
* **The registered gate itself** — strict, both backends: `status: pass`, `violations: 0`, exit 0
  (`raw/gate-strict-build.log`, `raw/gate-strict-build-serial.log`), and all seven cases also pass in
  `--expect-divergence` mode (`raw/gate-divergence-parallel.log`), i.e. the task's baseline claim is
  reproduced independently.

---

## 10. Limits of this verification (what I did not test, and why)

1. **`--sample-ploidy`** cannot be passed to native at all (measured exit 2), so GATK's own
   ploidy-flag code path was not compared; I substituted GT/PL-carried ploidy (all EQUAL).
2. **Dense-mode comparisons are record-start only** by instruction: native does not materialise every
   covered coordinate, so full row-set equality is not achievable and is out of scope. Where a dense
   fixture is reported as DIFFER I state whether the difference is at a shared record start.
3. **Symbolic non-`<NON_REF>` ALTs** (e.g. `<DEL>`) were probed only at the unit level
   (`symbolic-only`, `symbolic-plus-nonref`: REF clipped to one base, symbolic copied) — GATK's gVCF
   input contract demands `<NON_REF>`, so I did not build an end-to-end fixture for them.
4. **Not exercised**: `-L`/`-XL` intervals, multiple `-V` inputs, `--max-alternate-alleles`,
   `--genotype-assignment-method` variants, `--annotate-with-num-discovered-alleles`, `.gz` inputs and
   `--only-output-calls-starting-in-intervals` — all orthogonal to the trim by inspection, none tested.
5. **Single-run evidence.** Every fixture was executed once per tool; I did not repeat runs to bound
   run-to-run nondeterminism (native's QUAL/QD paths consume a deterministic RNG stream, and repeated
   identical outputs were observed incidentally for duplicated fixtures, e.g. `p4-trim` vs
   `p4-nop-tetraploid-gt`).
6. **The trim-first counterfactual in the placement test is reasoned, not measured.** I may not modify
   the production helper, so "the `*` at 7 would have been pruned had the trim run first" is an
   argument from the merged REF length plus the measured span arithmetic, not an execution.
7. **`INFO/DP=0` is a contrived, though well-formed, input shape.** I did not find a realistic gVCF
   producer that writes `INFO/DP=0` on a called record; the finding is about a missing precondition
   rather than about a common corpus. The reach argument is the mechanism, and it is source-verified.

---

## 11. Reproduction

Fixtures and drivers (raw copies in `fastgatk-native/evidence/2026-09-13-round49/raw/`, runnable
copies in `.diag/round49/`):

```
.diag/round49/probe.py   # 24 fixtures: axes a,b,c,d,e,f,g,h + reach controls (batch 1)
.diag/round49/probe2.py  # 19 fixtures: END, prune-to-REF-only, DP boundary, '*' chain, streaming (batch 2)
.diag/round49/probe3.py  #  7 fixtures: polyploid without the flag, merged-locus dense, attribution probe (batch 3)
.diag/round49/probe4.py  #  3 fixtures: placement controls (batch 4)
.diag/round49/wide.py    # 60 + 44 shapes in one GATK run each               (wide sweep)
.diag/round49/TrimProbe.java, P2.java   # unit-level oracle against the pinned jar
.diag/round49/adv_compare.py, adv_analyze.py   # row/field attribution over work/
.diag/round49/work/, .diag/round49/wide/       # every input, reference, output and command log
```

Commands:

```bash
cd /home/turing-agents/Documents/fast-gatk/.diag/round49
python3 -u probe.py            > run1.log 2>&1     # batch 1 (24 fixtures)
python3 -u probe2.py           > run2.log 2>&1     # batch 2
python3 -u probe3.py           > run3.log 2>&1     # batch 3
python3 -u probe4.py           > run4.log 2>&1     # batch 4
python3 -u wide.py --stage A   > wideA.log 2>&1
python3 -u wide.py --stage B   > wideB.log 2>&1
python3 adv_analyze.py         > attribution-table.txt
third_party/jdk17/bin/javac -cp <gatk-jar> -d . TrimProbe.java P2.java
third_party/jdk17/bin/java  -Xmx1g -cp <gatk-jar>:. TrimProbe > trimprobe.out
third_party/jdk17/bin/java  -Xmx1g -cp <gatk-jar>:. P2        > p2.out
# registered gate, both backends
python3 fastgatk-native/scripts/verify_genotype_gvcf_reverse_trim_gatk_oracle.py --expect-divergence
FASTGATK_NATIVE_BUILD=$PWD/fastgatk-native/build        python3 fastgatk-native/scripts/verify_genotype_gvcf_reverse_trim_gatk_oracle.py
FASTGATK_NATIVE_BUILD=$PWD/fastgatk-native/build-serial python3 fastgatk-native/scripts/verify_genotype_gvcf_reverse_trim_gatk_oracle.py
```

Raw rows for every fixture, with the exact `gatk`/`native` command line echoed above each row set:
`raw/all-rows.txt`; per-field attribution: `raw/attribution-table.txt`.

### Files written by this round

* `fastgatk-native/evidence/2026-09-13-round49/reverse-trim-adversarial-verification.md` (this report)
* `fastgatk-native/evidence/2026-09-13-round49/raw/` — `run1.log`…`run4.log`, `wideA.log`,
  `wideB.log`, `trimprobe.out`, `p2.out`, `all-rows.txt`, `attribution-table.txt`,
  `gate-*.log`, and copies of the four probe drivers, `wide.py`, `TrimProbe.java`, `P2.java`
* scratch: `.diag/round49/{probe,probe2,probe3,probe4,wide,adv_compare,adv_analyze}.py`,
  `.diag/round49/{work,wide}/` (inputs, references, outputs, manifests)

No production file, `CMakeLists.txt` or `third_party/` file was modified; no `git add`/`git commit`
was run; no binary was rebuilt.
