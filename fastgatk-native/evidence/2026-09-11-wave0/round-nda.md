# Round: the `--annotate-with-num-discovered-alleles` data-row contract (INFO/NDA)

Scope: the divergence the previous round **measured but deliberately left**
(`.diag/round-header-lines.md` §7, last bullet): with
`--annotate-with-num-discovered-alleles`, pinned GATK 4.6.2.0 writes `INFO/NDA`
on the output row while native wrote **no `NDA` key at all**.  That is why the
previous round could fix only the NDA *header* wording and could not gate the
flag.  Gate first, fix second.

## 0. Verdict

* **GATK truth (measured, pinned 4.6.2.0):** `NDA` is the number of ALT alleles
  of the **merged input record**, with `<NON_REF>` excluded and a symbolic `*`
  **included** — not the number of *published* alleles, and not the number of
  per-sample discoveries.  It therefore reads **1** on every single-concrete-ALT
  locus, which is exactly the value native omitted.  The key and its
  `##INFO` declaration appear **only** under the flag.  A locus whose merged ALT
  set is empty (a plain `<NON_REF>`-only block forced out by
  `--include-non-variant-sites`) carries **no `NDA`** on either side.
* **GATK rule with source:** `GenotypingEngine.composeCallAttributes()` stores
  `vc.getAlternateAlleles().size()` of the VariantContext handed to
  `calculateGenotypes()` (`GenotypingEngine.java:464-465`).  That `vc` is the
  *pre-pruning* record: `reducedVC` for `--max-alternate-alleles` is built at
  `:137-144` and `outputAlternativeAlleles` at `:155`, both **after** the value
  is read.  `<NON_REF>` is absent from `vc` because
  `GenotypeGVCFsEngine.java:136` calls
  `merge(..., removeNonRefSymbolicAllele = true, ...)` and
  `ReferenceConfidenceVariantContextMerger.java:343-345` only re-adds it when
  that flag is false, while `:340-342` re-adds `Allele.SPAN_DEL` for a spanning
  event.  The declaration is added behind the same boolean through
  `GenotypingEngine.getAppropriateVCFInfoHeaders()` (`:90-94`) from
  `GenotypeGVCFsEngine.setupVCFWriter()` (`:401-402`), text
  `GATKVCFHeaderLines.java:205`.
* **Fix kept in the tree:** one condition in one function in
  `fastgatk-native/src/genotype_gvcf_tool.cpp` — `record.alleles.size() < 3`
  became `< 2` in `annotate_num_discovered_alleles()` (1 functional token; the
  rest of the `+29/-6` is the provenance comment).  No restructuring, no kernel
  ABI, no `CMakeLists.txt`, no Mutect2, no other tool.
* **Gate:** `verify_genotype_gvcf_spandel_gatk_oracle.py` extended with six
  gated `nda-annotation-*` cases.  **exit 1 before** the fix (5 violations, one
  per failing row, all of them the missing `NDA` key — §2.1), **exit 0 after**
  (37 cases, 0 violations).
* Step 5: (a) exit 0; (b1) exit 0; (b2) 15/15 exit 0; (c) 19/19 exit 0;
  (d) 300/300 on **both** backends, exit 0.  Change **kept in the tree**.

## 1. STEP 1 — measured GATK truth per case

Reproducers: `.diag/nda_probe.py` (`.diag/nda-probe.log`),
`.diag/nda_probe2.py` (`.diag/nda-probe2.log`) and
`.diag/nda_gate_probe.py` (`.diag/nda-gate-probe-before.log`, which runs **both**
tools).  Every run is `GenotypeGVCFs -R <100 bp chr1> -V <fixture>` with and
without `--annotate-with-num-discovered-alleles`; native additionally gets
`--gatk-compatible-annotations`.

### 1.1 Which rows carry NDA, and with what value

| fixture (input ALT) | published ALT | `NDA` with flag | `NDA` without flag |
| --- | --- | --- | --- |
| `G,<NON_REF>` called 0/1 | `G` | **1** | absent |
| `G,T,<NON_REF>` called 1/2 | `T` only (G pruned) | **2** | absent |
| `*,G,<NON_REF>` called 0/2 | `G` (`*` pruned) | **2** | absent |
| `*,<NON_REF>` + `--include-non-variant-sites` | `.` (REF-only) | **1** | absent |
| `G,<NON_REF>`, two samples | `G` | **1** | absent |
| `G,T,<NON_REF>` called 1/2 | `G,T` | **2** | absent |
| `<NON_REF>` only + `--include-non-variant-sites` | `.` (REF-only) | **absent** | absent |
| `<NON_REF>` only, default mode | (no record) | — | — |
| two samples, `G` vs `T` at one locus | two rows `G`, `T` | **1** each | absent |

Literal rows (pinned GATK 4.6.2.0, `.diag/nda-gate-probe-before.log`):

```
# single concrete ALT, with the flag
chr1  2  .  A   G   92.64  .  AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;NDA=1;QD=4.63  GT:AD:DP:GQ:PL  0/1:0,20:20:99:100,0,100
# the same fixture, without the flag
chr1  2  .  A   G   92.64  .  AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63         GT:AD:DP:GQ:PL  0/1:0,20:20:99:100,0,100
# '*',G pruned to G: the count is the INPUT count
chr1  2  .  A   G   82.26  .  AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;NDA=2;QD=4.11  GT:AD:DP:PL     0/1:0,20:20:0,0,0
# REF-only dense row materialized from a covered '*': the row has ALT='.', GT only, and still NDA=1
chr1  3  .  A   .   Infinity  .  DP=20;MLEAC=.;MLEAF=.;NDA=1                                             GT              ./.
# a locus whose merged ALT set is EMPTY: no NDA at all (and no MLEAC/MLEAF either)
chr1  2  .  A   .   .         .  DP=20                                                                   GT:AD           ./.:20
```

### 1.2 Does pruning change the count?  Does it come from INPUT or PUBLISHED?

**Input.**  `*,G,<NON_REF>` publishes only `G` yet reports `NDA=2`; `G,T,<NON_REF>`
with `T` pruned reports `NDA=2`.  The source explains it: the value is read from
`vc` (`GenotypingEngine.java:464-465`), the argument of `calculateGenotypes()`,
while the published subset is `outputAlternativeAlleles`, computed 300 lines
earlier at `:155`.

Max-ALT reduction is the third variant of the same question, and it is already
measured by a **registered** gate:
`verify_genotype_gvcf_max_alternate_alleles_gatk_oracle.py:77-94` feeds 7
concrete ALTs with `--max-alternate-alleles 6` and asserts GATK writes
`ALT=C,G,T,AA,AC,AG` (6 published) **with `NDA=7`** — input count, not reduced
count.  That test still passes after this round's change (exit 0, §5), and it is
the regression evidence that the `< 2` edit did not disturb the multi-ALT path.

### 1.3 Only with the flag?

Yes, in both directions: without the flag neither the key nor
`##INFO=<ID=NDA,...>` appears (`compare [with-flag]`/`[without-flag]` columns of
§1.1 and `nda_header` in `.diag/nda-probe.log`); with the flag GATK emits the
key on every emitted row whose merged ALT set is non-empty.  Both halves are
behind the same boolean (`GenotypingEngine.java:464` for the key, `:90-94` for
the declaration).  The gate asserts this per case, not just via the row diff
(§2.2).

### 1.4 Where GATK computes it — quotes

```
GenotypeCalculationArgumentCollection.java:73   @Argument(fullName = "annotate-with-num-discovered-alleles", ...)
GenotypeCalculationArgumentCollection.java:74   public boolean ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED = false;

GenotypingEngine.java:464-465   if ( configuration.genotypeArgs.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED ) {
                                    attributes.put(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY,
                                                   vc.getAlternateAlleles().size());

GenotypingEngine.java:90-94     public Set<VCFInfoHeaderLine> getAppropriateVCFInfoHeaders() {
                                    ...
                                    if ( configuration.genotypeArgs.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED ) {
                                        headerInfo.add(GATKVCFHeaderLines.getInfoLine(
                                                GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY));

GATKVCFConstants.java:73        public static final String NUMBER_OF_DISCOVERED_ALLELES_KEY = "NDA";
GATKVCFHeaderLines.java:205     addInfoLine(new VCFInfoHeaderLine(NUMBER_OF_DISCOVERED_ALLELES_KEY, 1,
                                    VCFHeaderLineType.Integer, "Number of alternate alleles discovered
                                    (but not necessarily genotyped) at this site"));

GenotypeGVCFsEngine.java:136    merger.merge(variantsToProcess, loc, ref.getBase(), true, false);
GenotypeGVCFsEngine.java:401-402  headerLines.addAll(annotationEngine.getVCFAnnotationDescriptions(false));
                                  headerLines.addAll(genotypingEngine.getAppropriateVCFInfoHeaders());
GenotypeGVCFsEngine.java:233-250  addGenotypingAnnotations(): the NDA carry is at :238-239, in a
                                  LinkedHashMap whose insertion order is MLEAC, MLEAF, NDA, AS_QUAL
ReferenceConfidenceVariantContextMerger.java:340-342  if (sawSpanningDeletion && (...)) add(Allele.SPAN_DEL);
ReferenceConfidenceVariantContextMerger.java:343-345  if (!removeNonRefSymbolicAllele) add(Allele.NON_REF_ALLELE);
GenotypeGVCFsEngine.java:148/:154  regenotypeVC()'s `originalVC.isVariant()` test: the empty-ALT-set
                                  REF-only row is built on the other branch and never reaches
                                  composeCallAttributes()
```

The `NDA` key's position in the row (between `MLEAF` and `QD`) follows from the
`LinkedHashMap` order at `:233-250`; native already produced that order, which is
why the measured rows differ only by the missing key.

## 2. STEP 2 — the gate, run BEFORE the fix

**Extended** the already-registered
`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` rather than
adding a new script: it already runs pinned GATK and native on the same fixtures
with the same arguments, and it is already registered in
`fastgatk-native/CMakeLists.txt` under the step-5(c) filter (`:1556-1559`, not
edited this round).  Oracle diff `+234/-0`.

Added:

* a docstring section *"The `--annotate-with-num-discovered-alleles` data-row
  contract"* quoting the `file:line`s of §1.4;
* fixtures `NDA_SINGLE_ALT_RECORD`, `NDA_TWO_SAMPLES_RECORD` and measured row
  constants `GATK_NDA_SINGLE_ALT_ROW`, `GATK_NDA_CONTROL_ROW`,
  `GATK_NDA_TWO_SAMPLES_ROW`, `GATK_NDA_REDUCED_ALT_SET_ROW`,
  `GATK_NDA_COVERED_STAR_PLUS_CONCRETE_ROWS`,
  `GATK_NDA_COVERED_STAR_ONLY_DENSE_ROWS`, plus the literal
  `GATK_NDA_HEADER_LINE` (`GATKVCFHeaderLines.java:205`);
* six **GATED** cases (the existing row comparison is byte-exact, so the `NDA`
  key is inside the assertion by construction):
  `nda-annotation-single-alt`, `nda-annotation-control-without-flag`,
  `nda-annotation-reduced-alt-set`, `nda-annotation-covered-star-plus-concrete-alt`,
  `nda-annotation-covered-star-only-dense`, `nda-annotation-two-samples`;
* a per-case assertion that the `##INFO=<ID=NDA,...>` declaration is present iff
  `--annotate-with-num-discovered-alleles` is in the shared argument list — on
  **both** sides, so the case cannot pass because both tools dropped it.

The pre-existing `expect` check (GATK rows vs. the measured constant) makes each
new case fail if pinned GATK ever stops writing `NDA`, so the gate cannot pass
for the wrong reason.

### 2.1 Exit status and the literal rows, BEFORE the fix

**exit 1**, `"status": "divergence"`, **5 violations**, every one of them a row
diff whose only difference is the `NDA` key (`.diag/nda-gate-before.log`):

```
[nda-annotation-single-alt] VIOLATION: row 0 is not byte-identical:
  GATK  ='...MLEAC=1;MLEAF=0.500;NDA=1;QD=4.63\tGT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100'
  NATIVE='...MLEAC=1;MLEAF=0.500;QD=4.63\tGT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100'
[nda-annotation-covered-star-plus-concrete-alt] VIOLATION: row 0 ... (same diff, locus 2)
[nda-annotation-covered-star-only-dense] VIOLATION: row 0 ... (same diff, locus 2)
[nda-annotation-covered-star-only-dense] VIOLATION: row 1 is not byte-identical:
  GATK  ='chr1\t3\t.\tA\t.\tInfinity\t.\tDP=20;MLEAC=.;MLEAF=.;NDA=1\tGT\t./.'
  NATIVE='chr1\t3\t.\tA\t.\tInfinity\t.\tDP=20;MLEAC=.;MLEAF=.\tGT\t./.'
[nda-annotation-two-samples] VIOLATION: row 0 ... (same diff, two samples)
```

The two cases that do **not** violate before the fix are the control
(`nda-annotation-control-without-flag`, where both sides correctly omit the key)
and `nda-annotation-reduced-alt-set` (a 2-ALT input, which the old `>= 2 ALT`
guard already admitted) — they are carried anyway because they pin the two rules
the fix must not break.  Header content matched byte-for-byte as a multiset on
all six cases, including the flag-conditional `NDA` declaration.

## 3. STEP 3 — the fix

`git diff --stat` = **2 files changed, 257 insertions(+), 6 deletions(-)**
(oracle `+234`, tool `+29/-6`; most of the tool's lines are the provenance
comment).

`fastgatk-native/src/genotype_gvcf_tool.cpp`, `annotate_num_discovered_alleles()`:

```diff
 void annotate_num_discovered_alleles(const bcf_hdr_t* output_header,
                                      Record& record,
                                      const Options& options) {
     if (!options.annotate_with_num_discovered_alleles ||
-        record.alleles.size() < 3 ||
+        record.alleles.size() < 2 ||
         bcf_hdr_id2int(output_header, BCF_DT_ID, "NDA") < 0)
         return;
```

plus a rewritten leading comment stating the measured rule and its `file:line`
provenance (the only other change; no code line moved).

**Why that one token is the whole fix.**  The function already had the right
value and the right timing: it is the first thing the compute stage does, before
`apply_gatk_max_alternate_alleles()` and before
`apply_gatk_output_allele_subset()` (`genotype_gvcf_tool.cpp:5988` streaming,
`:6708` aggregate), i.e. exactly GATK's "read `vc`, not `reducedVC`, not
`outputAlternativeAlleles`", and its loop already skips `<NON_REF>` and counts
`*`.  `record.alleles` is the merged union with `<NON_REF>` removed, matching
GATK's merged `vc`, and `record.alleles.size()` includes the reference, so the
old guard `size() < 3` meant "fewer than two ALTs" and silently suppressed
`NDA=1`.  `size() < 2` means "no ALT at all", which is the genuinely-empty case
where GATK also writes nothing; the loop's own `if (discovered <= 0) return;`
remains as the second guard.

**Blast radius, explicitly.**

* The function returns immediately unless `options.annotate_with_num_discovered_alleles`
  is set, so **no run that does not pass the flag can change behaviour at all**.
* One caller shape only: `annotate_num_discovered_alleles()` has exactly two
  call sites, both inside `fastgatk-native/src/genotype_gvcf_tool.cpp`, and it is
  declared in no header.
* Not touched: the genotyping engine, the Kokkos kernel ABI, the FILTER logic,
  the AF/PL kernels, output-allele subsetting, the header builders (the previous
  round's NDA wording fix is untouched and still correct), the non-GATK
  diagnostic profile, `fastgatk-native/CMakeLists.txt`, every root `*.md`,
  Mutect2, every other tool.
* No stale assertion to correct (§4).

## 4. STEP 4 — stale assertions

**No registered test pinned the old behaviour; no test line was corrected.**

| File / check | Why it could have been stale | Outcome |
| --- | --- | --- |
| `verify_genotype_gvcf_max_alternate_alleles_gatk_oracle.py:94` (`assert "NDA=7" in gatk_row[7].split(";")`) | the only other place in the tree that mentions `NDA` | **not stale**: it asserts GATK's own value (8 ALTs incl. `<NON_REF>` → 7 concrete) and then requires native's rows to equal GATK's byte for byte (`:111-113`); the input has 7 ALTs, so the old guard already admitted it and the new one cannot change it.  **Re-run: exit 0** (`{"java_native_rows_exact": true, "stream_rows_exact": true}`) |
| `grep -rn 'NDA=' fastgatk-native/ --include=*.py --include=*.cpp --include=*.hpp` | any hard-coded expectation of "no NDA" / "NDA only when >= 2 ALTs" | only the line above; no such expectation exists |
| the 31 pre-existing cases of this oracle | they never pass the flag | not stale (see the blast-radius argument in §3); all 31 still pass in step 5(a) |
| `genotype_gvcf_tool.cpp` self-comment | it claimed "Pure REF/<NON_REF> reference-confidence blocks do not receive NDA" and `>= two ALTs` | **stale prose, corrected in this round** — replaced by the measured rule and its provenance.  No assertion depended on it |

The only test-script change in this round is the oracle itself (§2).

## 5. STEP 5 — gate results (mandatory order)

| Step | Command | Result |
| --- | --- | --- |
| a | `python3 fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` | **exit 0** — 37 cases, 0 violations (`status: pass`); **exit 1 before** the fix (5 violations, §2.1) |
| b1 | `python3 fastgatk-native/scripts/verify_genotype_gvcf.py` | **exit 0** — `{"status": "pass", "output_records": 1}` (`.diag/nda-verify-gvcf.log`) |
| b2 | `ctest --test-dir fastgatk-native/build -R 'genotype-gvcf' -V` | **exit 0** — **15/15 passed**, 514.27 s (`.diag/nda-ctest-genotype-gvcf.log`) |
| c | the 17-name strict-gate filter of the task | **exit 0** — **19/19 passed**, 1413.38 s (`.diag/nda-strict-gates.log`) |
| d | `fastgatk-native/scripts/run_regression.sh --label nda-annotation` | **exit 0** — **300/300 on BOTH backends** (`.diag/regression/20260912-114941/{omp,serial}.log`), 0 staleness warnings |
| extra | `python3 fastgatk-native/scripts/verify_genotype_gvcf_max_alternate_alleles_gatk_oracle.py` | **exit 0** (`.diag/nda-maxalt-oracle.log`) — the flag's multi-ALT/max-ALT path |

### 5.1 Dual-backend regression (d)

```
REG_EXIT=0
| backend | build dir                             | result | passed | time        |
| omp     | OpenMP (fastgatk-native/build)        | pass   | 300/300 | 1296.78 sec |
| serial  | Serial (fastgatk-native/build-serial) | pass   | 300/300 | 1315.81 sec |
```

Run window 11:49:41 → 12:11:37, label `nda-annotation`, git `d4faf64` + the two
uncommitted files of §8, `FASTGATK_REQUIRE_GATK_ORACLE=1` (the runner's default),
ctest parallelism 8.

### 5.2 What was edited after the suite, and the proof it changed nothing

`run_regression.sh` finished at 12:11:37.  Afterwards I corrected **ten
imprecise `file:line` citations** in the prose I had added (they pointed at
`GenotypingEngine.java:140-146` → `:137-144`,
`ReferenceConfidenceVariantContextMerger.java:341-343` → `:340-342`,
`GenotypeGVCFsEngine.java:139` → `:148`/`:154` for the non-variant branch, and
`GenotypeGVCFsEngine.java:234-241` → `:233-250`).  These are comment/docstring
`why`-text only.  Proof, per the task's rule:

```
before edit (tested by the suite):
  4f9943c601227c1f8c1e49513f55e762  fastgatk-native/build/fastgatk-genotype-gvcf
  b179c11d01fd6b154e18d250eb90f0ec  fastgatk-native/build-serial/fastgatk-genotype-gvcf
both trees rebuilt after the edit (12:11:58 / 12:12:09):
  4f9943c601227c1f8c1e49513f55e762  fastgatk-native/build/fastgatk-genotype-gvcf
  b179c11d01fd6b154e18d250eb90f0ec  fastgatk-native/build-serial/fastgatk-genotype-gvcf
  -> BYTE-IDENTICAL (diff of .diag/nda-binaries-{before,after}-citation-fix.md5 is empty)
```

The binaries the suite tested are byte-identical to the surviving ones.  As an
extra check on the script-side prose edit, step 5(a) was **re-run** against the
final tree: exit 0, 37 cases, 0 violations
(`.diag/nda-gate-after-final.log`).  Steps 5(b), 5(c) and 5(d) were **not**
re-run, because no compiled byte changed.

## 6. What remains unproven

* The gate pins the fixtures it runs (6 new cases over 5 fixture shapes, 37 cases
  total).  §1's probes cover 13 shapes, and the source rule is read, but this is
  a point measurement, not a proof for every input shape.  Untested by a gate:
  a locus where two input records at the same position carry *different*
  concrete ALTs (measured: GATK emits two rows, one per ALT group, each
  `NDA=1` — `.diag/nda-probe2.log`, `cross-sample-union`; native does not
  reproduce that grouping and it is a pre-existing divergence unrelated to
  `NDA`), `.bcf`/`.vcf.gz` inputs, a GenomicsDB source, and `--dbsnp`.
* The "empty merged ALT set ⇒ no `NDA`" half is measured only through the dense
  path (`.diag/nda-probe2.log`, `homref-nonref-only-dense`).  I read
  `regenotypeVC()`'s branch structure but did not trace every route by which a
  locus with no merged ALT can still reach `composeCallAttributes()`.
* `record.alleles` is assumed to be exactly GATK's merged `vc` allele list
  (`<NON_REF>` removed, `*` retained) at the point of the call.  That assumption
  is confirmed by the six gated cases (including the `*`-counting and REF-only
  dense ones) but not by reading every producer of `record.alleles` in the
  ~7000-line tool.
* The three `max-alt`-reduction statements rest on the source read plus the
  existing registered oracle (`NDA=7` with `--max-alternate-alleles 6`); my own
  probe case for it produced no records (a malformed PL width in the probe
  fixture), so it was not measured directly this round.
* `(d) order-only` header differences and the `##FILTER` group's absolute index
  remain the pre-existing whole-header ordering divergence recorded by the
  previous round; untouched and still not gated.
* "No registered test pinned the old behaviour" (§4) is a grep plus the step-5
  results, not a proof: a test could construct such an expectation dynamically.
  300/300 on both backends **with the flag-bearing max-ALT oracle included** is
  the empirical check.

## 7. Tree state

* Change **kept in the tree** (no commit, no branch, no branch switch).
* Modified — exactly two files:
  `fastgatk-native/src/genotype_gvcf_tool.cpp`
  (`md5 ef29462fdb29c11b1c3f022e500c0e2c`, +29/-6) and
  `fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`
  (`md5 f99c4e8369c1d48c0209a76ead116a2c`, +234/-0).
* NOT touched: `fastgatk-native/CMakeLists.txt`, every root `*.md`, Mutect2,
  every other tool, and every registered test script other than the oracle.
* Binaries (both rebuilt after the last edit, both byte-identical to the ones the
  suite ran):
  `4f9943c601227c1f8c1e49513f55e762  fastgatk-native/build/fastgatk-genotype-gvcf`,
  `b179c11d01fd6b154e18d250eb90f0ec  fastgatk-native/build-serial/fastgatk-genotype-gvcf`.
* Artifacts: probes `.diag/nda_probe.py`, `.diag/nda_probe2.py`,
  `.diag/nda_gate_probe.py` (+ `.diag/nda-probe{,2}.log`,
  `.diag/nda-gate-probe-before.log`); gate logs
  `.diag/nda-gate-{before,after-cases,after,after-final}.log`; build logs
  `.diag/nda-build-{omp,serial}{,-citation}.log`; step-5 logs
  `.diag/nda-{verify-gvcf,ctest-genotype-gvcf,strict-gates,regression}.log|out`
  and `.diag/nda-maxalt-oracle.log`; regression evidence
  `.diag/regression/20260912-114941/`; hash files
  `.diag/nda-binaries-{before,after}-citation-fix.md5`.
