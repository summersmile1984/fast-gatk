# Round: the dense-mode non-PASS row's FILTER column (`RGQ` where GATK writes `.`)

Scope: the residue the previous round reported without investigating
(`.diag/round-qual-infinity.md` §8): on the fixture of
`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`, native
writes `FILTER=RGQ` on the dense-mode (`--include-non-variant-sites`) REF-only
row whose source leaf is non-PASS, while pinned GATK 4.6.2.0 writes `.`.

## 0. Verdict

* **GATK truth (measured, pinned 4.6.2.0):** the row is

  ```
  chr1	3	.	A	.	Infinity	.	DP=7;MLEAC=.;MLEAF=.	GT	./.
  ```

  FILTER is the literal `.` (unfiltered), not `LowQual`, even though the source
  leaf is `FILTER=LowQual`.  GATK's output header for that run is
  `##FILTER=<ID=LowQual,Description="Low quality">`; native's compat header has
  **no** `##FILTER` line at all (§1.3, §7).
* **Mechanism (two layers, both measured):**
  1. **presence** — native never rebuilds the record's FILTER state, so the
     source leaf's FILTER is *inherited* by the output record, whereas GATK
     rebuilds the call from `new VariantContextBuilder(callSourceString(),
     vc.getContig(), vc.getStart(), vc.getEnd(), outputAlleles)`
     (`GenotypingEngine.java:181`), a constructor that copies nothing and leaves
     `filtersWereApplied == false`, so htsjdk writes `.`
     (`VCFEncoder.getFilterString()`);
  2. **name** — the inherited `d.flt[]` entry is an index into the ID dictionary
     of the header the record was **read** with, but it is resolved at
     `vcf_format1`/`bcf_write` time against the **output** header.  HTSlib keeps
     one dictionary shared by `##FILTER`/`##INFO`/`##FORMAT` and
     `vcf_parse_filter()` auto-registers an undeclared source FILTER *while the
     record is parsed*, i.e. after the output header was duplicated from the
     header-only first pass — so the index names the first id appended to the
     output header, `RGQ`.  The token is therefore **not** a leaked FORMAT key
     name in any semantic sense: it is arbitrary dictionary content at a stale
     index.
* **Previous round's guess: CONFIRMED in shape, corrected in detail.**  It
  guessed "native appears to carry the source record's FILTER **id** into a
  header whose FILTER dictionary has different indices".  The carried thing is
  an id/index as guessed, but (a) there is no separate "FILTER dictionary" —
  HTSlib has a single shared ID dictionary, which is *why* a FORMAT id (`RGQ`,
  and `GQ` in the isolation experiment) can appear in the FILTER column; (b)
  the index mismatch is not in the source of the divergence but in its
  *spelling*: with the source FILTER declared in the input header native writes
  `LowQual` — still wrong, because GATK writes `.` — so "dictionary index
  off-by-one" alone would not have fixed the row.
* **Reach:** measured in the dense path (single and two samples), in the
  ordinary non-dense variant path, through **both** writers (aggregate and
  `--stream-by-locus`), and in the non-compatibility diagnostic profile
  (out of scope for GATK parity, left unfixed, §7).
* **Fix kept in the tree:** one local change in
  `fastgatk-native/src/genotype_gvcf_tool.cpp`
  (`apply_gatk_annotation_compatibility`, `:4170-4200`): the "canonicalize a
  lone PASS to `.`" special case became "clear every inherited filter", citing
  `GenotypingEngine.java:181` and the measured aliasing.  39 changed lines in
  that file (`+28/-11`), of which 9 are code (a 6-line PASS special case
  replaced by a 3-line clear) and the rest is comment; no restructuring, no
  kernel ABI, no `CMakeLists.txt`, no Mutect2, no other tool.
* Oracle: **exit 1 before** the fix (5 FILTER-only violations, literal rows in
  §3.1), **exit 0 after** (24 cases, 0 violations).  Five new **GATED** cases in
  the already-registered oracle + one REPORTED-ONLY case pinning the half of the
  contract native does not implement (`LowQual`).
* Step 5: (a) exit 0; (b) B1 exit 0, B2 15/15 exit 0; (c) 19/19 exit 0;
  (d) **300/300 on both backends**, `REG_EXIT=0`.  Change **kept in the tree**.

## 1. STEP 1 — measured truth, byte for byte

### 1.1 Reproducer

`.diag/filter_rgq_probe.py` (log `.diag/filter-rgq-probe.log`, JSON
`.diag/filter-rgq-probe.json`) and `.diag/filter_rgq_probe2.py` (logs
`.diag/filter-rgq-probe2.log`, `.diag/filter-rgq-probe2-after.log`) build the
oracle's `STAR_ONLY_COVERED_NONPASS_RECORD` body with the oracle's own 100 bp
`chr1` reference, index it with `gatk IndexFeatureFile`, and run

```
GATK:   java -jar gatk-package-4.6.2.0-local.jar GenotypeGVCFs -R ref.fa -V in.g.vcf \
          --include-non-variant-sites -O out.vcf --create-output-variant-index false
native: fastgatk-genotype-gvcf -R ref.fa -V in.g.vcf --include-non-variant-sites \
          --gatk-compatible-annotations -O out.vcf
```

Each row is printed as **raw bytes** and each FILTER column as
`row.split("\t")[6]`, so the token itself is the evidence.

### 1.2 Literal rows, FILTER token, both headers' FILTER lines

| probe case | GATK row | native row | GATK `##FILTER` header | native `##FILTER` header |
| --- | --- | --- | --- | --- |
| **dense non-PASS, undeclared (the reported case)** | `chr1 3 . A . Infinity . DP=7;MLEAC=.;MLEAF=. GT ./.` | … `Infinity` **`RGQ`** `DP=7;…` | `##FILTER=<ID=LowQual,Description="Low quality">` | *(none)* |
| dense non-PASS, `##FILTER=<ID=LowQual,…>` declared in the input | same as above | … `Infinity` **`LowQual`** `DP=7;…` | both lines (GATK's standard line + the input's) | the input's declaration only |
| dense non-PASS, input token renamed `StrandBias` (undeclared) | same | … `Infinity` **`RGQ`** `DP=7;…` | `##FILTER=<ID=LowQual,Description="Low quality">` | *(none)* |
| dense **PASS** (control) | `chr1 3 . A . Infinity . DP=20;…` | identical | `##FILTER=<ID=LowQual,Description="Low quality">` | *(none)* |
| ordinary variant path, input `LowQual` undeclared | `chr1 2 . A G 82.26 . AC=1;… GT:AD:DP:PL 0/1:0,20:20:0,0,0` | … `82.26` **`RGQ`** `AC=1;…` | `##FILTER=<ID=LowQual,Description="Low quality">` | *(none)* |
| ordinary variant path, input `LowQual` declared | same | … `82.26` **`LowQual`** `AC=1;…` | GATK standard + the input's declaration | the input's declaration only |
| dense non-PASS, two samples (stream and aggregate) | `chr1 3 . A . Infinity . DP=7;MLEAC=.;MLEAF=. GT ./. ./.` | … `Infinity` **`RGQ`** … | as above | *(none)* |
| dense **weak locus** (QUAL 23.14) | `chr1 2 . A . 23.14 LowQual DP=20;MLEAC=.;MLEAF=. GT ./.` | … `23.14` **`.`** … | as above | *(none)* |

The last row is the *only* measured case where GATK's FILTER is not `.`: it is
GATK's own `LowQual`, and native writes `.` there because it does not implement
the rule at all (`LowQual` with an input token of `PASS` is already wrong today,
so this is a pre-existing gap, not a consequence of the inherited-filter
defect).  It is pinned as REPORTED ONLY in the gate.

### 1.3 Why `RGQ` specifically — the aliasing experiment (falsifiable prediction)

The prediction: the token is the ID that occupies the inherited index in the
**output** header; the output header is `bcf_hdr_dup()` of the header-only first
read plus the ids appended by `add_genotype_output_header_fields()`
(`genotype_gvcf_tool.cpp:5035`; `:5041` appends FORMAT/GQ when absent, `:5044`
FORMAT/RGQ), while
the source record was decoded with a header in which HTSlib had auto-registered
the undeclared `FILTER` (`third_party/htslib-build/htslib-src/vcf.c:3549-3556`:
`ksprintf(&tmp, "##FILTER=<ID=%s,Description=\"Dummy\">", t)`).  The fixture
header declares `FORMAT/GQ`, so the **first** id appended to the output header
is `RGQ` → the token is `RGQ`.  Remove the `FORMAT/GQ` declaration from the
input header and the first appended id becomes `GQ` → the token must become
`GQ`.

Measured, before the fix (`.diag/filter-rgq-probe2.log`):

```
=== dense-nonpass-header-without-GQ
    NATIVE FILTER: ['.', 'GQ']      <- prediction confirmed
```

and, decisively, an input FILTER with a completely different name still renders
as the same token:

```
=== dense-othertoken-undeclared   (input FILTER=StrandBias)
    GATK   FILTER column: ['.', '.']
    NATIVE FILTER column: ['.', 'RGQ']
```

Three facts together identify the mechanism:
1. the input token's *name* does not reach the output (`StrandBias` → `RGQ`);
2. the token is the first appended output-header id (`RGQ`, and `GQ` when `GQ`
   is the first appended id);
3. declaring the source FILTER in the input header removes the wrongness of the
   *name* (`LowQual`), which proves the index is resolved against a different
   dictionary.

### 1.4 GATK's rule for this row, with source

```java
// GenotypingEngine.calculateGenotypes(), GenotypingEngine.java:181-186
final VariantContextBuilder builder = new VariantContextBuilder(
        callSourceString(), vc.getContig(), vc.getStart(), vc.getEnd(), outputAlleles);
builder.log10PError(log10Confidence);
if ( ! passesCallThreshold(phredScaledConfidence) ) {
    builder.filter(GATKVCFConstants.LOW_QUAL_FILTER_NAME);
}
```

* **Nothing is copied from the source VariantContext.**  This is the 5-arg
  constructor, not `VariantContextBuilder(VariantContext)`: disassembled from
  the pinned jar, it initialises `filters = null` (and `ID = "."`,
  `log10PError = 1.0`, `attributes = null`) and then stores only
  source/contig/start/stop — so the built context is unfiltered and
  `filtersWereApplied` is false unless `filter(...)` is called at `:184-186`.
* **htsjdk renders `.` exactly for that state** — disassembled
  `htsjdk.variant.vcf.VCFEncoder.getFilterString()` from the pinned jar:

  ```
  if (vc.isFiltered())      -> join sorted filters with ";"
  else if (vc.filtersWereApplied()) -> "PASS"
  else                      -> "."
  ```

  so `filtersWereApplied == false` ⟹ `.`, which is what is measured.
* **The only filter GATK can apply here is its own `LowQual`**, and the test is
  `!passesCallThreshold(phredScaledConfidence)`
  (`GenotypingEngine.java:430-431`: `conf >= standardConfidenceForCalling`,
  default 30) on the **recomputed** confidence — never on the source FILTER.
  For this row the recomputed confidence is the `Infinity` QUAL
  (`.diag/round-qual-infinity.md`), which passes, so no filter is applied.
* Consequently `LowQual` in the input has **no** influence on the output row:
  the same fixture with `FILTER=PASS` and with `FILTER=LowQual` produce byte
  identical GATK output (measured, §1.2 rows 1 and 4).

### 1.5 The responsible code line

The **only** filter-writing call site in the whole tool is
`apply_gatk_annotation_compatibility()` in
`fastgatk-native/src/genotype_gvcf_tool.cpp` (there is exactly one
`bcf_update_filter` in the translation unit; the genotype engine never sets a
filter otherwise), and it only *removed* one:

```cpp
// before (was :4170-4182)
if (!options.gatk_annotation_compatibility) return;
// HTSJDK's VariantContext writer renders an unfiltered call as `.` even
// when the incoming gVCF encoded the equivalent state as FILTER=PASS.
// Preserve real filter labels, but canonicalize this PASS spelling at the
// final compatibility writer boundary.
bcf_unpack(record.value, BCF_UN_FLT);
if (record.value->d.n_flt == 1) {
    const auto* filter = bcf_hdr_int2id(output_header, BCF_DT_ID, record.value->d.flt[0]);
    if (filter != nullptr && std::strcmp(filter, "PASS") == 0 &&
        bcf_update_filter(output_header, record.value, nullptr, 0) != 0)
        throw ...("cannot canonicalize PASS filter");
}
```

So a non-PASS source filter survived to the writer, and the index was resolved
there against `output_header` — the `RGQ` token.  The record reaches the writer
through two encode stages, both of which call this function first
(`:5718`/`:5738` in the streaming path, `:6427`/`:6447` in the aggregate path).

## 2. STEP 2 — the gate, written and run BEFORE the fix

Extended the already-registered
`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` (rather
than adding a new script) because (i) its `CASES` structure already expresses
exactly these fixtures and compares data rows byte-for-byte, and (ii) it is
already registered in `fastgatk-native/CMakeLists.txt:1556-1560` under the
strict-gate filter of step 5(c) — `CMakeLists.txt` must not be edited in this
round, so a new file could not have been registered by me.

Added:

* `HEADER_WITH_LOWQUAL_FILTER` — the same header with
  `##FILTER=<ID=LowQual,Description="Low quality">` declared, so the source
  FILTER is registered while the header is *read*;
* `STAR_ONLY_COVERED_NONPASS_TWO_SAMPLES_RECORD`, `G_PLAUSIBLE_NONPASS_RECORD`,
  `WEAK_LOCUS_RECORD` fixtures;
* measured row constants `GATK_STAR_ONLY_DENSE_NONPASS_ROW`,
  `GATK_STAR_ONLY_DENSE_NONPASS_TWO_SAMPLES_ROW`, `GATK_G_NONPASS_ROW`,
  `GATK_WEAK_LOCUS_LOWQUAL_ROW`;
* **GATED** `covered-star-only-record-dense-non-pass` (the reported case),
  **GATED** `covered-star-only-record-dense-non-pass-two-samples`, **GATED**
  `covered-star-only-record-dense-non-pass-stream-by-locus` (the other writer),
  **GATED** `non-pass-variant-input-undeclared-filter` (ordinary variant path,
  wrong-name half), **GATED** `non-pass-variant-input-declared-filter` (ordinary
  variant path, right name — still must be `.`);
* REPORTED ONLY `weak-locus-lowqual-filter-not-implemented` (GATK's own
  `LowQual`, the unimplemented half);
* `run_case()` honours an optional per-case `"native_args"` key (native-only
  flags, used by the `--stream-by-locus` case); no pre-existing case sets it;
* a "The FILTER column" docstring section quoting `GenotypingEngine.java:181`
  and `:184-186`, htsjdk's `getFilterString()` rule and the measured aliasing
  experiment.

### 2.1 Exit status and literal rows

**BEFORE the fix** — **exit 1**, `"status": "divergence"`, **5 violations**,
every one of them the FILTER column and nothing else
(`.diag/filter-rgq-gate-before.log`):

```
[covered-star-only-record-dense-non-pass] row 1 is not byte-identical:
  GATK  ='chr1\t3\t.\tA\t.\tInfinity\t.\tDP=7;MLEAC=.;MLEAF=.\tGT\t./.'
  NATIVE='chr1\t3\t.\tA\t.\tInfinity\tRGQ\tDP=7;MLEAC=.;MLEAF=.\tGT\t./.'
[covered-star-only-record-dense-non-pass-two-samples] row 1 is not byte-identical:
  GATK  ='chr1\t3\t.\tA\t.\tInfinity\t.\tDP=7;MLEAC=.;MLEAF=.\tGT\t./.\t./.'
  NATIVE='chr1\t3\t.\tA\t.\tInfinity\tRGQ\tDP=7;MLEAC=.;MLEAF=.\tGT\t./.\t./.'
[covered-star-only-record-dense-non-pass-stream-by-locus] row 1 is not byte-identical:
  GATK  ='chr1\t3\t.\tA\t.\tInfinity\t.\tDP=7;MLEAC=.;MLEAF=.\tGT\t./.'
  NATIVE='chr1\t3\t.\tA\t.\tInfinity\tRGQ\tDP=7;MLEAC=.;MLEAF=.\tGT\t./.'
[non-pass-variant-input-undeclared-filter] row 0 is not byte-identical:
  GATK  ='chr1\t2\t.\tA\tG\t82.26\t.\tAC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11\tGT:AD:DP:PL\t0/1:0,20:20:0,0,0'
  NATIVE='chr1\t2\t.\tA\tG\t82.26\tRGQ\tAC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11\tGT:AD:DP:PL\t0/1:0,20:20:0,0,0'
[non-pass-variant-input-declared-filter] row 0 is not byte-identical:
  GATK  ='chr1\t2\t.\tA\tG\t82.26\t.\tAC=1;…\tGT:AD:DP:PL\t0/1:0,20:20:0,0,0'
  NATIVE='chr1\t2\t.\tA\tG\t82.26\tLowQual\tAC=1;…\tGT:AD:DP:PL\t0/1:0,20:20:0,0,0'
```

The 17 previously gated cases stayed green **before** the fix (including the
dense `Infinity` QUAL cases and the PASS controls), so the five violations
isolate the FILTER column.  The GATK-truth assertions (`expect`) also passed
before the fix — the measured rows did not move.

**AFTER the fix: exit 0**, `"status": "pass"`, 0 violations, 24 cases
(`.diag/filter-rgq-gate-after.log`): 22 gated + 2 reported-only, 31 GATK rows
compared byte for byte.

## 3. STEP 3 — the fix

`git diff --stat` = **2 files changed, 240 insertions(+), 11 deletions(-)**
(oracle +212, tool +28/-11 — of the tool delta only two lines are behaviour).
The whole behaviour change in `fastgatk-native/src/genotype_gvcf_tool.cpp`
(`apply_gatk_annotation_compatibility`, now `:4170-4200`):

```diff
     if (!options.gatk_annotation_compatibility) return;
-    // HTSJDK's VariantContext writer renders an unfiltered call as `.` even
-    // when the incoming gVCF encoded the equivalent state as FILTER=PASS.
-    // Preserve real filter labels, but canonicalize this PASS spelling at the
-    // final compatibility writer boundary.
+    // ... (comment quoting GenotypingEngine.java:181 and :184-186,
+    //      VCFEncoder.getFilterString(), and the measured RGQ/GQ aliasing)
     bcf_unpack(record.value, BCF_UN_FLT);
-    if (record.value->d.n_flt == 1) {
-        const auto* filter = bcf_hdr_int2id(output_header, BCF_DT_ID,
-                                            record.value->d.flt[0]);
-        if (filter != nullptr && std::strcmp(filter, "PASS") == 0 &&
-            bcf_update_filter(output_header, record.value, nullptr, 0) != 0)
-            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot canonicalize PASS filter");
-    }
+    if (record.value->d.n_flt > 0 &&
+        bcf_update_filter(output_header, record.value, nullptr, 0) != 0)
+        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear inherited FILTER");
```

Why this is the correct and minimal fix:

* It implements GATK's rule for the FILTER column as far as the contract is
  observable: **the output record carries no inherited filter**, so the column
  is `.` exactly when GATK's engine applies no filter.  `.` is what GATK writes
  for every measured row except the `LowQual` one it decides itself.
* It subsumes the old PASS canonicalization (PASS → `.` is still produced, by
  the same call), and it additionally fixes the multi-filter case
  (`q10;LowQual`) that the old special case left alone.
* It is local to the `--gatk-compatible-annotations` boundary: the function
  returns early otherwise, so the non-compatibility diagnostic profile, every
  other tool, the Kokkos kernels and the generated headers are untouched.
* **Not** a name-mapping patch: fixing only the index (translating the filter
  by name) would have produced `FILTER=LowQual` on a row where GATK writes `.`
  (measured, §1.2 row 2) — i.e. the "obvious" off-by-one repair is refuted by
  measurement.

Not restructured: the genotyping engine, the AF/PL Kokkos kernels and their
ABI, the union/merge stages, the output-allele-subset predicate, the dense
REF-only materialization, the QUAL/`Infinity` handling of the previous round,
and every other tool.  `fastgatk-native/CMakeLists.txt` and every root `*.md`
are unmodified.

### 3.1 Shared code / collateral, explicitly

* The changed function is in the **anonymous namespace** of
  `genotype_gvcf_tool.cpp` (`namespace {` at `:48`) and is called from exactly
  four places, all in that file's two encode stages (`:5718`, `:5738`, `:6427`,
  `:6447`).  No other tool and no other translation unit can be affected.
* Header dictionary: **no matching change was needed**.  Clearing a filter
  cannot make the output header inconsistent, and native never emits `LowQual`
  itself, so no `##FILTER` line has to exist for the fix to work.  (GATK's
  output header does carry `##FILTER=<ID=LowQual,Description="Low quality">`
  even when the input never mentions it; native's does not — a *header-level*
  divergence, measured and deliberately left alone, §7.)
* Collateral check on the compat surface: the only rows whose FILTER changes
  are rows with an inherited non-PASS filter, which were wrong before; rows
  whose source FILTER is `PASS` already rendered `.` and still do.

## 4. STEP 4 — stale assertions

**No registered test pinned the old behaviour; no test line needed correcting.**

| File / check | Why it could have been stale | Outcome |
| --- | --- | --- |
| every `fastgatk-native/scripts/*.py` | a test could pin `RGQ` as a FILTER token | 7 `RGQ` hits, **all FORMAT-field names**: `verify_gatk_genotype_gvcf.py:264,270` (FORMAT/RGQ), `verify_genotype_gvcf_include_non_variant_gatk_oracle.py:73` (`FORMAT == "GT:DP:RGQ"`), `verify_hc_gvcf_symbolic_prior_gatk_oracle.py:204` and `hc_symbolic_prior_fixture_lib.py:239` (field-name lists).  None is a FILTER assertion |
| every `*.py` / `*.cpp` under `fastgatk-native/` | a fixture could feed a non-PASS FILTER column to the genotype tool | `grep` for an escaped or literal-tab `LowQual`/`StrandBias` **FILTER** column in a genotype fixture returns **nothing**; every genotype fixture's FILTER column is `PASS` (except `verify_genotype_gvcf_spandel_gatk_oracle.py:351`/`:361`, whose low-quality leaf is dropped by both tools in the modes that were gated before this round) |
| the genotype oracles that compare whole rows (`verify_gatk_genotype_gvcf*.py`, `verify_hc_dense_gvcf_genotype_gatk_oracle.py`, `verify_genotype_gvcf_include_non_variant_gatk_oracle.py`, `..._legacy_qual.py`) | they compare data rows, so the FILTER column is implicitly compared | no `LowQual`/`filter`/`[6]` expectation in any of them |
| `fastgatk-native/src/reblock_gvcf_tool.cpp:1814-1815`, `genotype_gvcf_tool.cpp` RGQ sites | native's own FORMAT/RGQ is unrelated to the FILTER column | untouched, and step 5(d) re-ran their tests |

The only test-script change in this round is the oracle itself (§2).

## 5. Diff summary

```
 fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py | 212 +++++++++++++++++++++
 fastgatk-native/src/genotype_gvcf_tool.cpp                         |  39 ++--
 2 files changed, 240 insertions(+), 11 deletions(-)
```

`git status --short` (nothing else in the tree):

```
 M fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py
 M fastgatk-native/src/genotype_gvcf_tool.cpp
```

## 6. STEP 5 — gate results (mandatory order)

| Step | Command | Result |
| --- | --- | --- |
| a | `python3 fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` | **exit 0** (24 cases: 22 gated + 2 reported-only, 31 rows, 0 violations); **exit 1 before** the fix (5 violations, §2.1) |
| b1 | `python3 fastgatk-native/scripts/verify_genotype_gvcf.py` | `B1_EXIT=0` — `{"status": "pass", "output_records": 1}` (`.diag/filter-rgq-verify-gvcf.log`) |
| b2 | `ctest --test-dir fastgatk-native/build -R 'genotype-gvcf' -V` | `B2_EXIT=0` — **15/15 passed**, incl. `fastgatk-genotype-gvcf-spandel-gatk-oracle` (`.diag/filter-rgq-ctest-genotype-gvcf.log`) |
| c | the 17-name strict-gate filter from the task | `C_EXIT=0` — **19/19 matched tests passed**, incl. `fastgatk-genotype-gvcf-spandel-gatk-oracle` at the final revision, 208.30 s (`.diag/filter-rgq-strict-gates.log`) |
| d | `fastgatk-native/scripts/run_regression.sh --label filter-rgq` | see below |

`fastgatk-native/scripts/run_regression.sh --label filter-rgq`
(omp + serial, `FASTGATK_REQUIRE_GATK_ORACLE=1` by default, ctest parallelism 8):

```
REG_EXIT=0
| 后端   | 构建目录                              | 结果 | 通过/总数 | 耗时        |
| omp    | OpenMP (fastgatk-native/build)        | 通过 | 300/300   | 1319.95 sec |
| serial | Serial (fastgatk-native/build-serial) | 通过 | 300/300   | 1318.66 sec |
```

Evidence block `.diag/regression/20260912-060534/{omp.log,serial.log,evidence.md}`,
run window 06:05:34 → 06:27:34, label `filter-rgq`, git `b2d87bc` + the two
uncommitted files of §5, and the runner reports its staleness check as clean.

**Nothing was edited after the suite run**, so no rebuild or re-run was needed;
`md5sum -c .diag/filter-rgq-binaries.md5` re-confirms

```
9cafa54acc60c4f92c11ba6e06c42c37  fastgatk-native/build/fastgatk-genotype-gvcf        (omp)
be6faf9b3604b0bec1a02c682a2c8c2b  fastgatk-native/build-serial/fastgatk-genotype-gvcf (serial)
```

and both binaries were built after the last source edit and before every gate
above (logs `.diag/filter-rgq-build-{omp,serial}.log`).  Timestamps, all before
the suite window 06:05:34 → 06:27:34: oracle 05:27:58 < tool
`genotype_gvcf_tool.cpp` 05:31:54 < omp binary 05:32:16 < serial binary
05:32:27.  The only file touched afterwards is this report, which is not a test
input (no source, test script, header or build input changed after the suite).

## 7. How far the defect reaches, and what is deliberately left

Measured, with the tree at the fixed revision where noted:

| Surface | Before | After |
| --- | --- | --- |
| dense (`--include-non-variant-sites`) non-PASS REF-only row | `RGQ` | `.` ✅ |
| dense non-PASS, two samples | `RGQ` | `.` ✅ |
| dense non-PASS via `--stream-by-locus` (other writer) | `RGQ` | `.` ✅ |
| ordinary non-dense variant path, undeclared header | `RGQ` | `.` ✅ |
| ordinary non-dense variant path, declared header | `LowQual` | `.` ✅ |
| `FORMAT/GQ` already declared (isolation experiment) | `GQ` | `.` ✅ |
| dense PASS rows (control) | `.` | `.` ✅ |
| **non-compatibility native diagnostic profile** (no `--gatk-compatible-annotations`) | `RGQ` (and `PASS` kept) | **`RGQ` — unchanged, NOT fixed** |
| **GATK's own `LowQual`** (weak locus, QUAL 23.14) | `.` | **`.` — unchanged, NOT fixed** |
| **output header `##FILTER` lines** | no `##FILTER` at all vs GATK's `##FILTER=<ID=LowQual,Description="Low quality">` | **unchanged, NOT fixed** |

* **Non-compat profile (residual, reported).**  The aliasing still renders the
  inherited filter as `RGQ` there, because the fix is scoped to the
  GATK-compatibility writer.  That profile is explicitly a native diagnostic
  output, not a GATK parity surface, and it *intends* to carry filters through
  (it keeps `PASS`, which the compat boundary canonicalizes).  The proper repair
  there is name-based translation of the inherited filter (resolve
  `d.flt[]` against the reading header and re-register by name), which is a
  change in the record-merge/write path shared by the tool — **not minimal**,
  and not needed for the reported divergence.  Recommendation: do it in the
  round that adopts a per-record "reading header" handle, together with a gate
  on the non-compat profile.
* **GATK's own `LowQual` (residual, reported).**  Native never applies
  `GenotypingEngine.java:184-186`, so a locus whose recomputed confidence is
  below `--standard-min-confidence-threshold-for-calling` gets `.` instead of
  `LowQual`.  Pinned as REPORTED ONLY.  Implementing it needs the *pre-posterior*
  phred confidence at the writer boundary plus a `##FILTER=<ID=LowQual,…>` header
  line, i.e. it is a separate change with a much wider blast radius than this
  round's fix (every emitted record's FILTER column) — **not minimal**, so it
  was left out; the residual is now measured and gated as reported-only instead
  of being unknown.
* **Header `##FILTER` lines (residual, reported).**  GATK's output header
  always declares `LowQual` (its standard header line), even for a run whose
  input never mentions filters (`##FILTER=<ID=LowQual,Description="Low quality">`
  measured for the dense PASS control too); native's compat header has no
  `##FILTER` line unless the input declared one.  The gate compares data rows
  (as required), so this is not gated.  *Speculation:* adding the line would be
  a one-line header change, but it touches every compat output header and is
  unrelated to the reported column; I did not measure its collateral.

## 8. What remains unproven

* The gate asserts byte-identical rows on the fixtures it pins — including the
  FILTER column now — but it is a point measurement, not a proof that native's
  FILTER state equals GATK's for every input shape.  Untested shapes:
  multi-filter source rows (`q10;LowQual`), `a;b`-style filters whose names are
  declared, `--force-output-intervals`, `--max-alternate-alleles` reduction
  combined with a non-PASS leaf, multi-shard inputs, and `.vcf.gz` inputs.
  *Speculation:* the fix is a clear-all at the writer boundary, so the input's
  filter spelling cannot matter, but that is not measured for those shapes.
* The "no test pins an inherited non-PASS filter on the compat surface"
  conclusion in §4 is a grep + full-suite result, not a proof: a test could
  construct such a fixture dynamically.  Step 5(d) is the empirical check
  (300/300 on both backends with `FASTGATK_REQUIRE_GATK_ORACLE=1`), which is
  evidence, not proof.
* The mechanism (§1.3) is established by three measurements and by reading
  `vcf_parse_filter`/`bcf_update_filter` semantics; I did not instrument HTSlib
  to print the dictionary index itself.
* Whether GATK's `LowQual` and native's absence of it can be aligned without
  touching the posterior-QUAL path is **unproven** and explicitly out of scope.

## 9. Tree state

* Change **kept in the tree** (no commit, no branch, no branch switch).
* Modified — exactly two files (§5).  Scratch lives under `.diag/` only, which
  is git-ignored.
* Not touched: `fastgatk-native/CMakeLists.txt`, every root `*.md`, Mutect2,
  every other tool, and every registered test script other than the oracle.
* Artifacts: gate logs `.diag/filter-rgq-gate-{before,after}.log`, probes
  `.diag/filter_rgq_probe.py`, `.diag/filter_rgq_probe2.py` (+ their logs and
  JSON), measurement log `.diag/filter-rgq-measure-gatk.log`, build logs
  `.diag/filter-rgq-build-{omp,serial}.log`, step-5 logs
  `.diag/filter-rgq-{verify-gvcf,ctest-genotype-gvcf,strict-gates}.log`,
  binary hashes `.diag/filter-rgq-binaries.md5`.
