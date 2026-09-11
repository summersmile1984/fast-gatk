# Round: VariantFiltration `AS_FilterStatus` for `--apply-allele-specific-filters`
(native `LowASQD,PASS` vs pinned GATK 4.6.2.0 `SITE|SITE`)

Date: 2026-09-11. Working tree: `08ce73d` + uncommitted changes (no commit, no branches created).
Scope: `fastgatk-variant-filtration` only. Mutect2 and every other tool untouched.

## 1. The divergence

`fastgatk-native/scripts/verify_variant_filtration.py` (pre-round lines 504-505) asserted

```python
assert "AS_FilterStatus=LowASQD,PASS" in allele_records[0][7]
assert "AS_FilterStatus=PASS,LowASQD" in allele_records[1][7]
```

for that file's own fixture (lines 481-489) and these arguments:

```
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##INFO=<ID=QD,Number=1,Type=Float,Description=Quality by depth>
##INFO=<ID=AS_QD,Number=A,Type=Float,Description=Allele quality by depth>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM  POS ID REF ALT QUAL FILTER INFO              FORMAT S1
chr1    1   .  A   C,G 50   PASS   AS_QD=1.0,3.0     GT     0/1
chr1    2   .  C   T,G 50   PASS   AS_QD=3.0,1.0     GT     0/1

VariantFiltration --filter-expression 'vc.getAttribute("AS_QD") < 2' \
                  --filter-name LowASQD --apply-allele-specific-filters
```

The audit's claim reproduces exactly: pinned GATK emits `AS_FilterStatus=SITE|SITE`, joins
alleles with `|`, and applies **no** allele filter at all. Native emitted
`AS_FilterStatus=LowASQD,PASS` / `PASS,LowASQD` with a comma separator and `PASS` for the
unfiltered allele. The test pinned native-only behaviour.

## 2. Measured GATK truth (STEP 1)

Command actually run (gate case `as-qd-number-a`; GATK needs an unindexed-readable input, and
Python's `gzip` does not produce BGZF, so the gate uses an **uncompressed** `.vcf` for both
tools — GATK otherwise aborts with *"An index is required but was not found"*):

```
third_party/jdk17/bin/java -Xmx1g -jar \
  third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar VariantFiltration \
  -V allele-specific.vcf -O gatk.vcf \
  --filter-expression 'vc.getAttribute("AS_QD") < 2' --filter-name LowASQD \
  --apply-allele-specific-filters
```

Literal GATK data rows (exit 0), versus native before the fix:

```
GATK   chr1 1 . A C,G 50 PASS AS_FilterStatus=SITE|SITE;AS_QD=1.0,3.0  GT 0/1
GATK   chr1 2 . C T,G 50 PASS AS_FilterStatus=SITE|SITE;AS_QD=3.0,1.0  GT 0/1
NATIVE chr1 1 . A C,G 50 PASS AS_QD=1,3;AS_FilterStatus=LowASQD,PASS   GT 0/1   (pre-fix)
NATIVE chr1 2 . C T,G 50 PASS AS_QD=3,1;AS_FilterStatus=PASS,LowASQD   GT 0/1   (pre-fix)
```

Additional probes against the same jar (all exit 0) that pin the rule:

| probe (same fixture, flag on) | GATK |
| --- | --- |
| `vc.getAttribute("AS_QD") < 2` | `PASS` / `SITE\|SITE` — never fires |
| `vc.hasAttribute("AS_QD")` | `PASS` / `SITE\|SITE` — the attribute is **not visible** |
| `vc.hasAttribute("DP")` | `PASS` / `SITE\|SITE` — no INFO attribute is visible |
| `vc.getGenotype("S1").isHet()` (no flag) | fires; with the flag → `PASS` / `SITE\|SITE` + `JexlEngine - attempting to call method on null` |
| `true` (constant) | `Const` / `Const\|Const` — per-allele evaluation really runs |
| `QUAL < 60` | `LowQual` / `LowQual\|LowQual` — site FILTER = intersection |
| `vc.getAttribute("AS_QD") < 2` + `--missing-values-evaluate-as-failing` | `PASS` / `SITE\|SITE` — the missing-value policy does **not** apply |
| input `AS_FilterStatus=foo\|bar` | preserved verbatim |
| at 1 ALT (biallelic) | `PASS` / `SITE` |
| no flag at all, `vc.getAttribute("AS_QD") < 2` | GATK **crashes**, exit 3, `NumberFormatException: For input string: "1.0"` |

The audit's `SITE|SITE` claim is confirmed; it does not need a different starting point.

## 3. GATK rule, with file:line

* `gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/filters/VariantFiltration.java:359-365`
  — `apply()` takes the allele path only when `--apply-allele-specific-filters` (`applyForAllele`, `:240`) is set:
  `splitMultiAllelics(variant).stream().map(vc -> filter(vc, ...))`; the non-allele path (`:367`) filters the original record.
* `VariantFiltration.java:371-378` — `splitMultiAllelics` rebuilds each single-ALT context with
  `new VariantContextBuilder("SimpleSplit", vc.getContig(), vc.getStart(), vc.getEnd(), Arrays.asList(vc.getReference(), Allele.NO_CALL)).alleles(Arrays.asList(vc.getReference(), allele)).make(true)`.
  That 5-arg constructor of **htsjdk 4.2.0** sets `attributes = Collections.emptyMap()` and
  `genotypes = GenotypesContext.NO_GENOTYPES` (verified in the pinned
  `htsjdk-4.2.0.jar` with `javap -c`: constructor stores `Collections.emptyMap()` into
  `attributes` and `NO_GENOTYPES` into `genotypes`), so the per-allele context carries **no
  INFO and no FORMAT/genotypes** — that is why `vc.getAttribute("AS_QD")` is JEXL `null` there
  and why the rule cannot fire, while `vc.hasAttribute("AS_QD")` is false.
* `VariantFiltration.java:400-438` — `filter()` evaluates the site expressions against that
  bare context (`matchesFilter` → `VariantContextUtils.match`, `:421`).
* `.../walkers/mutect/filtering/AlleleFilterUtils.java:94-122` — `addAlleleAndSiteFilters`
  assembles the result:
  * `:98-108` decode the record's own `AS_FilterStatus` (`decodeASFilters`, `:24-28`, splits on
    `|` then `,` and trims) and merge new labels through `addAlleleFilters` (`:69-82`);
  * `:99-102` a pre-existing vector of the wrong arity makes GATK return the record untouched;
  * `:104-106` otherwise seed every allele slot with `GATKVCFConstants.SITE_LEVEL_FILTERS`;
  * `:112-114` `--invalidate-previous-filters` clears site FILTER (and reseeds);
  * `:115-120` site FILTER = **intersection** of the allele filter sets (`retainAll`), PASS when empty.
* `.../utils/variant/GATKVCFConstants.java:201` — `SITE_LEVEL_FILTERS = "SITE"`; this is the
  literal placeholder for an allele no filter selected — **not** `PASS`.
* `.../walkers/annotator/AnnotationUtils.java:21` — `ALLELE_SPECIFIC_RAW_DELIM = "|"`;
  `encodeAnyASListWithRawDelim` (`:57-59`) joins the per-allele entries with `|` and
  `encodeStringList` (`:48-50`) joins the labels of one allele with `,` (`LIST_DELIMITER`).

Answers to the three questions asked:

1. **SITE vs allele name**: `SITE` is written for an ALT whose new filter set is empty
   (`AlleleFilterUtils.java:104-106`, `:69-74`); a real filter name replaces it. `|` separates
   the per-ALT entries (`AnnotationUtils.java:21`).
2. **A JEXL expression that references an allele-level annotation** is evaluated *per split
   ALT*, but against a context that has no INFO at all, so in practice it never applies.
3. **When no filter fires**: FILTER becomes `PASS` (empty intersection, `:118-120`) while
   `AS_FilterStatus` is still written as `SITE|SITE` (`:104-109`).

## 4. Strict gate (STEP 2)

New file: `fastgatk-native/scripts/verify_variant_filtration_asfilterstatus_gatk_oracle.py`
(main() → int, argparse `--native` / `FASTGATK_VARIANT_FILTRATION_BINARY`, `TemporaryDirectory`
scratch, `oracle_guard` + `FASTGATK_REQUIRE_GATK_ORACLE` handling, final JSON status line,
`--expect-divergence`).

It builds the fixture in a temp dir, runs pinned GATK and native with identical arguments and
gates `CHROM..FILTER` byte identity plus the `AS_FilterStatus` value, and additionally asserts
the literal GATK expectation measured above (so fixture/CLI drift fails closed). Cases:
`as-qd-number-a` (the test fixture), `preexisting-as-filter-status`,
`invalidate-previous-filters`, `missing-values-evaluate-as-failing` are **gated**;
`invert-filter-expression` is **reported only** (see §9).

Exit status before the fix (STEP 2 requirement, literal):

```
VIOLATION: row 0: AS_FilterStatus differs: GATK='SITE|SITE' NATIVE='LowASQD,PASS'
VIOLATION: row 1: AS_FilterStatus differs: GATK='SITE|SITE' NATIVE='PASS,LowASQD'
EXIT=1
```

Exit status after the fix: `EXIT=0`, `"status": "pass"`, `"violations": []`.

## 5. The fix (STEP 3)

One function plus its call site in `fastgatk-native/src/variant_filtration_tool.cpp`
(+105/-14 lines including comments; no header, no shared kernel, no other tool):

* `apply_allele_specific_filters` now (a) decodes an existing `AS_FilterStatus` from the input
  record and only merges labels for alleles that fired (`AlleleFilterUtils.java:98-108`),
  returning the record untouched on an arity mismatch, and reseeding when
  `--invalidate-previous-filters` is set; (b) seeds every slot with the literal `SITE`
  placeholder instead of `PASS`; (c) evaluates each `AS_*` rule against a **scratch copy of the
  record with `n_info = 0`, `n_fmt = 0`, `n_sample = 0`** — i.e. GATK's INFO-less/genotype-less
  split context (the record is `bcf_dup`'d and unpacked first, so htslib never rebuilds those
  counts: `htslib/vcf.c:3994-4059` guards `bcf_unpack` on `b->unpacked`, and
  `htslib/vcf.c:5817-5827` then reports every INFO lookup as absent) — and with
  `missing_fails = false`, because the split context holds a *present but null* JEXL value,
  not an undefined property; (d) joins the labels of one allele with `,` and the alleles with
  `|`; (e) replaces the `PASS` status string with `SITE`.
* New file-local helpers `split_as_filter_component` / `decode_as_filter_status` (guarded to
  `apply_allele_specific_filters`).
* Call site: passes `options.invalidate_previous_filters`; the now-unused `missing_fails`
  parameter was dropped from the signature, because the split-context evaluation always forces
  the missing-values policy to false (GATK's JEXL null is not an undefined property).

Not touched: the filtering engine, the mask/cluster paths, the genotype-filter path, the
`AS_*`-requires-the-flag CLI guard, `VerifyVariantFiltration`-style header text, and every
other tool.

## 6. Test lines corrected (STEP 4)

`fastgatk-native/scripts/verify_variant_filtration.py` (old → new):

| pre-round line | assertion | new line | assertion |
| --- | --- | --- | --- |
| 504 | `assert "AS_FilterStatus=LowASQD,PASS" in allele_records[0][7]` | 516 | `assert allele_records[0][7] == "AS_QD=1,3;AS_FilterStatus=SITE\|SITE"` |
| 505 | `assert "AS_FilterStatus=PASS,LowASQD" in allele_records[1][7]` | 517 | `assert allele_records[1][7] == "AS_QD=3,1;AS_FilterStatus=SITE\|SITE"` |
| 508 | `assert ...["allele_filtered"] == 2` | 522 | `assert ...["allele_filtered"] == 0` |

plus a comment block at lines 482-491 recording that the previous expectation was native-only
and citing the GATK rule (`VariantFiltration.java:359-378`,
`AlleleFilterUtils.java:94-122`, `GATKVCFConstants.java:201`). The assertions were made
**stricter** (whole INFO field equality instead of substring), not weaker.

Dependency check for the old value (grep over the tree, excluding `build*` and `work/`):
`AS_FilterStatus`/`LowASQD`/`allele_filtered` appear only in that file (and in
`fastgatk-native/README.md`, updated to state the GATK contract);
`verify_variant_filtration_gatk_oracle.py`, `..._missing_boolean_gatk_oracle.py` and
`..._set_nocall_gatk_oracle.py` never pass `--apply-allele-specific-filters`, so no other
registered gate depended on the old encoding. The neighbouring assertion at line 515
(`FILTER == "PASS"` for both records) and lines 531-532 (native still refuses `AS_*` without
the flag) are unaffected and still pass.

## 7. STEP 5 results

| step | command | result |
| --- | --- | --- |
| 5a | `python3 fastgatk-native/scripts/verify_variant_filtration_asfilterstatus_gatk_oracle.py` | **exit 0** (`status: pass`, 5 cases, 4 gated + 1 reported-only) |
| 5b | `python3 fastgatk-native/scripts/verify_variant_filtration.py` | **exit 0** (`{"status": "pass", "gatk_oracle_exact": true}`) |
| 5b | `ctest --test-dir fastgatk-native/build -R 'filtration' -V` | **exit 0**, 4/4 passed (189 s) |
| 5c | `ctest --test-dir fastgatk-native/build -R 'window-invariance-gatk-oracle\|alleles-overlap-gate-oracle\|span-del-qual-gatk-oracle\|gvcf-symbolic-prior-gatk-oracle\|af-zero-format-gatk-oracle\|arbitrary-ploidy-span-del-prior-oracle\|polyploid-gvcf-span-del-prior-oracle\|spanning-prior-genotype-gq-oracle\|multialt-owner-annotation-oracle\|mutect2-recheck\|gvcf-indel-end-gatk-oracle\|culprit-gatk-oracle\|select-variants-refonly-gatk-oracle' -V` | **exit 0**, 15/15 passed (976 s) |
| 5d | rebuild `fastgatk-variant-filtration` in `fastgatk-native/build` and `fastgatk-native/build-serial` | both linked clean (no warnings) |
| 5d | `fastgatk-native/scripts/run_regression.sh --label variant-filtration-asfilterstatus` | **exit 0 — OpenMP `build` 296/296 (1278.76 s), Serial `build-serial` 296/296 (1238.46 s)**, zero staleness warnings (re-checked manually: no `*.cpp`/`*.hpp` newer than the newest artifact in either tree). Evidence: `.diag/regression/20260911-182144/{omp,serial}.log`, `evidence.md` |

## 8. Tree state

The change **is in the working tree** (no commit, no branch):

```
 M fastgatk-native/README.md                                       (+7/-1)
 M fastgatk-native/scripts/verify_variant_filtration.py            (+17/-3)
 M fastgatk-native/src/variant_filtration_tool.cpp                 (+105/-14)
?? fastgatk-native/scripts/verify_variant_filtration_asfilterstatus_gatk_oracle.py  (372 lines)
```

`fastgatk-native/CMakeLists.txt` was **not** edited (the new gate is not registered in CTest;
it is run directly, as the orchestrator registers gates).

## 9. What remains unproven / deliberately not fixed

Measured on the same jar, same fixture, flag on — **not** covered by this round's fix, each
reported with its localization (all are the *same* key in other invocations, i.e. the allele
path itself only runs when an `AS_*` expression is present):

1. **Flag on, no `AS_*` expression** (e.g. `--filter-expression 'QUAL < 60'`, or a `--mask`):
   GATK writes `AS_FilterStatus=LowQual|LowQual` / `Mask|Mask`; native writes no
   `AS_FilterStatus` at all. Localization: `variant_filtration_tool.cpp:2252-2259`
   (`allele_specific_rules_present` guard), `:2302-2305` (header line), `:2398` (call),
   `:2449-2450` (manifest flag). The FILTER columns already agree.
2. **Site FILTER from the allele intersection** (`AlleleFilterUtils.java:115-120`) is not
   implemented: with `--invert-filter-expression` (or a composite rule that matches on the
   split context) GATK writes `FILTER=LowASQD` while native still writes `PASS`; the
   `AS_FilterStatus` value already matches. This is gate case `invert-filter-expression`,
   deliberately **reported, not gated**.
3. **Mask/cluster inside the allele path**: GATK applies the mask to each split context (which
   is why GATK's own `testMask4.vcf` shows `foo|SITE` on an indel where trimming changes the
   interval); native applies masks at site level only. Untested here.
4. **Allele methods in the split context**: GATK's per-allele context has exactly two alleles
   (`getAlternateAllele(0)` is the ALT under consideration); native keeps the full multi-allelic
   record in the scratch copy, so an expression mixing `getAlternateAllele(i)` with an `AS_*`
   field could still resolve a different allele. Untested.
5. **Untested by any gate**: `--mask` + flag, `--cluster-window-size` + flag, `--set-filtered-genotype-to-no-call`
   + flag, symbolic/spanning-deletion records, and INFO keys other than `AS_QD` in allele mode.
6. **Speculation (labelled as such)**: that GATK's INFO-less split context is an unintended
   bug rather than a design choice. The evidence for intent is only indirect (the GATK
   integration tests use `--apply-allele-specific-filters` exclusively with `--mask`), so this
   round matched the measured behaviour bug-for-bug instead of guessing. Similarly, this round
   did not verify that `bcf_hdr_remove`-style header stripping is unusable for this purpose
   beyond the observation that it did not work here (htslib's `bcf_hdr_remove` clears the hrec
   list and the `hrec` dict slot but not `bcf_idinfo_t.info[]`, `htslib/vcf.c:1348-1362`,
   `:1005-1010`, so lookups kept succeeding) — measured, then replaced by the record-copy
   approach.
