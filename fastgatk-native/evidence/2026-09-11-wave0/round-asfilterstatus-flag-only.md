# Round: VariantFiltration `--apply-allele-specific-filters` with **no** `AS_*` rule
(the trigger-only gap reported by the previous round, re-measured)

Date: 2026-09-11. Working tree: `8d54a80` + uncommitted changes (no commit, no branch).
Scope: `fastgatk-variant-filtration` only. Mutect2 and every other tool untouched.
`fastgatk-native/CMakeLists.txt` and the repo-root `*.md` files were **not** edited.

## 0. Headline: the reported framing was too narrow

The previous round localised this gap to `allele_specific_rules_present`
(`variant_filtration_tool.cpp:2252-2259`, `:2302-2305`, `:2398`, `:2449-2450`) and
described it as "native writes no `AS_FilterStatus`". Re-measuring against the
pinned jar shows the divergence is **not confined to that trigger, and not
confined to the `AS_FilterStatus` column**:

* with the flag set, GATK never filters the record itself
  (`VariantFiltration.java:359-365`): **every** site expression, the mask and the
  cluster test are evaluated once per ALT against a split context that carries no
  INFO, no genotypes, no filters and the *default* QUAL;
* therefore the **site FILTER column** also diverges, e.g.
  `--filter-expression 'DP > 5' --apply-allele-specific-filters` → GATK `PASS`,
  native (pre-fix) `HighDP`; and `'QUAL < 60'` on a record whose real QUAL is 80 →
  GATK `LowQual`, native (pre-fix) `foo`/`PASS`.

A trigger-only change (emit `AS_FilterStatus` when the flag is set) can therefore
**never** produce GATK-identical rows for the case set STEP 1 asks for (records
that fail *and* pass the expression). The fix that does is the one applied here:
route the whole flag path through GATK's per-ALT split-context filtering.

## 1. Measured GATK truth (STEP 1)

Command shape (identical for both tools; GATK needs a plain, unindexed-readable
input, so the fixtures are uncompressed `.vcf`):

```
third_party/jdk17/bin/java -Xmx1g -jar \
  third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar VariantFiltration \
  -V in.vcf -O gatk.vcf <args>
```

Fixture `in.vcf` (all probes below):

```
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##FILTER=<ID=foo,Description=A pre-existing input filter>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM  POS ID REF ALT QUAL FILTER INFO   FORMAT S1
chr1    1   .  A   C,G 50   PASS   DP=10  GT     0/1
chr1    2   .  C   T,G 80   foo    DP=10  GT     0/1
chr1    3   .  G   A   30   PASS   DP=10  GT     1/1
```

Literal GATK data rows (`CHROM POS … FILTER INFO`, exit 0 in every case):

```
A) flag + --filter-expression 'QUAL < 60' --filter-name LowQual
   chr1 1 . A C,G 50 LowQual     AS_FilterStatus=LowQual|LowQual;DP=10  GT 0/1
   chr1 2 . C T,G 80 LowQual;foo AS_FilterStatus=LowQual|LowQual;DP=10  GT 0/1   <- real QUAL 80!
   chr1 3 . G A   30 LowQual     AS_FilterStatus=LowQual;DP=10          GT 1/1

B) flag + --filter-expression 'DP > 5' --filter-name HighDP
   chr1 1 . A C,G 50 PASS     AS_FilterStatus=SITE|SITE;DP=10  GT 0/1
   chr1 2 . C T,G 80 foo      AS_FilterStatus=SITE|SITE;DP=10  GT 0/1
   chr1 3 . G A   30 PASS     AS_FilterStatus=SITE;DP=10       GT 1/1

C) flag + --mask mask.vcf --mask-name Mask          (mask covers chr1:1 only)
   chr1 1 . A C,G 50 Mask     AS_FilterStatus=Mask|Mask;DP=10  GT 0/1
   chr1 2 . C T,G 80 foo      AS_FilterStatus=SITE|SITE;DP=10  GT 0/1
   chr1 3 . G A   30 PASS     AS_FilterStatus=SITE;DP=10       GT 1/1

D) no flag + 'QUAL < 60'  (control)
   chr1 1 . A C,G 50 LowQual  DP=10  GT 0/1          <- no AS_FilterStatus at all
   chr1 2 . C T,G 80 foo      DP=10  GT 0/1
   chr1 3 . G A   30 LowQual  DP=10  GT 1/1

E) flag alone, no expression, no mask (GATK exit 0; native exit 2, see §8.1)
   chr1 1 . A C,G 50 PASS     AS_FilterStatus=SITE|SITE;DP=10  GT 0/1
   chr1 2 . C T,G 80 foo      AS_FilterStatus=SITE|SITE;DP=10  GT 0/1
   chr1 3 . G A   30 PASS     AS_FilterStatus=SITE;DP=10       GT 1/1
```

Answers to the questions STEP 1 asks:

| question | measured answer |
| --- | --- |
| per-allele entries = filter name repeated once per allele? | yes, **for every** ALT (`LowQual\|LowQual`, `Mask\|Mask`); an ALT no filter selected keeps the literal `SITE`; a biallelic record has a single entry (`LowQual`) |
| is the site FILTER set too? | yes. It is the **intersection** of the per-ALT filter sets, added to the record's existing FILTER (`LowQual;foo`); PASS when the intersection is empty |
| does the flag change *which* records the rule matches? | yes. The split context's QUAL is the constant `-10.0`, so `QUAL < 60` matches **every** record, including the QUAL=80 one; an INFO rule (`DP > 5`) matches **none** |
| no flag ⇒ AS_FilterStatus? | no. The control (D) writes no `AS_FilterStatus` and filters the record itself |
| flag with no expression at all? | `AS_FilterStatus=SITE\|SITE` is still written (E) |

Supporting probes against the same jar that pin the split-context values
(`--filter-name F --apply-allele-specific-filters`, `FILTER` column shown):

| expression | GATK | conclusion |
| --- | --- | --- |
| `QUAL < 0` | `F` | QUAL is below zero |
| `QUAL > 0` | `PASS` | QUAL is `-10.0`, not the record's QUAL |
| `QUAL == 0` | `PASS` | |
| `DP > 5` / `DP < 5` / `DP == 10` / `vc.getAttribute("DP") < 5` / `vc.hasAttribute("DP")` | `PASS` | no INFO attribute exists in the split context |
| `DP > 5` + `--missing-values-evaluate-as-failing` | `F` | a **bare identifier** that is absent makes the whole JEXL expression null, so the missing-value policy applies |
| `vc.getAttribute("DP") < 5` + `--missing-values-evaluate-as-failing` | `PASS` | a **method call** returning null compares as `false`, so the policy does *not* apply (this is the distinction native cannot currently express — §8.3) |
| `vc.isSNP()` | `F` | allele/type methods work in the split context |
| `vc.getFilters().isEmpty()` (record FILTER=foo) | `F` | the split context carries **no** filters |
| `NALLELES == 2` (3-allele record) | `F` | the split context has exactly two alleles |
| `START < 3` | `PASS` | `START` is not a JEXL attribute at all (native-only alias — §8.4) |

### 1.1 Extra differential probes (not part of the gate)

GATK and native were also compared on the same fixture for 21 further
invocations: **19 matched exactly** on `(FILTER, AS_FilterStatus)`, including
`flag + two rules`, `flag + AS_QD + QUAL`, `flag + --invalidate-previous-filters`,
`flag + --mask` (both directions), `flag + vc.isSNP()/isIndel()/isBiallelic()`,
`flag + vc.getFilters().isEmpty()`, `flag + vc.isFiltered()`,
`flag + vc.getAlternateAllele(0).length() > 0`, `flag + NALLELES > 100` and the
three unflagged controls.  One probe (`N_ALLELES > 100`) could not be compared
(native rejects the JEXL name `N_ALLELES`), and the only mismatch was
`START < 3` (§8.4, a native-only alias that also diverges without the flag).

Edge-shape smoke run (same jar, flag on, `QUAL < 60` and `DP > 5`): a gVCF
`<NON_REF>` block, a 3-ALT record, a deletion, a spanning-deletion `*` and a
biallelic record all match exactly (`LowQual`, `LowQual|LowQual|LowQual`, `SITE`,
`SITE|SITE|SITE`), with no crash in either tool.  An **unflagged** run of the
same fixture exposed one further pre-existing divergence, unrelated to the flag:
for a record with `QUAL=.`, GATK's `QUAL` is still `-10 * NO_LOG10_PERROR` so
`QUAL < 60` fires, while native's `read_field` reports a missing QUAL and the
rule does not fire.  Measured, **not** fixed (site path, out of scope).

## 2. GATK rule, with file:line

* `gatk-source/.../walkers/filters/VariantFiltration.java:359-365` — `apply()`
  takes the allele path whenever `applyForAllele` (`:240`, the
  `--apply-allele-specific-filters` flag) is set. **Not** conditioned on any
  `AS_*` expression:
  ```java
  if (applyForAllele) {
      final List<VariantContext> filtered = splitMultiAllelics(variant).stream()
          .map(vc -> filter(vc, new FeatureContext(featureContext,
                     new SimpleInterval(vc.getContig(), vc.getStart(), vc.getEnd()))))
          .collect(Collectors.toList());
      final List<Set<String>> alleleFilters = filtered.stream()
          .map(filteredvc -> filteredvc.getFilters()).collect(Collectors.toList());
      final VariantContext filteredVC = AlleleFilterUtils.addAlleleAndSiteFilters(
          variant, alleleFilters, invalidatePreviousFilters);
      writer.add(filteredVC);
  } else {
      writer.add(filter(variant, featureContext));
  }
  ```
  That is *why* `AS_FilterStatus` appears even when no allele-level rule matched:
  the flag alone forces the allele path, and `AlleleFilterUtils.java:104-106`
  then seeds the placeholder unconditionally.
* `VariantFiltration.java:371-378` — `splitMultiAllelics` rebuilds each ALT with
  `new VariantContextBuilder("SimpleSplit", contig, start, end, [ref, NO_CALL])`.
  That 5-argument htsjdk constructor sets `attributes = Collections.emptyMap()`
  and `genotypes = GenotypesContext.NO_GENOTYPES` and leaves `log10PError` at the
  builder default, verified in the pinned jar with `javap -constants`:
  `htsjdk.variant.variantcontext.VariantContext.NO_LOG10_PERROR = 1.0d`.
* htsjdk `VariantJEXLContext` (decompiled from the pinned jar with
  `javap -p -c`) maps the JEXL names `CHROM/POS/TYPE/QUAL/ALLELES/N_ALLELES/
  FILTER` to getters and `QUAL` is `lambda$static$2` =
  `-10.0 * vc.getLog10PError()` → `-10.0` in every split context.  Every other
  name is resolved by `vc.hasAttribute(name) ? vc.getAttribute(name)` and by
  filter-name membership; an unknown name yields JEXL `null`.
* `VariantFiltration.java:400-438` — `filter()` adds the mask
  (`:379-392`, `maskVariants.isEmpty() == filterRecordsNotInMask`), runs
  `areClusteredSNPs` (`:415-417`, `CLUSTERED_SNP_FILTER_NAME` at `:530`) and then
  the expressions (`:419-424`, `matchesFilter(vc, null, exp, invert)`), i.e. all
  of it per split ALT.
* `.../mutect/filtering/AlleleFilterUtils.java:94-122` —
  `addAlleleAndSiteFilters` decodes the record's own `AS_FilterStatus`
  (`decodeASFilters`, `:24-28`), returns the record untouched on an arity
  mismatch (`:99-102`), seeds every ALT with `SITE` (`:104-106`,
  `GATKVCFConstants.SITE_LEVEL_FILTERS` at `GATKVCFConstants.java:201`), merges
  the labels through `addAlleleFilters` (`:69-82`), joins with
  `AnnotationUtils.encodeAnyASListWithRawDelim` (`|`,
  `AnnotationUtils.java:21`) and sets the site FILTER to the **intersection**
  (`:115-117`), adding it to the record's existing filters and falling back to
  PASS/`.` (`:118-120`).

## 3. Strict gate (STEP 2)

New file:
`fastgatk-native/scripts/verify_variant_filtration_flag_only_gatk_oracle.py`
(`main()` → int, argparse `--native` /
`FASTGATK_VARIANT_FILTRATION_BINARY` env default, `TemporaryDirectory` scratch,
`oracle_guard` + `FASTGATK_REQUIRE_GATK_ORACLE`, final JSON status line,
`--expect-divergence`), modelled on
`verify_variant_filtration_asfilterstatus_gatk_oracle.py`.

It builds the fixtures in a temp dir, indexes the mask with GATK's own
`IndexFeatureFile` (GATK needs random access: *"Input mask.vcf must support random
access to enable queries by interval"*), runs pinned GATK and native with
identical arguments, and gates `CHROM..FILTER` **and the whole INFO map** (no
floating-point INFO is used, so an exact map comparison is available; key order
is reported, not gated) plus the literal GATK expectation measured in §1.
Gated cases: `flag-qual-expression`, `flag-info-expression`, `flag-mask`,
`no-flag-control`.  Reported only: `flag-without-expression` (§8.1).

Exit status before the fix (recorded literally):

```
[flag-qual-expression] row 0: INFO differ: GATK={'AS_FilterStatus': 'LowQual|LowQual', 'DP': '10'} NATIVE={'DP': '10'}
[flag-qual-expression] row 1: CHROM..FILTER differ: GATK=[... 'LowQual;foo'] NATIVE=[... 'foo']
[flag-qual-expression] row 1: INFO differ: GATK={'AS_FilterStatus': 'LowQual|LowQual', 'DP': '10'} NATIVE={'DP': '10'}
[flag-qual-expression] row 2: INFO differ: GATK={'AS_FilterStatus': 'LowQual', 'DP': '10'} NATIVE={'DP': '10'}
[flag-info-expression] row 0: CHROM..FILTER differ: GATK=[... 'PASS']   NATIVE=[... 'HighDP']
[flag-info-expression] row 0: INFO differ: GATK={'AS_FilterStatus': 'SITE|SITE', 'DP': '10'} NATIVE={'DP': '10'}
[flag-info-expression] row 1: CHROM..FILTER differ: GATK=[... 'foo']    NATIVE=[... 'HighDP;foo']
[flag-info-expression] row 1: INFO differ: GATK={'AS_FilterStatus': 'SITE|SITE', 'DP': '10'} NATIVE={'DP': '10'}
[flag-info-expression] row 2: CHROM..FILTER differ: GATK=[... 'PASS']   NATIVE=[... 'HighDP']
[flag-info-expression] row 2: INFO differ: GATK={'AS_FilterStatus': 'SITE', 'DP': '10'} NATIVE={'DP': '10'}
[flag-mask] row 0: INFO differ: GATK={'AS_FilterStatus': 'Mask|Mask', 'DP': '10'} NATIVE={'DP': '10'}
[flag-mask] row 1: INFO differ: GATK={'AS_FilterStatus': 'SITE|SITE', 'DP': '10'} NATIVE={'DP': '10'}
[flag-mask] row 2: INFO differ: GATK={'AS_FilterStatus': 'SITE', 'DP': '10'} NATIVE={'DP': '10'}
{"status": "divergence", ...}
EXIT=1
```

(`no-flag-control` passed before the fix — the control is the contract that must
not move.)

Exit status after the fix: `EXIT=0`, `"status": "pass"`, `"violations": []`.

## 4. The fix (STEP 3)

`fastgatk-native/src/variant_filtration_tool.cpp` only (+258/−154 including
comments, mostly the re-indented `else` branch).  No header, no shared kernel, no
other tool, no CMake change.

* `apply_allele_specific_filters` (which only ever ran `AS_*` rules) is replaced
  by `apply_allele_path_filters` (`AllelePathOutcome apply_allele_path_filters(...)`,
  defined after `load_cluster_positions`), which reproduces GATK's allele path:
  1. early return when the record has no ALT (`AlleleFilterUtils.java:95-97`);
  2. decode/reseed `AS_FilterStatus` exactly as before (arity guard `:99-102`,
     `SITE` seeds `:104-106`, `--invalidate-previous-filters` reseed `:112-114`);
  3. build **one** scratch record per `bcf_dup(record)` that emulates
     `splitMultiAllelics`'s context: `n_info = n_fmt = n_sample = 0`,
     `d.n_flt = 0`, `qual = kSplitContextQual` (`-10.0F`,
     `VariantContext.NO_LOG10_PERROR`), and, inside the allele loop,
     `n_allele = 2` with `d.allele = [reference, this ALT]`;
  4. per ALT: add the mask name when
     `overlaps_mask(...) != filterRecordsNotInMask` and the cluster name when the
     position is clustered, then evaluate **every** rule against the scratch
     record with the allele index and `missing_fails = false`;
  5. merge the per-ALT sets into the decoded ones with GATK's `addAlleleFilters`
     rules (`:69-82`), encode with `|`/`,` and write `AS_FilterStatus`;
  6. set the site FILTER to the **intersection** of the per-ALT sets, added to the
     record's existing filters (`:115-120`).
* `run_tool`: the per-record body branches on
  `options.apply_allele_specific_filters` — the existing site rule loop + mask +
  cluster blocks are now the `else` branch (byte-identical logic, only
  re-indented), and the flag branch calls the allele path and feeds the
  telemetry (`allele_filtered`, `masked_records`, `clustered_records`,
  `filtered_records`).
* Header: `##INFO=<ID=AS_FilterStatus,…>` is emitted whenever the flag is set
  (previously only when an `AS_*` rule was present); the
  `allele_specific_filters` manifest field is likewise driven by the flag.
* `allele_specific_rules_present` is gone; the `AS_*`-requires-the-flag CLI guard
  is unchanged.
* New file-scope constants `kClusterFilterName` (`VariantFiltration.java:530`)
  and `kSplitContextQual` replace the local `cluster_name` string.

Preserved from the previous round (verified, not assumed): the `SITE`
placeholder, the `|` separator, the decode/merge of a pre-existing
`AS_FilterStatus`, the arity guard, and `missing_fails = false` inside the split
context.

## 5. Test lines corrected (STEP 4)

| file | line(s) | change | why |
| --- | --- | --- | --- |
| `fastgatk-native/scripts/verify_variant_filtration_asfilterstatus_gatk_oracle.py` | 67-71 (docstring), 162-170 (case `invert-filter-expression`) | `"gated": False` → `True`; the stale "native still writes PASS / separately scoped divergence" wording replaced | the site FILTER intersection (`AlleleFilterUtils.java:115-120`) is now implemented, and the case measures identical on both tools; a reported-only case that passes is a weaker gate than a gated one |
| `fastgatk-native/scripts/verify_variant_filtration.py` | 478-496 (comment only) | "AS rules do not alter the site FILTER column" → the measured contract (with the flag, *every* rule is evaluated per ALT and the site FILTER is the intersection) + a reference to the new gate | the old wording pinned the pre-fix model; the assertions themselves (now lines 519-526) were already corrected in the previous round and are unchanged |
| `fastgatk-native/README.md` | 3135-3147 | the `--apply-allele-specific-filters` paragraph now says the flag takes the allele path for **all** expressions, that the per-ALT context has no INFO/genotypes/filters and default QUAL, and names both gates | documentation of the contract, not at the repo root |

No assertion was weakened; the only assertion changes are in the sibling oracle
(`invert-filter-expression` moved from reported to gated, with its measured
expectation already pinned at `("LowASQD", "LowASQD|LowASQD")`).

Other `verify_variant_filtration_*_oracle.py` siblings
(`_gatk_oracle`, `_missing_boolean_gatk_oracle`, `_set_nocall_gatk_oracle`) never
pass `--apply-allele-specific-filters` (grep over `scripts/`), so no assertion of
theirs depended on the old behaviour; all of them pass unchanged (§6, 5b).

## 6. STEP 5 results

| step | command | result |
| --- | --- | --- |
| 5a | `python3 fastgatk-native/scripts/verify_variant_filtration_flag_only_gatk_oracle.py` | **exit 0** (`status: pass`, `violations: []`, 4 gated cases + 1 reported-only) — was `exit 1` before the fix; log `.diag/gate-5a-flag-only.log` (its single printed `VIOLATION:` line belongs to the non-gated `flag-without-expression` case, §8.1) |
| 5b | `python3 fastgatk-native/scripts/verify_variant_filtration.py` | **exit 0** (`{"status": "pass", "filtered_records": 2, "gatk_oracle_exact": true}`) |
| 5b | `ctest --test-dir fastgatk-native/build -R 'filtration' -V` | **exit 0**, 5/5 passed (203.26 s): asfilterstatus-gatk-oracle, contract, gatk-oracle, missing-boolean-gatk-oracle, set-nocall-gatk-oracle — log `.diag/ctest-5b.log` |
| 5c | `ctest --test-dir fastgatk-native/build -R 'window-invariance-gatk-oracle\|…\|asfilterstatus-gatk-oracle' -V` | **exit 0**, 16/16 passed (997.87 s) — log `.diag/ctest-5c.log` |
| 5d | rebuild `fastgatk-variant-filtration` in `fastgatk-native/build` **and** `fastgatk-native/build-serial` | both linked clean, exit 0 (no warnings); both binaries newer than the newest `*.cpp`/`*.hpp` (no staleness) |
| 5d | `fastgatk-native/scripts/run_regression.sh --label asfilterstatus-flag-only` | **exit 0 — OpenMP `build` 297/297 (1310.66 s), Serial `build-serial` 297/297 (1245.97 s)**, zero staleness warnings. Evidence: `.diag/regression/20260911-195723/{omp,serial}.log`, `summary.txt` |

## 7. Tree state

The change **is in the working tree** (no commit, no branch):

```
 M fastgatk-native/README.md
 M fastgatk-native/scripts/verify_variant_filtration.py
 M fastgatk-native/scripts/verify_variant_filtration_asfilterstatus_gatk_oracle.py
 M fastgatk-native/src/variant_filtration_tool.cpp
?? fastgatk-native/scripts/verify_variant_filtration_flag_only_gatk_oracle.py
```

(That is the `未提交变更 5 项` the regression runner recorded for
`8d54a80`; nothing was committed and no branch was created.)

`fastgatk-native/CMakeLists.txt` was not edited (the new gate is run directly, as
the orchestrator registers gates).

## 8. What remains unproven / deliberately not fixed

1. **`--apply-allele-specific-filters` with no expression and no mask** is
   accepted by GATK (case E above, exit 0, `AS_FilterStatus=SITE|SITE`) and
   rejected by native during argument validation
   (`variant_filtration_tool.cpp` parse guard: *"at least one site/genotype filter
   expression or --mask is required"*, exit 2).  That is a CLI-contract
   difference, not an encoding difference, so it is **reported, not gated** in the
   new gate.  Localisation: the `parse()` validation; recommendation: relax the
   guard only together with a dedicated CLI oracle, since it also decides what an
   empty rule list means for the manifest.
2. **Mask/cluster interval fidelity inside the split context.** GATK queries the
   mask with the *trimmed* split interval
   (`new SimpleInterval(vc.getContig(), vc.getStart(), vc.getEnd())`); native
   reuses the record-level overlap query (`overlaps_mask`) for every ALT.  Equal
   for SNVs (measured: case C and the `--filter-not-in-mask` inverse both match);
   **unproven** for indels/symbolic records where `trimAlleles` shortens the span.
   Likewise the split context's mask/cluster labels are not visible to
   FILTER-set-dependent expressions (see 3), while GATK's are.
3. **Bare-identifier vs null-returning-method JEXL semantics.** Measured in the
   split context: `DP > 5` with `--missing-values-evaluate-as-failing` *fires* in
   GATK (the absent identifier makes the expression null, so the policy applies)
   while `vc.getAttribute("DP") > 5` does *not* (a null-valued method call
   compares as `false`).  Native's parser canonicalises both spellings into the
   same node, so the allele path keeps the previous round's choice
   (`missing_fails = false`, which is what the registered
   `missing-values-evaluate-as-failing` case pins).  Result: bare-INFO + that flag
   + this flag is a **known remaining divergence** in the flag path; it is not in
   either gate.  A faithful fix needs a parser-level marker and a per-leaf
   missing-value rule — deliberately not attempted here.
4. **Pre-existing expression-language gaps, unrelated to this change** (they
   differ with *and* without the flag): native treats `START` as an alias for the
   1-based position while GATK's JEXL has no `START` attribute (it has `POS` and
   `vc.getStart()`), and GATK's `CHROM`/`FILTER` attributes are not implemented in
   native.  Exposed by the same probes; left alone.
5. **Not gated**: `--cluster-size`/`--cluster-window-size` + flag
   (GATK needs a random-access input to cluster; not measured at all),
   symbolic/spanning deletions + flag (smoke-tested in §1.1 and matching, but not
   part of the gate fixture), `--set-filtered-genotype-to-no-call` + flag (GATK
   applies the genotype rewrite to the genotype-less split contexts, i.e. it is a
   no-op there, while native still applies genotype rules to the record —
   unchanged by this round), `AS_FilterStatus` arity-mismatch inputs, INFO keys
   other than `AS_QD`, and multiple `--filter-expression`s against a mask.
6. **Speculation (labelled)**: that GATK's INFO-less, genotype-less split context
   is an unintended consequence of the 5-argument `VariantContextBuilder` rather
   than a design choice.  This round matched the measured behaviour
   bug-for-bug; no intent claim is made.  The gate's mask case is SNV-only, so the
   "trimmed interval" question in 2 stays open by construction.
