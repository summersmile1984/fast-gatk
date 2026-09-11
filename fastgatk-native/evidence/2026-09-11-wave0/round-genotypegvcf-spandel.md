# Round: GenotypeGVCFs orphan spanning-deletion fixture (`*` in ALT)

Scope: `fastgatk-native/scripts/verify_genotype_gvcf.py:391`, which asserted
`len(star_records) == 1 and "*" in ALT` for the fixture built at lines 373-381.

## 0. Verdict

* The audit's report is **confirmed by measurement**: pinned GATK 4.6.2.0 emits
  **zero** records for that fixture by default, and exactly one REF-only no-call
  row with `ALT='.'` / `GT='./.'` under `-all-sites`
  (`--include-non-variant-sites`).
* The native **GATK-compatible profile already matched GATK for the default
  options** (zero rows). Only the stale `verify_genotype_gvcf.py` assertion,
  which pinned the *native diagnostic* profile (no
  `--gatk-compatible-annotations`), disagreed.
* The `-all-sites` option **did** diverge: native **crashed** (`exit 2`,
  `invalid genotype PL remap dimensions`) where GATK emits the REF-only row.
  This was a **general** defect, not `*`-specific: it also fires for a plain
  `A G,<NON_REF>` gVCF record whose only ALT fails the standard-confidence
  threshold (measured, see §1.3).
* Fix applied and kept in the tree: a bounded, single-file change in
  `fastgatk-native/src/genotype_gvcf_tool.cpp` that materializes GATK's
  monomorphic REF-only call instead of projecting a 1-allele output through the
  Number=G remap kernel.
* New strict oracle added; it fails before the fix (exit 1) and passes after
  (exit 0).

## 1. STEP 1 — measured GATK truth

### 1.1 Fixture (byte-identical to `verify_genotype_gvcf.py`)

Header = `multi_header` of `verify_genotype_gvcf.py:328-337` with the sample
column renamed to `STAR`; body = line 378-381:

```
chr1	2	.	A	*,G,<NON_REF>	.	PASS	DP=20	GT:DP:AD:PL	0/1:20:12,8,0,0:0,0,100,100,100,100,100,100,100,100
```

Input is a **plain** `.vcf` (as instructed, GATK refuses an unindexed `.vcf.gz`);
it is indexed with `gatk IndexFeatureFile`, which writes a tribble `.idx`
(GATK's `--include-non-variant-sites` group-by-locus traversal requires it).
The reference is a generated 100 bp `chr1` (`>chr1` + `AAAA…`) with a hand-written
`.fai` and a `CreateSequenceDictionary` / hand-written `.dict`; the fixture
header's `##contig=<ID=chr1,length=100>` matches it.

Command (GATK side):

```
third_party/jdk17/bin/java -Xmx1g -jar \
  third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
  GenotypeGVCFs -R <chr1.fa> -V star.g.vcf -O out.vcf \
  --create-output-variant-index false [OPTIONS]
```

### 1.2 Literal rows per option (measured, not assumed)

| OPTIONS | GATK 4.6.2.0 data rows |
| --- | --- |
| *(default)* | **0 rows** |
| `-all-sites` | 1 row: `chr1 2 . A . 127.78 . DP=20;MLEAC=.;MLEAF=. GT ./.` |
| `--include-non-variant-sites` | identical to `-all-sites` |
| `-all-sites --include-non-variant-sites` | identical |
| `-stand-call-conf 20` (default conf) | **0 rows** |
| `-stand-call-conf 10/20/50 -all-sites` | the same single row (QUAL 127.78) |
| `--standard-min-confidence-threshold-for-calling 50 --include-non-variant-sites` | the same single row |

`-all-sites` and `--include-non-variant-sites` are two spellings of one
`@Argument` (`GenotypeGVCFs.java:116-117` and `:126-127`), and the
standard-confidence value cannot change the outcome here because the ALT set is
already empty before any per-allele threshold is consulted.

Literal row, untabulated for readability — the real bytes are tab-separated:

```
chr1	2	.	A	.	127.78	.	DP=20;MLEAC=.;MLEAF=.	GT	./.
```

Note the four properties the old assertion got wrong: `ALT` is `.` (no `*`),
`GT` is `./.`, `FORMAT` is `GT` alone, and by default the record does not exist
at all.

### 1.3 Control fixtures (same measurement round)

| Fixture ALT | OPTIONS | GATK row |
| --- | --- | --- |
| `*,<NON_REF>` (PL `0,0,100,100,100,100`) | default | 0 rows |
| `*,<NON_REF>` | `-all-sites` | `chr1 2 . A . Infinity . DP=20;MLEAC=.;MLEAF=. GT ./.` |
| `G,<NON_REF>` (PL `0,0,100,100,100,100`) | default | 0 rows |
| `G,<NON_REF>` | `-all-sites` | `chr1 2 . A . 23.14 LowQual DP=20;MLEAC=.;MLEAF=. GT ./.` |
| `*,G,<NON_REF>`, G best (PL idx4 = 0) | default | `chr1 2 . A G 82.26 . AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11 GT:AD:DP:PL 0/1:0,20:20:0,0,0` |

The third row-pair is the important one for scope: a fixture with **no `*` at
all** collapses to the same REF-only output, so the native crash was never a
spanning-deletion bug.

## 2. The GATK rule (read from `gatk-source/`, then confirmed by measurement)

1. `GenotypingEngine.calculateOutputAlleleSubset()` publishes an ALT only when
   it passes `standardConfidenceForCalling`, and drops a symbolic `*` outright
   when no emitted deletion covers the locus — `GenotypingEngine.java:311`
   (`isPlausible`), `:314` (`isSpuriousSpanningDeletion`), `:316` (`toOutput`),
   `:318` (`siteIsMonomorphic`), `:327` (returns the subset).
   For this fixture both `*` and `G` are pruned, so the ALT set is empty.
2. With an empty ALT set the site is "monomorphic", so
   `passesEmitThreshold()` is false even at very high QUAL
   (`GenotypingEngine.java:425-427`: `(outputMode == EMIT_ALL_CONFIDENT_SITES || !bestGuessIsRef) && passesCallThreshold`),
   and `calculateGenotypes()` returns `null` — `GenotypingEngine.java:167-169`.
   That is why the **default** traversal writes nothing.
3. `--include-non-variant-sites` / `-all-sites` sets
   `OutputMode.EMIT_ALL_ACTIVE_SITES` (`GenotypeGVCFsEngine.java:381-383`), so
   `emitAllActiveSites()` is true (`GenotypingEngine.java:421-423`) and both
   `null` returns are skipped (`:167-169` and the spanning-deletion-only guard
   at `:173-175`).
4. The rebuilt VariantContext is REF-only:
   `GenotypingEngine.java:190` (`outputAlleles.size() == 1` →
   `GATKVariantContextUtils.subsetToRefOnly(vc, defaultPloidy)`), then
   `GenotypeGVCFsEngine.regenotypeVC()` calls
   `cleanupGenotypeAnnotations(result, createRefGTs=true, keepSB=false)`
   (`GenotypeGVCFsEngine.java:191-194`; defined at `:440`). That helper keeps a
   hom-ref call only when the projected genotype still has `depth > 0` **and**
   a GQ — moving the GQ to RGQ (`:474-484`) — and otherwise publishes `./.`
   with `noGQ()/noDP()` and `noPL()` (`:485-491`).
5. QUAL for a monomorphic site uses the *complementary* posterior:
   `log10Confidence = ... : AFresult.log10ProbVariantPresent()`
   (`GenotypingEngine.java:158-160`), i.e.
   `MathUtils.log10OneMinusPow10(log10PosteriorOfNoVariant)`
   (`AFCalculationResult.java:114-120`).
6. `GenotypeGVCFs.apply()` writes the record because it is not
   "spanning-deletion only" — it has no ALT at all:
   `GenotypeGVCFs.java:323-330` and
   `GATKVariantContextUtils.isSpanningDeletionOnly` (`:2089-2091`).

**Which option produces a surviving record:** only the dense mode
(`-all-sites` / `--include-non-variant-sites`, identical spellings). Its ALT is
`.` and its GT is `./.` (or `0/0` with RGQ when the source genotype carried a
positive GQ). The default option produces **no record**.

## 3. STEP 2 — the strict gate

`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`
(modelled on `verify_reblock_gvcf_droplowqual_gatk_oracle.py`: `main() -> int`,
argparse with `--native` / `$FASTGATK_GENOTYPE_BINARY`, `TemporaryDirectory`,
`oracle_guard`, final JSON status line, `--expect-divergence`).

It runs pinned GATK and native with identical arguments on the same generated
reference and the same plain VCF input, and asserts the data rows are
byte-identical. Native additionally receives `--gatk-compatible-annotations`
(the house convention of every registered genotype-gvcf oracle: the default
native profile is an explicitly non-GATK diagnostic output that publishes
RCQ/RCP and skips GATK's AF-based allele pruning).

Cases: `default`, `include-non-variant-sites`,
`stand-call-conf-50-include-non-variant-sites` (all **gated**), plus
`surviving-concrete-alt` (**reported only**, see §7).

Exit status:

* **before the fix: exit 1** (`status: divergence`). Literal evidence:
  * `[default]` GATK 0 rows / native 0 rows — already matching.
  * `[include-non-variant-sites]` `gatk_exit=0 native_exit=2`,
    `native: invalid genotype PL remap dimensions`;
    `VIOLATION: record count differs: GATK=1 native=0`.
  * `[stand-call-conf-50-include-non-variant-sites]` same two violations.
* **after the fix: exit 0** (`status: pass`); all three gated cases report
  `GATK rows == NATIVE rows` byte-for-byte, including
  `chr1	2	.	A	.	127.78	.	DP=20;MLEAC=.;MLEAF=.	GT	./.`.

Evidence logs: `.diag/spandel-gate-before.log`, `.diag/spandel-gate-after.log`.

## 4. STEP 3 — the fix

Localization (before the fix):
`apply_gatk_output_allele_subset()` dropped a REF-only output only when
`!include_non_variant_sites`; in dense mode it fell through to
`remap_record_to_allele_union(..., output_alleles)` with a **1-allele** target,
and `fastgatk-kernels/src/genotype.cpp:795-800`
(`remap_genotype_pl_kokkos`) rejects `target_allele_count < 2` with
`invalid genotype PL remap dimensions`. Backtrace confirmed with gdb:
`remap_genotype_pl_kokkos` ← `remove_non_ref_allele` ← `run_tool`
(`genotype_gvcf_tool.cpp:2606-2617` region, kernels/genotype.cpp:800).

Change (all in `fastgatk-native/src/genotype_gvcf_tool.cpp`, +210/-11):

1. `Record`: added `source_has_gq`, `finalized_monomorphic_ref`,
   `cohort_log10_p_no_variant`, `cohort_quality_available`.
2. `extract_materialized_fields()`: records whether the **input** FORMAT
   carried GQ (native later synthesizes a PL-derived GQ, which GATK does not
   have — this is what selects GATK's `hasGQ()` branch).
3. `merge_sample_fields()`: OR the flag across the merged group.
4. `update_cohort_af_annotations()`: keeps
   `AlleleFrequencyResult.log10_p_no_variant` on the record.
5. New `materialize_gatk_monomorphic_ref_call()`: builds GATK's REF-only row —
   `ALT='.'`, per-sample GT `./.` (or `0/0` + `RGQ` when the source had a
   positive GQ and depth), FORMAT reduced to `GT` (or `GT:DP:RGQ`), INFO
   `MLEAC`/`MLEAF` written as missing Number=A values, INFO `DP` preserved, and
   QUAL recomputed with `gatk_log10_one_minus_pow10()` (the monomorphic
   formula).
6. `apply_gatk_output_allele_subset()`: the `output_alleles.size() == 1` branch
   now materializes the REF-only call in dense mode instead of falling through
   to the 2-allele-only remap.
7. Both compute lambdas (stream and aggregate) short-circuit the remaining
   annotation stages when `record.finalized_monomorphic_ref` is set — GATK
   emits this row straight out of `regenotypeVC` with no site QD/FS/SOR,
   no AC/AF/AN, and no ExcessHet.

Not restructured: the genotyping engine, the AF kernel, the PL remap kernel and
all other tools are untouched. The change is confined to one translation unit.

The full diff is preserved in `.diag/spandel-fix.diff`.

## 5. STEP 4 — corrected test lines

`fastgatk-native/scripts/verify_genotype_gvcf.py` (old line numbers refer to the
pre-change file):

| Old line(s) | Change |
| --- | --- |
| 369-372 (comment) | replaced with the measured GATK rule and `file:line` citations, naming the option each expectation belongs to |
| 383-386 | the star run now passes `--gatk-compatible-annotations`, so it exercises the GATK-comparable profile |
| **391** `assert len(star_records) == 1 and "*" in star_records[0][4]` | **replaced** by `assert star_records == []` — the GATK contract for the **default** option |
| 394-395 (new) | added a `--include-non-variant-sites --gatk-compatible-annotations` run |
| 416-417 (new) | `assert star_nonvariant_rows == ["chr1\t2\t.\tA\t.\t127.78\t.\tDP=20;MLEAC=.;MLEAF=.\tGT\t./."]` — exact literal row for the **`-all-sites`** option |
| 396-399 | the `--gp-qual` run now also passes `--include-non-variant-sites --gatk-compatible-annotations` |
| 405 (old) `float(star_posterior_records[0][5]) <= float(star_records[0][5])` | compared against the non-variant row instead of the now-empty default run |

Line numbers in the new file: comment block 369-380, compat flag at 392,
`assert star_records == []` at 400, the dense-mode run at 407-413, the exact row
assertion at 416-417, the `--gp-qual` dense run at 419-424, the QUAL comparison
at 433-435.

`verify_genotype_gvcf.py` now passes (`{"status": "pass", "output_records": 1}`,
exit 0).

### Sibling `verify_genotype_gvcf_*.py` files

Checked each registered sibling for assertions that depend on the old
`*`-survives behaviour; none needed a change, and all pass (see §6):

* `verify_genotype_gvcf_spanning_deletion_gatk_oracle.py` — its fixtures use
  `C,*,<NON_REF>` where `C` is plausible or owned by a deletion, so GATK keeps
  `C` and the native contract there is unchanged (`orphan_star_removed`,
  `supported_star_preserved` both verified against GATK).
* `verify_genotype_gvcf_max_alternate_alleles_gatk_oracle.py`,
  `..._multisample_reference_confidence_oracle.py`,
  `..._assignment_...`, `..._exclude_intervals_...`, `..._gp_input_...`,
  `..._include_non_variant_...`, `..._starts_in_intervals_...`,
  `..._multisample_...`, `..._multiallelic_...`, `..._inbreeding_...`,
  `..._legacy_qual_...` — no `*`-in-ALT expectation and no dependence on the
  collapsed-allele-set path; all pass unchanged.

The non-`_gatk_oracle` files (`verify_gatk_genotype_gvcf*.py`) compare GATK
header/annotation sets for HC-produced gVCFs and do not touch this fixture.

## 6. STEP 5 — gate results

| Step | Command | Result |
| --- | --- | --- |
| a | `python3 fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` | **exit 0** (exit 1 before the fix) |
| b1 | `python3 fastgatk-native/scripts/verify_genotype_gvcf.py` | **exit 0** |
| b2 | `ctest --test-dir fastgatk-native/build -R 'genotype-gvcf' -V` | **exit 0 — 14/14 passed** |
| c | the 16-name strict-gate filter from the task | **exit 0 — 18/18 matched tests passed** |
| d | `run_regression.sh --label genotype-gvcf-spandel` | see §6.1 |

Details of b2: contract, gatk-oracle, legacy-qual, multisample, multiallelic,
inbreeding, include-non-variant, assignment, multisample-reference-confidence,
max-alternate-alleles, gp-input, exclude-intervals, starts-in-intervals,
spanning-deletion.

Details of c: window-invariance, ploidy-window-invariance, alleles-overlap,
span-del-qual, gvcf-symbolic-prior, arbitrary-ploidy-span-del-prior,
polyploid-gvcf-span-del-prior, spanning-prior-genotype-gq, af-zero-format,
multialt-owner-annotation, mutect2-recheck ×2, gvcf-indel-end, culprit,
select-variants-refonly, asfilterstatus, variant-filtration-flag-only,
droplowqual.

### 6.1 Full double-backend regression

`fastgatk-native/scripts/run_regression.sh --label genotype-gvcf-spandel`
(omp + serial, `FASTGATK_REQUIRE_GATK_ORACLE=1` by default).

```
| 后端   | 构建目录                        | 结果 | 通过/总数 | 耗时        |
| omp    | OpenMP (fastgatk-native/build)  | 通过 | 299/299   | 1288.34 sec |
| serial | Serial (fastgatk-native/build-serial) | 通过 | 299/299 | 1327.48 sec |
```

`REG_EXIT=0`. Evidence block:
`.diag/regression/20260911-231503/{evidence.md,omp.log,serial.log,omp.status,serial.status}`.
The runner's staleness check reported no warning for either build directory.

Caveat, stated plainly: during the *comment-only* line-reference corrections
(the `file:line` citations inside comments, applied while the suite was
already running) the source mtime moved ahead of the binaries. No statement,
type or control flow changed. Both trees were rebuilt afterwards and re-run:
the new oracle gate is exit 0 again, and the `genotype-gvcf` ctest subset is
14/14 again on the rebuilt binaries, which is the state the tree now holds.
The 299/299 numbers above come from the run on the pre-correction binaries,
whose object code is identical.

## 7. What remains unproven / open

* **`surviving-concrete-alt` (reported only, still divergent).** With
  `*,G,<NON_REF>` and G owning the best genotype, GATK emits
  `A G 82.26 . AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11 GT:AD:DP:PL 0/1:0,20:20:0,0,0`;
  native's compat profile emits the same ALT/QUAL/FILTER but
  `AC=0;AF=0.00;AN=0;…;GT:AD:DP:GQ:PL ./.:0,20:20:0:0,0,0`. That is a
  **PREFER_PLS genotype re-derivation** difference after an orphan `*` is
  pruned (GATK's projected PL row makes it call `0/1`; native derives a
  hom-ref tie, zero GQ, and its existing "uninformative hom-ref → no-call"
  rule turns it into `./.`). It is *not* part of this round's fixture and is
  deliberately left unfixed and un-gated.
* **`FILTER=LowQual` is not implemented in native GenotypeGVCFs at all**
  (measured: `grep -c LowQual fastgatk-native/src/genotype_gvcf_tool.cpp` = 0).
  The measured GATK rule is `QUAL < standardConfidenceForCalling → LowQual`
  (`GenotypingEngine.java:183-186`). The gated fixtures were chosen so that
  QUAL is 127.78 ≥ 30 and the filter is `.`; the collapsed case with a low
  QUAL (`A G,<NON_REF> -all-sites` → `23.14 LowQual`) therefore still differs
  in the FILTER column. Speculation: adding it globally could shift other
  registered expectations, so it was not attempted inside this round.
* **`QUAL=Infinity`** (the `*,<NON_REF>`-only fixture) is rendered by htslib
  as `inf`/`0` rather than GATK's literal `Infinity`
  (`gatk_qual_output` maps non-finite to 0). Not gated; separate formatting
  concern.
* **`GT:DP:RGQ` branch.** The GQ-present variant of
  `cleanupGenotypeAnnotations` is implemented from source but is **not
  measured** — none of the fixtures in this round carried FORMAT/GQ on a
  record whose ALTs all collapse. The registered
  `include-non-variant` oracle (HC-produced gVCFs, which do carry GQ) passes,
  but it is not known to exercise this exact branch.
* **`source_has_gq` is a per-record flag, not per-sample.** GATK tests
  `g.hasGQ()` per genotype; multi-shard inputs whose shards differ in FORMAT
  could therefore differ. No fixture or registered gate exercises that, and the
  affected code path previously crashed, so this cannot silently change an
  existing result — but it is a deliberate simplification.
* **Other INFO fields** are carried through untouched by
  `materialize_gatk_monomorphic_ref_call()`. Only `DP` was measured here; an
  input gVCF carrying other INFO keys at a collapsing locus has not been
  checked against GATK.
* The `-all-sites` **short name** is not accepted by native
  (`unknown option: -all-sites`); the gate therefore uses the long name on both
  sides. Adding the alias is a CLI-surface change outside this round's scope.

## 8. Tree state

* Change **kept in the tree** (no commit, no branch).
* Modified: `fastgatk-native/src/genotype_gvcf_tool.cpp`,
  `fastgatk-native/scripts/verify_genotype_gvcf.py`.
* Added: `fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`.
* Not touched: `fastgatk-native/CMakeLists.txt`, root `*.md`, Mutect2 and every
  other tool.
* Both build trees (`fastgatk-native/build`, `fastgatk-native/build-serial`)
  rebuilt with the change.
