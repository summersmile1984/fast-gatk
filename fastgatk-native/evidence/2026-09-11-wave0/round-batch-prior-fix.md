# Round: batch fix of the three remaining `*` (spanning-deletion) prior sites

Date: 2026-09-11 (local). Repo `/home/turing-agents/Documents/fast-gatk`.
Base revision `d2fe345` (working tree clean at start). **No commit, no branch**; the three fixes are
left in the working tree. `fastgatk-native/CMakeLists.txt` was **not** touched, no root `*.md` was
touched, no oracle was edited or weakened, Mutect2 untouched, no scratch outside
`tempfile.TemporaryDirectory`.

## 0. Bottom line

Shared root rule (GATK `AlleleFrequencyCalculator.java:175-176`, mirrored for priors by
`GenotypePriorCalculator.java:139-153`):

```
a.length() == refLength ? snpPseudocount : indelPseudocount     // refLength = vc.getReference().length()
```

`htsjdk Allele.SPAN_DEL` is `create("*", false)`: `isCalled() == true`, `isSymbolic() == false`,
`length() == 1`. So on a **1 bp REF** record `*` is an SNP-class allele (SNP pseudocount;
`log10(snpHet) - log10(3)` for the genotype prior). Only a REF longer than 1 bp gives it the
indel class.

| # | Target | Status | Oracle before | Oracle after |
|---|---|---|---|---|
| 1 | `hc_call.cpp` ordinary-VCF arbitrary-ploidy `*` AF pseudocount | **FIXED, CONTAINED** | exit **1** | exit **0** (`status: pass`) |
| 2 | `hc_call.cpp` polyploid / hidden-spanning gVCF `*` AF pseudocount | **FIXED, CONTAINED** | exit **1** | exit **0** (`status: pass`) |
| 3 | `hc_call.cpp:682` `*` genotype prior (`GenotypePriorCalculator` path) | **FIXED, CONTAINED** | exit **1** | exit **0** (`status: pass`) |
| — | optional: native rejects `--heterozygosity 0` / `--indel-heterozygosity 0` / `--heterozygosity-stdev 0` | **NOT ATTEMPTED** (investigated, see §4) | — | — |

Nothing was reverted: every target's oracle went from red to green on the first application of the
minimal fix, and the full dual-backend suite stayed at 286/286. All three fixes live in one file,
`fastgatk-native/src/hc_call.cpp` (+45/−4 lines, comments included).

Note on the brief's literal values: the brief lists target 1's pre-fix pair as
"368.90 (native) vs 369.37 (GATK)". The observed direction is the opposite — **GATK 368.90 /
native 369.37** — which is confirmed independently by the pre-fix oracle run below and is the
direction the defect predicts (native's `*` gets the *smaller* indel pseudocount, so P(no variant)
and therefore QUAL come out *higher*). The brief's two numbers are the right pair with the two
labels swapped.

---

## 1. Target 1 — ordinary-VCF arbitrary-ploidy `*` AF pseudocount

Site (before): `fastgatk-native/src/hc_call.cpp:3358-3360` in `vcf_text(...)`.
Oracle: `fastgatk-native/scripts/verify_hc_arbitrary_ploidy_span_del_prior_fixture_oracle.py`
(strict mode, all 5 cases: 3 divergence-gated + 2 byte-parity controls).

### Diff applied (now at `hc_call.cpp:3382-3390`)

```diff
             if (has_spanning_deletion)
-                prior_pseudocounts.back() = result.genotype_indel_heterozygosity *
-                    ref_pseudocount;
+                // AlleleFrequencyCalculator selects every allele's Dirichlet
+                // pseudocount from its length:
+                // `a.length() == refLength ? snpPseudocount : indelPseudocount`
+                // (AlleleFrequencyCalculator.java:175-176, refLength =
+                // vc.getReference().length()).  htsjdk's symbolic spanning
+                // deletion is a 1 bp allele (`Allele.SPAN_DEL.length() == 1`),
+                // so on a 1 bp REF record it takes the SNP pseudocount and only
+                // a longer REF gives it the indel pseudocount.  The
+                // unconditional indel pseudocount here under-weighted `*` and
+                // therefore changed P(no variant) = the record's QUAL.
+                prior_pseudocounts.back() =
+                    (ref.size() == 1 ? result.genotype_snp_heterozygosity
+                                     : result.genotype_indel_heterozygosity) *
+                    ref_pseudocount;
```

`ref` is the record's own REF allele (`std::get<2>(key)`, i.e. `candidate_reference(...)` of the
group key), which is exactly `vc.getReference()` for the emitted record. This is the same ternary
shape as the two previously landed instances (`calling_pipeline.cpp:17320-17324` and
`hc_call.cpp:4876-4878`), which are the only other `*` AF-prior sites that exist.

### Oracle exit status before / after

```
timeout 2400 python3 fastgatk-native/scripts/verify_hc_arbitrary_ploidy_span_del_prior_fixture_oracle.py
# BEFORE: exit 1   (.diag/batch_t1_prefix.log)
# AFTER : exit 0   (.diag/batch_t1_final.log)  "status": "pass"
```

### Literal GATK vs native values

BEFORE (`violations`, `.diag/batch_t1_prefix.log`):

```
ord-ploidy3-cigar-del3: record 0 (chr1:698): QUAL '368.90' -> '369.37'
ord-ploidy3-cigar-del3: record 0 (chr1:698): INFO QD '8.20' -> '8.21'
ord-ploidy4-cigar-del2: record 0 (chr1:698): QUAL '391.84' -> '391.94'
```

AFTER (byte-identical, `.diag/batch_t1_final.log`):

```
GATK   chr1 698 . C A 368.90 . AC=1;AF=0.333;AN=3;...;QD=8.20;...  GT:AD:DP:GQ:PL  0/0/1:33,12:45:58:424,0,58,1244
NATIVE chr1 698 . C A 368.90 . AC=1;AF=0.333;AN=3;...;QD=8.20;...  GT:AD:DP:GQ:PL  0/0/1:33,12:45:58:424,0,58,1244
GATK   chr1 698 . C A 391.84 . AC=1;AF=0.250;AN=4;...;QD=8.91;...  GT:AD:DP:GQ:PL  0/0/0/1:32,12:44:19:426,0,19,90,1241
NATIVE chr1 698 . C A 391.84 . AC=1;AF=0.250;AN=4;...;QD=8.91;...  GT:AD:DP:GQ:PL  0/0/0/1:32,12:44:19:426,0,19,90,1241
```

Controls unchanged and byte-identical: `ord-ploidy3-no-deletion` 421.02/421.02 (no `*` at all),
`ord-ploidy2-cigar-del3` 384.86/384.86 (diploid path), `ord-ploidy1-snp30-del6` 409.04/409.04.

---

## 2. Target 2 — polyploid / hidden-spanning gVCF `*` AF pseudocount

Site (before): `fastgatk-native/src/hc_call.cpp:4479-4483` in `gvcf(...)`.
Oracle: `fastgatk-native/scripts/verify_hc_polyploid_gvcf_span_del_prior_fixture_oracle.py`
(strict mode; hard-gated cases `bp-ploidy3-cigar-del3`, `bp-ploidy3-frameshift-del3` plus the two
controls; the two `del2` cases are `report_only` probes whose flanking `chr1:697` differences are
notes only).

### Diff applied (now at `hc_call.cpp:4505-4527`)

```diff
                     if (include_spanning_deletion || hidden_spanning_deletion) {
                         const auto spanning = group.candidates.size() + 1U;
+                        // AlleleFrequencyCalculator selects every allele's
+                        // Dirichlet pseudocount from its length:
+                        // `a.length() == refLength ? snpPseudocount :
+                        // indelPseudocount` (AlleleFrequencyCalculator.java:175-176,
+                        // refLength = vc.getReference().length()).  htsjdk's
+                        // symbolic spanning deletion is a 1 bp allele
+                        // (`Allele.SPAN_DEL.length() == 1`), so on a 1 bp REF
+                        // record it takes the SNP pseudocount and only a longer
+                        // REF gives it the indel pseudocount.  The unconditional
+                        // indel pseudocount here under-weighted `*` and moved
+                        // P(no variant), i.e. the reference-confidence QUAL, of
+                        // every polyploid gVCF record that carries `*`.
                         prior_pseudocounts[spanning] =
-                            result.genotype_indel_heterozygosity * ref_pseudocount;
+                            (group.reference.size() == 1
+                                 ? result.genotype_snp_heterozygosity
+                                 : result.genotype_indel_heterozygosity) *
+                            ref_pseudocount;
                     }
```

`group.reference` is the same expression the already-fixed sibling in the diploid gVCF writer uses
(`hc_call.cpp:4876-4878`), so the three gVCF/ordinary writers now agree by construction. The
`<NON_REF>` allele's unconditional indel pseudocounts (`hc_call.cpp:4530-4531` polyploid,
`hc_call.cpp:4883-4884` diploid) are **correct and untouched**: `Allele.NON_REF_ALLELE.length() == 0`,
which is never equal to a real REF length.

### Oracle exit status before / after

```
timeout 2400 python3 fastgatk-native/scripts/verify_hc_polyploid_gvcf_span_del_prior_fixture_oracle.py
# BEFORE: exit 1   (.diag/batch_t2_prefix.log)
# AFTER : exit 0   (.diag/batch_t2_final.log)  "status": "pass"
```

### Literal GATK vs native values (QUAL at `chr1:698`)

| case | ploidy | deletion reads | GATK before | native before | GATK after | native after |
|---|---|---|---|---|---|---|
| `bp-ploidy3-cigar-del3` (gated) | 3 | 3 (real `138M1D161M`) | 368.90 | **369.37** | 368.90 | 368.90 |
| `bp-ploidy3-frameshift-del3` (gated) | 3 | 3 (299M/sequence removal) | 368.90 | **369.37** | 368.90 | 368.90 |
| `bp-ploidy3-cigar-del2` (probe) | 3 | 2 | 411.25 | **416.91** | 411.25 | 411.25 |
| `bp-ploidy4-cigar-del2` (probe) | 4 | 2 | 391.84 | **391.94** | 391.84 | 391.84 |
| `bp-ploidy2-cigar-del3` (control) | 2 | 3 | 384.86 | 384.86 | 384.86 | 384.86 |
| `bp-ploidy3-no-deletion` (control) | 3 | 0 | 421.02 | 421.02 | 421.02 | 421.02 |

BEFORE `violations`:

```
bp-ploidy3-cigar-del3: record 198 (chr1:698): QUAL '368.90' -> '369.37'
bp-ploidy3-frameshift-del3: record 198 (chr1:698): QUAL '368.90' -> '369.37'
```

The report-only probe notes dropped from 13 to 11: the two `chr1:698 QUAL` entries disappeared and
only the pre-existing (unrelated) `chr1:697` flanking-record differences remain — GATK prints a
zero/known QUAL and rank-sums there, native does not retain that site. Line execution is proven
directly by the emitted ALT list: native emits `C A,*,<NON_REF>` at `chr1:698`
(`native_emits_spanning_deletion: True` in all four `*` cases before and after).

---

## 3. Target 3 — `*` genotype prior in `genotype_priors_for_group`

Site (before): `fastgatk-native/src/hc_call.cpp:682-688`. Reached only with
`--genotype-assignment-method USE_POSTERIOR_PROBABILITIES`, because
`result.genotype_priors_used = options.use_genotype_priors && options.use_posterior_genotype_assignment`
and HaplotypeCaller's default `USE_PLS_TO_ASSIGN` never reads the priors.
Oracle: `fastgatk-native/scripts/verify_hc_spanning_prior_genotype_gq_fixture_oracle.py`
(strict mode; 4 divergence-gated cases + 2 controls).

### Diff applied (now at `hc_call.cpp:682-704`)

```diff
     if (include_spanning_deletion) {
-        const auto spanning_heterozygosity = result.genotype_indel_heterozygosity;
+        // GenotypePriorCalculator.assumingHW classifies every called,
+        // non-symbolic allele from its length: `allele.length() == refLength ?
+        // SNP : INDEL` (GenotypePriorCalculator.java:139-153) and both SNP
+        // values carry the log10(3) normalization (het = log10(snpHet) -
+        // log10(3), hom-var = 2*log10(snpHet) - log10(3), i.e. NOT twice the
+        // normalized het value).  htsjdk's symbolic spanning deletion is a 1 bp
+        // allele (`Allele.SPAN_DEL.length() == 1`), so on a 1 bp REF record it
+        // is a SNP event; only a longer REF makes it an INDEL event.
+        const bool spanning_is_snp = !candidates.empty() &&
+            candidate_reference(*candidates.front()).size() == 1;
+        const auto spanning_heterozygosity = spanning_is_snp
+            ? result.genotype_snp_heterozygosity
+            : result.genotype_indel_heterozygosity;
         if (!(spanning_heterozygosity > 0.0) || !std::isfinite(spanning_heterozygosity))
             return {};
         const auto log10_spanning_het = std::log10(std::max(spanning_heterozygosity, 1e-300));
-        allele_priors.emplace_back(log10_spanning_het, 2.0 * log10_spanning_het);
+        const auto spanning_snp_normalization = spanning_is_snp ? std::log10(3.0) : 0.0;
+        allele_priors.emplace_back(log10_spanning_het - spanning_snp_normalization,
+                                   2.0 * log10_spanning_het - spanning_snp_normalization);
     }
```

The `log10(3)` normalization matches the concrete-allele path exactly
(`calling_pipeline.cpp:16621-16626`: `prior_het = log10_het - snp_normalization`,
`prior_hom_alt = 2.0*log10_het - snp_normalization`), i.e. the same shape as
`GenotypePriorCalculator`'s `homValues[SNP] = snpHom - LOG10_SNP_NORMALIZATION_CONSTANT` with
`snpHom = snpHet*2` from `assumingHW`. The `include_non_ref` OTHER branch just below is left
alone (GATK's `hetValues[OTHER] = max(snpHet, indelHet)`, no normalization — already correct).

### Oracle exit status before / after

```
timeout 2400 python3 fastgatk-native/scripts/verify_hc_spanning_prior_genotype_gq_fixture_oracle.py
# BEFORE: exit 1   (.diag/batch_t3_prefix.log)
# AFTER : exit 0   (.diag/batch_t3_final.log)  "status": "pass"
```

### Literal GATK vs native values (GQ at `chr1:698`)

| case | ploidy | del reads | GATK GQ | native GQ before | native GQ after | GT (both sides) |
|---|---|---|---|---|---|---|
| `priors-p3-cigar-del3` | 3 | 3 | 29 | **33** | 29 | `REF/REF/A` (`0/0/2` vs `0/0/1` after ALT reorder) |
| `priors-p3-cigar-del2` | 3 | 2 | 61 | **65** | 61 | `REF/REF/A` |
| `priors-p3-cigar-del6` | 3 | 6 | 67 | **63** | 67 | `0/1/2` on both sides |
| `priors-p4-cigar-del2` | 4 | 2 | 27 | **31** | 27 | `REF/REF/REF/A` |
| `priors-no-spanning-deletion` (control) | 3 | 0 | 84 | 84 | 84 | `REF/REF/A` |
| `priors-disabled-span-del` (control) | 3 | 3, no `USE_POSTERIOR_PROBABILITIES` | 6 | 6 | 6 | `0/1/2` |

BEFORE `violations`:

```
priors-p3-cigar-del3: record 198 (chr1:698): FORMAT GQ '29' -> '33'
priors-p3-cigar-del2: record 198 (chr1:698): FORMAT GQ '61' -> '65'
priors-p3-cigar-del6: record 198 (chr1:698): FORMAT GQ '67' -> '63'
priors-p4-cigar-del2: record 198 (chr1:698): FORMAT GQ '27' -> '31'
```

AFTER: `violations: []` — GQ identical on all six cases, including the non-one-signed `del6` case
(native 63→67). This is the predicted `10*(log10(1e-3/1.25e-4) - log10(3)) = 10*(0.903-0.477) ≈
4.26` phred, seen as ±4 integer GQ in the fixture — the direction flips with the called genotype,
which is what a prior (not a likelihood) defect looks like.

The two unrelated FORMAT gaps this oracle reports as notes are unchanged and still open:
native emits no `GP`/`PG` under `USE_POSTERIOR_PROBABILITIES` (`.diag/round-fixture-construction.md`
§D.2), and the `AF`/`MLEAF` zero-precision note (§D.1) is unchanged (see §5).

---

## 4. Optional item: zero-valued heterozygosity options — investigated, **NOT** changed

Measured directly (scratch inside `TemporaryDirectory`, same fixture family, `-ERC BP_RESOLUTION
--sample-ploidy 3`, `-L chr1:500-780 --min-pruning 1`), GATK 4.6.2.0 vs native:

| option | GATK | native | GATK `chr1:698` record |
|---|---|---|---|
| `--heterozygosity 0` | exit **0**, 281 records | exit **2** `error: calling thresholds must be positive` | `C <NON_REF> . . . GT:AD:DP:GQ:PL 0/0/0:33,12:45:0:0,0,0,940` (both ALTs collapse to `<NON_REF>`) |
| `--indel-heterozygosity 0` | exit **0**, 281 records | exit **2** `error: calling thresholds must be positive` | `C *,A,<NON_REF> 368.90 ...` — **identical to the default run** |
| `--heterozygosity-stdev 0` | exit **0**, 281 records | exit **2** `error: call-confidence prior/threshold options are invalid` | `C <NON_REF> . . . GT:AD:DP:GQ:PL 0/0/0:33,12:45:0:0,0,0,940` |

So GATK genuinely accepts all three (the fixture-construction report's §D.3 is confirmed
empirically). I nevertheless did **not** relax native's validation, because it is neither minimal
nor safely gateable in this round:

1. **`--heterozygosity 0` / `--heterozygosity-stdev 0` are degenerate, not merely accepted.** In
   GATK the SNP heterozygosity is both the SNP prior and the numerator of
   `refPseudocount = snpHet / stdev^2`, and `stdev` is its denominator. Setting either to 0 makes
   GATK drop both alternate alleles at `chr1:698` (record above). Reproducing that output requires
   reproducing `log10(0) = -Infinity` and infinitely-large/zero Dirichlet pseudocount semantics
   through the AF/prior kernels — native currently clamps with `std::max(het, 1e-300)` and skips the
   AF block when `ref_pseudocount` is not finite/positive. That is a pipeline change, not a
   validation relaxation.
2. **It would activate an already-known, unfixed divergence.** With `heterozygosity <= 0` native's
   `ref_pseudocount` guard fails, so the gVCF/ordinary writers fall into the
   `else if (!include_spanning_deletion)` fallback `estimate_mle_allele_counts(..., include_non_ref=true)`
   (`hc_call.cpp:2938-2940`), which uses `max(snpHet, indelHet)` for `<NON_REF>` where GATK uses the
   indel prior — sweep item #11 / fixture-round §C.3, explicitly not fixed. Relaxing the gate would
   start exercising a *known-bad* path.
3. **It cannot be gated the way the other gates are.** Registering a new CTest oracle requires
   editing `fastgatk-native/CMakeLists.txt`, which this round forbids (hard constraint 1), and the
   brief's final gate lists fixed test names. An ungated exit-code change on a shared validation
   would be a strictly worse trade than leaving the divergence documented.

For `--indel-heterozygosity 0` alone the change looks *behaviourally* inert on this fixture (GATK's
output is byte-identical to the default run, and after target 1/2 `*` no longer consumes the indel
pseudocount on a 1 bp REF at all), so a future round could relax just that one value with a fixture
— but that is a separate, separately-gated change. **Speculation:** whether a concrete indel allele
(spliced into an AF matrix) would stay identical under `--indel-heterozygosity 0` is untested; I did
not build that fixture.

---

## 5. Final gate

All commands run against the **final** working tree (all three fixes, git `d2fe345` + 1 modified
file), OpenMP binary `fastgatk-native/build/fastgatk-hc-call` and Serial binary
`fastgatk-native/build-serial/fastgatk-hc-call`, both rebuilt.

### (a) three site oracles (final binary)

| command | exit | status |
|---|---|---|
| `python3 fastgatk-native/scripts/verify_hc_arbitrary_ploidy_span_del_prior_fixture_oracle.py` | **0** | `"status": "pass"` (`.diag/batch_t1_final.log`) |
| `python3 fastgatk-native/scripts/verify_hc_polyploid_gvcf_span_del_prior_fixture_oracle.py` | **0** | `"status": "pass"` (`.diag/batch_t2_final.log`) |
| `python3 fastgatk-native/scripts/verify_hc_spanning_prior_genotype_gq_fixture_oracle.py` | **0** | `"status": "pass"` (`.diag/batch_t3_final.log`) |

### (b) previously registered strict gates

```
timeout 2400 third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest --test-dir fastgatk-native/build \
  -R 'window-invariance-gatk-oracle|alleles-overlap-gate-oracle|span-del-qual-gatk-oracle|gvcf-symbolic-prior-gatk-oracle' -V
```

**exit 0** (`.diag/batch_ctest_registered.log`):

```
1/5 Test #83: fastgatk-hc-window-invariance-gatk-oracle ..........   Passed   60.47 sec
2/5 Test #84: fastgatk-hc-ploidy-window-invariance-gatk-oracle ...   Passed   44.10 sec
3/5 Test #85: fastgatk-hc-alleles-overlap-gate-oracle ............   Passed   53.00 sec
4/5 Test #86: fastgatk-hc-span-del-qual-gatk-oracle ..............   Passed  120.55 sec
5/5 Test #87: fastgatk-hc-gvcf-symbolic-prior-gatk-oracle ........   Passed   77.43 sec
100% tests passed, 0 tests failed out of 5
```

### (c) alleles oracle

```
timeout 3000 python3 fastgatk-native/scripts/verify_hc_alleles_gatk_oracle.py
# exit 0, "status": "pass", "multi_alt_and_anchored_indel_gatk_exact": true   (.diag/batch_alleles_oracle.log)
```

### (d) full dual-backend regression

```
timeout 5400 fastgatk-native/scripts/run_regression.sh --label 'batch symbolic prior fix'
# exit 0
```

Evidence block (`.diag/batch_full_regression.out`, outputs in
`.diag/regression/20260911-081700/{omp,serial}.log`):

```
- 标签：batch symbolic prior fix
- git：`d2fe345`（未提交变更 1 项）
- 过滤：`<全量>`　并行度：8

| 后端 | 构建目录 | 结果 | 通过/总数 | 耗时 |
| --- | --- | --- | --- | --- |
| omp | `OpenMP (fastgatk-native/build)` | 通过 | 286/286 | 1068.20 sec |
| serial | `Serial (fastgatk-native/build-serial)` | 通过 | 286/286 | 1058.98 sec |
```

`100% tests passed, 0 tests failed out of 286` on both backends, including all five registered
strict gates (omp `#83`–`#87`, serial `#83`–`#87`, all `Passed`).

### (e) nothing had to be reverted

No target was reverted. The "if the full suite fails" branch of the brief was not entered, so
§targets-reverted is empty.

---

## 6. Which targets are CONTAINED, and which are NOT

**CONTAINED in the working tree** (present, individually gated, and consistent with the full
286/286 dual-backend run):

* target 1 — arbitrary-ploidy ordinary-VCF `*` AF pseudocount now length-derived
  (`hc_call.cpp:3382-3390`).
* target 2 — polyploid/hidden-spanning gVCF `*` AF pseudocount now length-derived
  (`hc_call.cpp:4505-4527`).
* target 3 — `*` genotype prior is SNP-class with the `log10(3)` normalization on a 1 bp REF, and
  INDEL-class otherwise (`hc_call.cpp:682-704`).

**NOT contained** (still open, unchanged by this round, and none of them is reached by the three
fixtures):

* `calling_pipeline.cpp:405-424` `prior_allele_type` still classifies `"*"` as symbolic → `OTHER`
  (`max(snp,indel)`, no `log10(3)`). Fixture-round site #5/#9: not reached by any input tried. The
  three new targets do not cover this path, and this round did not attempt it (out of the brief's
  three targets).
* `calling_pipeline.cpp:6059-6074` `gvcf_candidate_is_non_monomorphic` `<NON_REF>` `max()`
  (site #6/#10) — not reached.
* `hc_call.cpp:2938-2940` `estimate_mle_allele_counts(..., include_non_ref=true)` `max()`
  (site #7/#11) — fallback-only, not reached; see §4.2 for why relaxing the `heterozygosity` gate
  would make it reachable.
* `genotype_gvcf_tool.cpp:1585-1615` `'*'`→Other in GenotypeGVCFs — out of HC scope, untouched.
* The two separately-tracked FORMAT/serialization gaps (native emits no `GP`/`PG` under
  `--genotype-assignment-method USE_POSTERIOR_PROBABILITIES`; `AF`/`MLEAF` zero formatting on
  multi-ALT records). Both are still reported as notes by the new oracles and are unchanged.
* The `chr1:697` flanking-record divergence in the `del2` probes (native does not retain that site:
  zero-QUAL vs missing-QUAL plus missing rank-sum/`RAW_MQandDP`/`SB`). Still report-only.
* The zero-valued `--heterozygosity` / `--indel-heterozygosity` / `--heterozygosity-stdev` exit-code
  divergence (§4).

No new gate was registered: editing `fastgatk-native/CMakeLists.txt` is forbidden in this round, so
`ctest` still prints `Total Tests: 286` and the three site oracles must be run directly (they are
**not** protected by `ctest` alone).

---

## 7. What remains unproven

* **Precision, not bit-equality.** The fixes are validated only to the printed precision of the
  compared fields (QUAL/QD at 2 decimals, GQ as an integer). The internal `double` posteriors are
  not proven bit-identical to GATK's, and target 3's improvement is measured as ±4 integer GQ, not
  as a bitwise posterior comparison.
* **Target 1's line execution is still inferred, not instrumented** (as in the fixture round): the
  binary has no DWARF line information and the ordinary-path helper's `has_spanning_deletion`
  argument was constant-propagated. The evidence remains behavioural plus the cross-writer QUAL
  identity with target 2 (both writers now produce 368.90 at the same site from the same
  `REF/concrete/*` matrix).
* **Ploidy coverage.** Target 3 is proven for ploidy 3 (del 2/3/6) and ploidy 4 (del 2) under
  posterior assignment, and the fix also affects the diploid gVCF/ordinary prior calls
  (`hc_call.cpp:4881`) which the fixtures do **not** exercise with posterior assignment — the
  diploid behaviour change is source-derived only.
* **The `REF length > 1` branch of all three ternaries is not covered by any fixture.** Every
  fixture uses a 1 bp REF, so only the SNP side of each new ternary is observable here; the
  (unchanged) indel side rests on the source rule and on the previously landed siblings, not on a
  new measurement.
* **No unit-level test of `genotype_priors_for_group`** was added, so a future edit could silently
  revert target 3 for the diploid call sites without any registered gate failing.
* **The optional zero-valued-option work** is an investigation only (§4). Whether GATK's degenerate
  `--heterozygosity 0` / `--heterozygosity-stdev 0` output can be reproduced by native's clamped
  priors without reworking the AF-block guards is **speculation** and untested.
* **Nothing here proves the remaining `NOT contained` sites are harmless** — only that these
  fixtures do not reach them.
