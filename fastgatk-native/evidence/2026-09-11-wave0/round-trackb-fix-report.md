# Round: Track B — `--alleles` overlapping forced-feature emission fix

Date: 2026-09-11. Repo: `/home/turing-agents/Documents/fast-gatk`, git `bd0e8f4`
(two uncommitted working-tree changes, listed in §1). **No `git commit`, no branch,
no edit to `fastgatk-native/CMakeLists.txt`, no edit to any root `*.md`.**
Binaries: `fastgatk-native/build/fastgatk-hc-call` (OpenMP), `fastgatk-native/build-serial/fastgatk-hc-call`
(Serial). Oracle: pinned GATK 4.6.2.0
(`third_party/jdk17/bin/java -Xmx1g -jar third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar HaplotypeCaller`).

**Bottom line: the defect is localized to a single statement, fixed with a minimal
change, and every gate passes (new strict oracle exit 0; both window-invariance
gates pass; the existing `--alleles` oracle passes; dual-backend regression
283/283 on OpenMP *and* Serial). The tree CONTAINS the fix.**

---

## 1. Changes in the working tree

```
$ git status --short
 M fastgatk-native/src/calling_pipeline.cpp
?? fastgatk-native/scripts/verify_hc_alleles_overlap_gate_oracle.py
```

`git diff` (production change, verbatim):

```diff
diff --git a/fastgatk-native/src/calling_pipeline.cpp b/fastgatk-native/src/calling_pipeline.cpp
index 482550f..e276f09 100644
--- a/fastgatk-native/src/calling_pipeline.cpp
+++ b/fastgatk-native/src/calling_pipeline.cpp
@@ -13757,6 +13757,23 @@ void accumulate_independent_region_telemetry(Result& merged, const Result& part)
         merged.pairhmm_skip_reason = part.pairhmm_skip_reason;
 }
 
+// GATK MathUtils.log10OneMinusPow10 (via NaturalLogUtils.log1mexp): the
+// log10 of `1 - 10^log10_value` computed without cancelling away the tiny
+// complement.  GenotypingEngine uses this to turn an
+// AlleleFrequencyCalculator `log10ProbOnlyRefAlleleExists` posterior into the
+// complementary `log10ProbVariantPresent` site confidence that a monomorphic
+// locus reports as its QUAL.
+double log10_one_minus_pow10(const double log10_value) {
+    if (log10_value > 0.0) return std::numeric_limits<double>::quiet_NaN();
+    if (log10_value == 0.0) return -std::numeric_limits<double>::infinity();
+    constexpr double kLog10 = 2.30258509299404568402;
+    constexpr double kLog2 = 0.69314718055994530942;
+    const double natural = log10_value * kLog10;
+    const double log1mexp = natural > -kLog2
+        ? std::log(-std::expm1(natural)) : std::log1p(-std::exp(natural));
+    return log1mexp / kLog10;
+}
+
 }  // namespace
 
 std::optional<GenotypeCall::Annotations> calculate_output_variant_annotations(
@@ -17305,8 +17322,40 @@ Result run(const io::ReadBatch& reads,
                     {}, {}, 2);
                 if (candidate < confidence_qual.size() &&
                     spanning_af.samples_with_likelihoods != 0 &&
-                    std::isfinite(spanning_af.qual))
-                    confidence_qual[candidate] = spanning_af.qual;
+                    std::isfinite(spanning_af.qual)) {
+                    // GenotypingEngine picks the site-confidence branch from this
+                    // same locus AF result: `log10ProbOnlyRefAlleleExists` (the
+                    // kernel's `qual`) when any alternate allele is plausible,
+                    // and the complementary `log10ProbVariantPresent` when every
+                    // alternate allele fails that test
+                    // (GenotypingEngine.calculateGenotypes,
+                    //  `!outputAlternativeAlleles.siteIsMonomorphic`).
+                    // A forced GenotypeGivenAlleles allele is retained in the
+                    // output even when it is not plausible, so a forced event
+                    // that an already-assembled deletion spans stays
+                    // monomorphic and must keep the complementary QUAL that the
+                    // forced-allele policy above just computed.  Overwriting it
+                    // with the polymorphic-branch value (~0 on a reference-only
+                    // pileup) dropped the forced row below
+                    // --standard-confidence-for-calling, which is the
+                    // overlapping-feature emission defect.
+                    bool loci_monomorphic = true;
+                    const double plausibility_log10 = -0.1 * options.standard_confidence_for_calling;
+                    for (int allele = 1;
+                         allele < static_cast<int>(spanning_af.log10_p_allele_absent.size());
+                         ++allele)
+                        loci_monomorphic = loci_monomorphic &&
+                            !(spanning_af.log10_p_allele_absent[static_cast<std::size_t>(allele)] +
+                              1.0e-10 < plausibility_log10);
+                    const double complement = loci_monomorphic
+                        ? log10_one_minus_pow10(spanning_af.log10_p_no_variant)
+                        : std::numeric_limits<double>::quiet_NaN();
+                    if (result.candidates[candidate].forced_by_alleles_feature &&
+                        std::isfinite(complement))
+                        confidence_qual[candidate] = -10.0 * complement;
+                    else
+                        confidence_qual[candidate] = spanning_af.qual;
+                }
             }
         }
     }
```

`fastgatk-native/scripts/verify_hc_alleles_overlap_gate_oracle.py` is new (strict
gate, §2a). Nothing else was added, moved or deleted: the temporary
`std::cerr` instrumentation used in §3 is **not** in the tree
(`grep -n "TEMP_PROBE\|TEMP_INSTRUMENT" fastgatk-native/src/calling_pipeline.cpp` → no match;
the scratch fixture/helper under `.diag/` were deleted).

Diff size: ` 1 file changed, 51 insertions(+), 2 deletions(-)`.

---

## 2. Gates (literal output, in the mandated order)

### 2a) New strict oracle — `fastgatk-native/scripts/verify_hc_alleles_overlap_gate_oracle.py`

Self-contained: it builds the fixture in a `tempfile.TemporaryDirectory`
(1500 bp `chr1`, 30×300M reference-only reads over the `8×CAG` tandem repeat at
900, plus 6 reads carrying an assembled 1 bp deletion at 700), runs the pinned
GATK jar and the native binary with identical arguments for each case **and**
for its `--drop-alleles` attribution control, and asserts the native data rows
are byte-identical to GATK's (header/provenance ignored). `--expect-divergence`
inverts it into a diagnostic that exits 0. Final JSON status line; `main() -> int`.

**Pre-fix run (before any production change was made) — exit 1, status `fail`:**

```
[overlapping-features / alleles] gatk_rows=2 native_rows=1 data_rows_byte_identical=False divergent_positions=['missing-in-native@1:POS=904']
    feature VCF rows: chr1	900	.	CAGCAG	C	.	PASS	.; chr1	904	.	A	T	.	PASS	.
    GATK   chr1	900	.	CAGCAG	C	126.81	.	AC=0;AF=0.00;AN=2;DP=30;ExcessHet=0.0000;FS=0.000;MLEAC=0;MLEAF=0.00;MQ=60.00;SOR=0.001	GT:AD:DP:GQ:PL	0/0:30,0:30:90:0,90,1350
    GATK   chr1	904	.	A	T	117.78	.	AC=0;AF=0.00;AN=2;DP=30;ExcessHet=0.0000;FS=0.000;MLEAC=0;MLEAF=0.00;MQ=60.00;SOR=0.001	GT:AD:DP:GQ:PL	0/0:30,0:30:90:0,90,1343
    NATIVE chr1	900	.	CAGCAG	C	126.81	.	AC=0;AF=0.00;AN=2;DP=30;ExcessHet=0.0000;FS=0.000;MLEAC=0;MLEAF=0.00;MQ=60.00;SOR=0.001	GT:AD:DP:GQ:PL	0/0:30,0:30:90:0,90,1350
[overlapping-features / drop-alleles] gatk_rows=0 native_rows=0 data_rows_byte_identical=True divergent_positions=[]
[assembled-event-overlap / alleles] gatk_rows=2 native_rows=1 data_rows_byte_identical=False divergent_positions=['missing-in-native@1:POS=700']
    GATK   chr1	697	.	TC	T	98.60	.	AC=1;...;QD=2.74;ReadPosRankSum=3.417;SOR=0.050	GT:AD:DP:GQ:PL	0/1:30,6:36:99:106,0,961
    GATK   chr1	700	.	CA	C	144.77	.	AC=0;AF=0.00;AN=2;DP=36;...;MLEAC=0;MLEAF=0.00;MQ=60.00;SOR=0.001	GT:AD:DP:GQ:PL	0/0:36,0:36:99:0,108,1481
    NATIVE chr1	697	.	TC	T	98.60	.	AC=1;...;QD=2.74;ReadPosRankSum=3.417;SOR=0.050	GT:AD:DP:GQ:PL	0/1:30,6:36:99:106,0,961
[assembled-event-overlap / drop-alleles] gatk_rows=1 native_rows=1 data_rows_byte_identical=True divergent_positions=[]
# 2 violation(s):
#   - overlapping-features [alleles]: data rows differ at ['missing-in-native@1:POS=904'] (gatk_rows=2 native_rows=1)
#   - assembled-event-overlap [alleles]: data rows differ at ['missing-in-native@1:POS=700'] (gatk_rows=2 native_rows=1)
STATUS: fail            (exit 1)
```

The two `--drop-alleles` controls reproduce the audit's attribution result
(`0` rows vs `0` rows for D2; `1` identical row for D1), so both lost rows are
`--alleles`-driven and the fixture is a valid forced-allele-emission probe. This
run is also the gate's **sensitivity check**: on the unfixed tree it fails.

**Post-fix run — exit 0, status `pass`** (`timeout 900 python3 fastgatk-native/scripts/verify_hc_alleles_overlap_gate_oracle.py`):

```
[overlapping-features / alleles] gatk_rows=2 native_rows=2 data_rows_byte_identical=True divergent_positions=[]
    GATK   chr1	900	.	CAGCAG	C	126.81	.	AC=0;AF=0.00;AN=2;DP=30;ExcessHet=0.0000;FS=0.000;MLEAC=0;MLEAF=0.00;MQ=60.00;SOR=0.001	GT:AD:DP:GQ:PL	0/0:30,0:30:90:0,90,1350
    GATK   chr1	904	.	A	T	117.78	.	AC=0;AF=0.00;AN=2;DP=30;ExcessHet=0.0000;FS=0.000;MLEAC=0;MLEAF=0.00;MQ=60.00;SOR=0.001	GT:AD:DP:GQ:PL	0/0:30,0:30:90:0,90,1343
    NATIVE chr1	900	.	CAGCAG	C	126.81	.	AC=0;AF=0.00;AN=2;DP=30;ExcessHet=0.0000;FS=0.000;MLEAC=0;MLEAF=0.00;MQ=60.00;SOR=0.001	GT:AD:DP:GQ:PL	0/0:30,0:30:90:0,90,1350
    NATIVE chr1	904	.	A	T	117.78	.	AC=0;AF=0.00;AN=2;DP=30;ExcessHet=0.0000;FS=0.000;MLEAC=0;MLEAF=0.00;MQ=60.00;SOR=0.001	GT:AD:DP:GQ:PL	0/0:30,0:30:90:0,90,1343
[overlapping-features / drop-alleles] gatk_rows=0 native_rows=0 data_rows_byte_identical=True divergent_positions=[]
[assembled-event-overlap / alleles] gatk_rows=2 native_rows=2 data_rows_byte_identical=True divergent_positions=[]
    NATIVE chr1	697	.	TC	T	98.60	.	AC=1;...;QD=2.74;ReadPosRankSum=3.417;SOR=0.050	GT:AD:DP:GQ:PL	0/1:30,6:36:99:106,0,961
    NATIVE chr1	700	.	CA	C	144.77	.	AC=0;AF=0.00;AN=2;DP=36;...;MLEAC=0;MLEAF=0.00;MQ=60.00;SOR=0.001	GT:AD:DP:GQ:PL	0/0:36,0:36:99:0,108,1481
[assembled-event-overlap / drop-alleles] gatk_rows=1 native_rows=1 data_rows_byte_identical=True divergent_positions=[]
STATUS: pass            (exit 0)
```

Every field matches, including the recovered rows' QUAL (`117.78`, `144.77`,
GATK's exact values), FILTER (`.`), INFO and `GT:AD:DP:GQ:PL`.

Same oracle on the **Serial** backend (`--native fastgatk-native/build-serial/fastgatk-hc-call`): exit 0,
all four comparisons `data_rows_byte_identical=True`.
`--expect-divergence`: exit 0, `"status": "diagnostic"`.

### 2b) Pre-existing window-invariance gates

```
$ third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest --test-dir fastgatk-native/build -R 'window-invariance-gatk-oracle' -V
1/2 Test #83: fastgatk-hc-window-invariance-gatk-oracle ..........   Passed   60.22 sec
2/2 Test #84: fastgatk-hc-ploidy-window-invariance-gatk-oracle ...   Passed   34.61 sec
100% tests passed, 0 tests failed out of 2
```
Exit 0. Raw log: `.diag/ctest_window_after.log`.

### 2c) Existing `--alleles` oracle

```
$ timeout 1800 python3 fastgatk-native/scripts/verify_hc_alleles_gatk_oracle.py
{"empty_pileup_feature_no_call_gatk_exact": true, "filtered_feature_default_ignored": true,
 "force_call_filtered_alleles_and_legacy_alias_gatk_exact": true, "forced_feature_active": true,
 "gvcf_forced_alt_and_non_ref_gatk_exact": true, "java_native_records_exact": true,
 "multi_alt_and_anchored_indel_gatk_exact": true, "release": "GATK 4.6.2.0", "status": "pass"}
```
Exit 0. Raw log: `.diag/alleles_oracle_after.log`.

### 2d) Full dual-backend regression

```
$ fastgatk-native/scripts/run_regression.sh --label 'Track B overlapping-feature fix'
| 后端 | 构建目录 | 结果 | 通过/总数 | 耗时 |
| omp    | OpenMP (fastgatk-native/build)        | 通过 | 283/283 | 1047.50 sec |
| serial | Serial (fastgatk-native/build-serial) | 通过 | 283/283 | 1022.39 sec |
证据目录：.diag/regression/20260911-042712
```
Exit 0. Both `omp.log` and `serial.log` contain
`100% tests passed, 0 tests failed out of 283`. The driver's staleness pre-check
reported no "source newer than binary" warning, i.e. both HC binaries were built
from the fixed source. Raw driver output: `.diag/regression_after.log`.

### 2e) Extra (non-mandatory) regression probe

`verify_hc_alleles_deep_boundary.py` is not registered in CTest, so its
previously-passing cases were re-run explicitly with the fixed binary:

```
$ timeout 1800 python3 fastgatk-native/scripts/verify_hc_alleles_deep_boundary.py \
      --case homopolymer-insertion-left --case homopolymer-insertion-right \
      --case homopolymer-deletion --case tandem-repeat-deletion-anchored \
      --case adjacent-snps --case region-boundary-span
... "diverged_cases": [], "release": "GATK 4.6.2.0", "status": "pass"
```
Exit 0; all six cases `status: match` (raw log `.diag/deep_boundary_subset_after.log`).
The deliberately-failing cases of that script (`tandem-repeat-one-unit-literal` = the
out-of-scope D3 exit-code divergence, and `assembled-event-overlap` /
`overlapping-features` which D1/D2 fixed and which the new gate now covers strictly)
were not part of this subset.

---

## 3. Localization and mechanism

**Function: `Result run(...)` at `fastgatk-native/src/calling_pipeline.cpp:13975` (pre-fix line number for the
defect: 17309, i.e. `confidence_qual[candidate] = spanning_af.qual;` inside
`if (!options.somatic_mode && !result.spanning_deletion_read_likelihoods.empty())`).
Fix site now at lines 17342–17357; new helper at line 13766.**

Evidence chain (temporary `std::cerr` under `FASTGATK_TEMP_PROBE_ALLELES_QUAL`,
removed from the tree; `FASTGATK_DEBUG_CALL_GATES=1` output is pre-existing):

```
===== overlapping-features (chr1:800-1080) =====
[TEMP_PROBE_QUAL_PRE] i=0 pos=899 ref=CAGCAG alt=C forced=1 default_qual=9.04539e-13 log10_p_alt_absent=-9.04539e-14 log10_p_variant_present=-12.6814
[TEMP_PROBE_QUAL_PRE] i=1 pos=903 ref=A alt=T forced=1 default_qual=7.23824e-12 log10_p_alt_absent=-7.23824e-13 log10_p_variant_present=-11.7782
[TEMP_PROBE_QUAL_POST_FORCED] i=0 pos=899 forced=1 qual=126.814
[TEMP_PROBE_QUAL_POST_FORCED] i=1 pos=903 forced=1 qual=117.782
[TEMP_PROBE_QUAL_SPANNING] i=1 pos=903 forced=1 previous_qual=117.782 spanning_qual=7.23824e-12 spanning_log10_p_no_variant=-7.23824e-13 samples=1
[FASTGATK_CALL_GATE] pos=899 ref=CAGCAG alt=C support=0 somatic_pair_evidence=0 genotype=0 qual=126.814 decision=emit
[FASTGATK_CALL_GATE] pos=903 ref=A alt=T support=0 somatic_pair_evidence=0 genotype=0 qual=7.23824e-12 threshold=30 decision=confidence_suppressed

===== assembled-event-overlap (chr1:500-780) =====
[TEMP_PROBE_QUAL_PRE] i=2 pos=699 ref=CA alt=C forced=1 default_qual=1.44649e-14 log10_p_alt_absent=-1.44649e-15 log10_p_variant_present=-14.4775
[TEMP_PROBE_QUAL_POST_FORCED] i=2 pos=699 forced=1 qual=144.775
[TEMP_PROBE_QUAL_SPANNING] i=2 pos=699 forced=1 previous_qual=144.775 spanning_qual=1.44649e-14 spanning_log10_p_no_variant=-1.44649e-15 samples=1
[FASTGATK_CALL_GATE] pos=699 ref=CA alt=C support=0 somatic_pair_evidence=0 genotype=0 qual=1.44649e-14 threshold=30 decision=confidence_suppressed
```

**Mechanism (measured, not hypothesised).** The forced-allele policy loop
(`calling_pipeline.cpp:17286-17296`) correctly replaces the candidate's QUAL with
GATK's *complementary* `P(variant present)` site confidence — 117.782 for the
`904 A>T` SNP and 144.775 for the `700 CA>C` deletion, exactly GATK's emitted
117.78 / 144.77. Immediately afterwards, the spanning-deletion override block
unconditionally overwrote `confidence_qual[candidate]` with the
**polymorphic-branch** AF value `spanning_af.qual` (7.23824e-12 / 1.44649e-14 —
in this fixture numerically identical to the plain biallelic default). That value
is ~0 for a reference-only pileup, so the subsequent
`qual + 1e-10 < options.standard_confidence_for_calling` gate
(`calling_pipeline.cpp:17521`) discarded the forced row
(`decision=confidence_suppressed`) and the VCF lost it. The candidate itself was
never lost earlier: `[FASTGATK_EVENTMAP_INJECT] path=2 … appended=1`,
`[FASTGATK_EVENTMAP_REPLAY] … reason=already-present` and
`[FASTGATK_EVENTMAP_REGION] … calling_class=1` all fire for both events; the loss
is exactly at the QUAL clobber.

Why the override runs for these candidates at all: the second forced event lies
inside the REF span of the first (an injected/assembled deletion), so
`derive_pairhmm_spanning_deletion_pl()` succeeds for it and the block *re-runs*
the AF calculation with an explicit `*` allele.

**GATK's rule that the fix restores** (`gatk-source/.../genotyper/GenotypingEngine.java:158-160`
and `296-322`, `afcalc/AFCalculationResult.java:113-138`):

```java
final double log10Confidence =
    !outputAlternativeAlleles.siteIsMonomorphic || configuration.annotateAllSitesWithPLs
        ? AFresult.log10ProbOnlyRefAlleleExists() + 0.0 : AFresult.log10ProbVariantPresent() + 0.0;
```
where `siteIsMonomorphic &= !isPlausible(alt)` and
`isPlausible(alt) = getLog10PosteriorOfAlleleAbsent(alt) + EPSILON < qualToErrorProbLog10(standardConfidenceForCalling)`.
Because a given allele is retained in the output even when it is not plausible
(`toOutput = isPlausible || … || forcedAlleles.contains(allele)`), a forced event
that no other allele makes plausible keeps the **monomorphic** branch, i.e. one
AF result per locus and the correct branch of *that* result. The fix applies
exactly that branch selection to the spanning (`REF, concrete ALT, *`) AF result,
and — to stay minimal and avoid changing any non-`--alleles` behaviour — takes
the complementary value only for `forced_by_alleles_feature` candidates; every
other candidate still receives `spanning_af.qual` exactly as before.

**Post-fix call gate for the counterexample** (both rows now emitted):

```
[FASTGATK_CALL_GATE] pos=899 ref=CAGCAG alt=C support=0 somatic_pair_evidence=0 genotype=0 qual=126.814 decision=emit
[FASTGATK_CALL_GATE] pos=903 ref=A alt=T support=0 somatic_pair_evidence=0 genotype=0 qual=117.782 decision=emit
```

---

## 4. Does the tree contain the change?

**Yes — the tree CONTAINS the fix.** `fastgatk-native/src/calling_pipeline.cpp`
is modified (working tree, uncommitted) and both HC binaries were rebuilt from
it (OpenMP before the new-oracle and ctest gates, Serial before the regression);
`fastgatk-native/scripts/verify_hc_alleles_overlap_gate_oracle.py` is the new
untracked gate. No temporary instrumentation and no scratch fixture remain.
`fastgatk-native/CMakeLists.txt` and all root `*.md` files are untouched
(the new oracle is therefore *not* registered in CTest — the orchestrator
registers oracles, per the task constraints).

---

## 5. What remains unproven / not fixed

1. **The `loci_monomorphic` branch is verified only on this fixture family.**
   Observed: for a reference-only pileup the spanning AF result is monomorphic
   and the complement is taken; that reproduces GATK's `117.78` / `144.77`
   exactly. The *polymorphic* leg (a forced allele at a locus where the
   spanning-deletion `*` allele itself is plausible) was probed once and native
   keeps `spanning_af.qual` there, which matches GATK's branch choice in that
   probe — but in that probe **GATK also emits a `*,A` multi-allelic row with
   `QUAL=0`, `FILTER=LowQual` that native does not emit at all.** That is a
   *different, pre-existing* divergence family (see (2)), so the polymorphic leg
   is only weakly exercised, not proven.
2. **A newly observed, NOT-fixed divergence of the same `--alleles` area**
   (out of this round's minimal scope; reported, not fixed). Probe: the D1
   fixture plus a single forced SNP record at the base the assembled deletion
   removes (`chr1 698 . C A`, interval `chr1:500-780`):
   ```
   GATK   chr1 697 . TC T  98.60 . AC=1;…                GT:AD:DP:GQ:PL  0/1:30,6:36:99:106,0,961
   GATK   chr1 698 . C  *,A 0     LowQual AC=1,0;…       GT:AD:DP:GQ:PL  0/1:30,6,0:36:99:106,0,961,207,1055,1478
   NATIVE chr1 697 . TC T  98.60 . AC=1;…                GT:AD:DP:GQ:PL  0/1:30,6:36:99:106,0,961
   ```
   Native's candidate at 697 (0-based) is injected, classified and reaches the
   call gate with `qual=0 → decision=confidence_suppressed`. Two GATK rules are
   involved and neither is native's: (i) `GenotypingEngine.calculateGenotypes`
   bypasses the emit-threshold null return when `forcedAlleles` is non-empty
   (`&& forcedAlleles.isEmpty()`), so a low-confidence forced row is still
   emitted and merely gets `FILTER=LowQual`; (ii) the emitted record carries the
   `*` spanning-deletion ALT. **This divergence is pre-existing and is NOT
   caused by this fix**: by case analysis on the patched statement, when
   `spanning_af.qual == 0` (i.e. `log10_p_no_variant == 0`) the assignment is
   `spanning_af.qual` on both the old and the new code path (either
   `loci_monomorphic` is false, or `log10_one_minus_pow10(0) = -inf` is not
   finite and the `else` branch runs), so that candidate's confidence is
   bit-identical before and after. It is a *different* defect from D1/D2 (it is
   not a `--alleles` QUAL-branch bug, and the D1/D2 fix would be wrong for it:
   GATK's D1/D2 rows have `FILTER="."` and QUAL above threshold, so a generic
   "forced alleles always emit + LowQual" rule would have diverged from GATK on
   D1/D2). It is **not covered by the new gate** and is left open on purpose:
   fixing it needs a `*`-ALT emission path plus a forced-allele emission bypass,
   which is neither minimal nor safe inside this round. **Recommendation:** open
   a separate track for "forced allele emission bypasses `passesEmitThreshold`
   and is only `LowQual`-filtered" (+ `*` ALT output), with its own strict gate.
3. **The out-of-scope D3 exit-code divergence** (GATK exit 3 /
   `Null alleles are not supported` on a `--alleles` record that minimises to an
   empty ALT, native exit 0) was deliberately not touched; `tandem-repeat-one-unit-literal`
   still diverges and is not part of the new gate.
4. **Unchanged, still-untested from the prior audit:** base-haplotype
   score-ranking tie behaviour (`NUM_HAPLOTYPES_TO_INJECT_FORCE_CALLING_ALLELES_INTO`),
   `--max-genotype-count` × NaN injected-haplotype scores, multi-sample/joint
   genotyping, symbolic feature ALTs and Mutect2's separate
   `restrict_to_forced_alleles` replay. All 283 registered tests pass on both
   backends, and the six previously-passing deep-boundary cases still match, but
   those do not cover the above.
5. **The new oracle is not registered in CTest** (CMakeLists.txt is off-limits
   this round), so it is not part of the 283/283 count; it must be run directly
   or registered by the orchestrator.
