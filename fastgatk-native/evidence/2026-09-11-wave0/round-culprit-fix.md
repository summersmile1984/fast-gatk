# Round: VariantRecalibrator `culprit` — GATK contract vs native provenance string

Label: `culprit-fix`
Scope: `fastgatk-native/src/variant_recalibrator_tool.cpp`,
`fastgatk-native/scripts/verify_variant_recalibrator.py`,
`fastgatk-native/scripts/verify_variant_recalibrator_gatk_model_oracle.py`,
new `fastgatk-native/scripts/verify_variant_recalibrator_culprit_gatk_oracle.py`.
No Mutect2 or unrelated tool was touched; no `CMakeLists.txt` or root `.md` was
touched; no commit/branch was created.

---

## 1. GATK truth (measured, not taken from the audit)

### 1.1 Jar grep (both strings, verbatim command)

```
$ unzip -p third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
      | strings | grep -c full-covariance-gmm
0
$ unzip -p third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
      | strings | grep -c culprit
13
```

The audit's grep claim is confirmed: `full-covariance-gmm` does not exist in the
pinned jar; `culprit` does (13 hits, including `AS_culprit` and `culpritString`).
Note `serialized-gmm` — the other native-only string — is likewise absent (it is
not a GATK string either; see §4).

### 1.2 Pinned GATK run on the exact `verify_variant_recalibrator.py` fixture

Fixture: 6 SNPs at `chr1:1..6` with `QD/MQ = (30,60),(28,58),(25,55),(4,25),(3,20),(2,15)`,
training = rows 1-3, known = rows 1-2. Command (verbatim, `WORK` = a
`tempfile.TemporaryDirectory`):

```
third_party/jdk17/bin/java -Xmx1g -jar third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
  VariantRecalibrator -V $WORK/input.vcf \
  --resource:truth,training=true,truth=true,known=false,prior=15.0 $WORK/training.vcf \
  --resource:known,training=false,truth=false,known=true $WORK/known.vcf \
  -an QD -an MQ --mode SNP --max-gaussians 1 --max-attempts 20 --k-means-iterations 20 \
  --bad-lod-score-cutoff 100.0 --dont-run-rscript true --create-output-variant-index false \
  --sites-only-vcf-output true -O $WORK/gatk.recal.vcf --tranches-file $WORK/gatk.tranches \
  --truth-sensitivity-tranche 100 --output-model $WORK/gatk.model     # exit 0
```

Exact data rows GATK wrote (all six; fields are tab-separated in the file and
shown whitespace-separated here, every token is verbatim):

```
chr1 1 . N <VQSR> . . END=1;NEGATIVE_TRAIN_SITE;POSITIVE_TRAIN_SITE;VQSLOD=1.4860;culprit=MQ
chr1 2 . N <VQSR> . . END=2;NEGATIVE_TRAIN_SITE;POSITIVE_TRAIN_SITE;VQSLOD=1.4860;culprit=MQ
chr1 3 . N <VQSR> . . END=3;NEGATIVE_TRAIN_SITE;POSITIVE_TRAIN_SITE;VQSLOD=1.4860;culprit=MQ
chr1 4 . N <VQSR> . . END=4;VQSLOD=-0.2329;culprit=MQ
chr1 5 . N <VQSR> . . END=5;VQSLOD=-0.2329;culprit=MQ
chr1 6 . N <VQSR> . . END=6;VQSLOD=-0.2329;culprit=MQ
```

**GATK emits `culprit=MQ` on this fixture. The audit's reading is confirmed**
(and it is *not* "something other than MQ", so the stop-condition in the brief
does not apply).

Native, before the fix (same inputs, same arguments minus GATK-only flags):

```
chr1 1 . N <VQSR> . . VQSLOD=1.4860;culprit=full-covariance-gmm;POSITIVE_TRAIN_SITE;NEGATIVE_TRAIN_SITE;END=1
...
chr1 6 . N <VQSR> . . VQSLOD=-0.2329;culprit=full-covariance-gmm;END=6
```

Note the `VQSLOD` column already matched GATK bit-for-bit (`1.4860` / `-0.2329`);
only the `culprit` token diverged.

### 1.3 Which log10 column the culprit corresponds to (GATK source, quoted)

```java
// gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/vqsr/VariantRecalibratorEngine.java:80-93
public void calculateWorstPerformingAnnotation( final List<VariantDatum> data, final GaussianMixtureModel goodModel, final GaussianMixtureModel badModel ) {
    for( final VariantDatum datum : data ) {
        int worstAnnotation = -1;
        double minProb = Double.MAX_VALUE;
        double worstValue = -1;
        for( int iii = 0; iii < datum.annotations.length; iii++ ) {
            final Double goodProbLog10 = goodModel.evaluateDatumInOneDimension(datum, iii);
            final Double badProbLog10 = badModel.evaluateDatumInOneDimension(datum, iii);
            if( goodProbLog10 != null && badProbLog10 != null ) {
                final double prob = goodProbLog10 - badProbLog10;
                if(prob < minProb) { minProb = prob; worstAnnotation = iii; worstValue = datum.annotations[iii];}
            }
        }
        datum.worstAnnotation = worstAnnotation;
        datum.worstValue = worstValue;
    }
}
```

```java
// .../vqsr/GaussianMixtureModel.java:208-222  -- the "log10 likelihood column"
public Double evaluateDatumInOneDimension( final VariantDatum datum, final int iii ) {
    if(datum.isNull[iii]) { return null; }
    final double[] pVarInGaussianLog10 = new double[gaussians.size()];
    int gaussianIndex = 0;
    for( final MultivariateGaussian gaussian : gaussians ) {
        pVarInGaussianLog10[gaussianIndex] = gaussian.pMixtureLog10;
        if (gaussian.pMixtureLog10 != Double.NEGATIVE_INFINITY) {
            pVarInGaussianLog10[gaussianIndex] += MathUtils.normalDistributionLog10(gaussian.mu[iii], gaussian.sigma.get(iii, iii), datum.annotations[iii]);
        }
        gaussianIndex++;
    }
    return nanTolerantLog10SumLog10(pVarInGaussianLog10); // Sum(pi_k * p(v|n,k))
}
```

```java
// .../vqsr/VariantDataManager.java:485  -- the emitted value
builder.attribute(GATKVCFConstants.CULPRIT_KEY, (datum.worstAnnotation != -1 ? annotationKeys.get(datum.worstAnnotation) : "NULL"));
```

So the culprit is **the name of the annotation dimension** whose
*good-minus-bad per-dimension log10 mixture likelihood* is smallest. In
`MathUtils.normalDistributionLog10(mean, sd, x)` the VQSR call site passes the
covariance diagonal `sigma.get(iii,iii)` — i.e. a **variance** in the `sd` slot
(`VariantDataManager.java:29` has the matching "this is really the standard
deviation" comment about GATK's VQSR naming). The comparison is strict
(`prob < minProb` with `minProb = Double.MAX_VALUE`), so ties keep the first
dimension of GATK's information-ordered `annotationKeys`.

### 1.4 Independent confirmation that the rule is per-datum, not a constant

On the 6-record fixture GATK's positive and negative models coincide (both are
fit on the three rows that survive `selectWorstVariants`, since the other three
fail `STD_THRESHOLD = 10`), so the contrast is exactly `0 - 0` for both
dimensions and the tie-break yields the first information-ordered annotation,
`MQ`. That makes the 6-record fixture a *weak* witness for the rule, so a second,
non-degenerate 60-record fixture was measured (40 training records with a wide
QD/MQ spread, 10 records extreme in QD, 10 extreme in MQ; GATK default gaussian
counts, `--bad-lod-score-cutoff 100.0`):

* GATK wrote `culprit=QD` on **26** records and `culprit=MQ` on **34** records —
  the value demonstrably varies per datum, so it can never be a model name.
* A Python re-implementation of `evaluateDatumInOneDimension` +
  `calculateWorstPerformingAnnotation`, fed **GATK's own `--output-model`
  GATKReport tables**, reproduced GATK's written culprit column on **60/60**
  records (0 mismatches). The same script fed **native's** model report
  reproduced native's column on **60/60** records.

This fixture is case `rich-per-datum` of the new oracle.

---

## 2. The new strict gate — `verify_variant_recalibrator_culprit_gatk_oracle.py`

`fastgatk-native/scripts/verify_variant_recalibrator_culprit_gatk_oracle.py`
(house conventions: `main() -> int`, `argparse` with `--native` /
`$FASTGATK_VARIANT_RECALIBRATOR_BINARY`, `--expect-divergence`,
`tempfile.TemporaryDirectory` scratch, `oracle_guard`, final JSON status line;
modelled on `verify_hc_af_zero_format_gatk_oracle.py`).

Three checks per case, each reported separately:

1. `cross_tool_byte_identity` — the `culprit=...` token of every data row,
   keyed by POS, must be byte-identical between pinned GATK and native.
2. `self_rule_fidelity` — each side's written culprit column must equal the
   `calculateWorstPerformingAnnotation` argmin re-derived *inside the oracle*
   from that side's own `--output-model` artifact.
3. `gatk_culprit_multiset` — the GATK side must still produce the measured
   annotation-name multiset (fixture-validity guard: no vacuous pass on a
   constant or on `NULL`).

Cases: `tiny-truth-known`, `tiny-truth-known-annotation-order-swapped` (both
gate 1 byte-identity, measured `{'MQ': 6}`), `rich-per-datum` (measured
`{'QD': 26, 'MQ': 34}`).

### Exit status before / after the fix

| run | command | exit | literal result |
|---|---|---|---|
| before (rev. A of the script, all three cases byte-identity-gated) | `python3 fastgatk-native/scripts/verify_variant_recalibrator_culprit_gatk_oracle.py` | **1** (3 violations) | GATK `MQ`/`QD` vs native `full-covariance-gmm` on 6/6, 6/6 and 60/60 records |
| after (rev. B, same command) | as above | **0** | 0 violations |

Rev. A→B difference: on `rich-per-datum`, gate 1 is reported rather than gated —
see §3.3 for the measured, quantified reason; gates 2 and 3 run on every case in
both revisions.

`--expect-divergence` was also run before the fix (exit 0, 3 violations
reported) to record the literal GATK-vs-native columns quoted in §1.2 and §3.3.

---

## 3. The fix (minimal, in `variant_recalibrator_tool.cpp`)

### 3.1 Diff applied

```diff
@@ struct Entry                              (variant_recalibrator_tool.cpp:211)
+    // Java VariantDatum.worstAnnotation: ...
+    int worst_annotation = -1;
@@ after mixture_log10_probability          (new; ~line 1869)
+double gatk_normal_distribution_log10(double mean, double sd, double x)      // MathUtils.normalDistributionLog10
+double mixture_log10_probability_one_dimension(const GaussianMixture&, double value, std::size_t dimension)
+void calculate_worst_annotation(std::vector<Entry>&, const GaussianMixture& good, const GaussianMixture& bad)
@@ write_recal_vcf signature                 (:2657)
+                     const std::unordered_map<std::string, std::string>& culprits,
@@ write_recal_vcf body                      (:2691, :2736, :2757)
+    const auto culprit_value = [&culprits](const std::string& key) -> std::string {...};  // "NULL" fallback
-                    culprits << model;              // AS_culprit
+                    culprit_list << culprit_value(keys[index]);
-    bcf_update_info_string(..., "CULPRIT", model.c_str());
+    bcf_update_info_string(..., "CULPRIT", culprit_value(keys.front()).c_str());
-    bcf_update_info_string(..., "culprit", model.c_str())
+    bcf_update_info_string(..., "culprit", culprit_value(keys.front()).c_str())
@@ main                                      (:2982)
+        calculate_worst_annotation(entries, good, bad);
+        std::unordered_map<std::string, std::string> culprits;   // key -> annotation name (or "NULL")
```

Nothing was restructured: scoring, model fitting, tranches and the recal VCF's
other columns are untouched. The `model` string still exists and is still used
for native-only provenance (manifest `"model"`, `##fastgatk_variant_recalibrator_status`,
VQSLOD header description) — those are not GATK fields and were deliberately left
alone (see §4, items 3-5).

### 3.2 Why it is a naming fix and not a numerical one

`calculate_worst_annotation` only *reads* the already-fitted models and the
already-normalized datum values; it feeds nothing back into VQSLOD, tranches or
the model report. Fidelity was checked by the oracle's gate 2 on all three
cases: **0 disagreements** between each side's written column and the rule
re-derived from that side's own model artifact (GATK 6/6, 6/6, 60/60; native
6/6, 6/6, 60/60). The pre-existing `VQSLOD` values are unchanged
(`1.4860` / `-0.2329` before and after).

### 3.3 Known residual divergence, quantified (reported, not hidden)

On `rich-per-datum`, 5 of 60 records get a different annotation from the two
tools: `POS 70 GATK=QD native=MQ`, `80 QD/MQ`, `120 MQ/QD`, `250 MQ/QD`,
`350 QD/MQ`. This is **not** a culprit-naming defect: both sides satisfy gate 2,
and the cause is the already-documented model-fidelity boundary (native's
multi-component VBEM fit is not bit-identical to Java's; cf.
`verify_variant_recalibrator_gatk_model_oracle.py`, which states that model
floating-point values remain a separate semantic boundary). Measured
`max |GATK − native|` over the shared model tables of that run:

| table | max abs delta |
|---|---|
| AnnotationMeans | 3.47e-07 |
| AnnotationStdevs | 1.51e-07 |
| PositiveModelMeans | **2.73** |
| NegativeModelMeans | **3.04** |
| PositiveModelCovariances | 0.25 |
| NegativeModelCovariances | **2.93** |

(On the 6-record fixture the same deltas are exactly `0.0` for all six tables,
which is why byte-identity is gated there.) The oracle prints these deltas and
the differing records in `cross_tool_differences` /
`model_max_abs_delta` on every run. *Speculation (labelled as such):* narrowing
this residual would require making the multi-component VBEM fit bit-identical to
Java's, which is outside a minimal culprit-naming fix.

---

## 4. Test lines corrected (all of them)

Every stale assertion that depended on the old string, and nothing else:

| file:line (after edit) | old value | new value | note |
|---|---|---|---|
| `verify_variant_recalibrator.py:84-98` (assert at `:98`) | `culprit=full-covariance-gmm` | `culprit=MQ` | measured GATK value on this fixture; comment cites `VariantRecalibratorEngine.java:80-93`, `GaussianMixtureModel.java:208-222`, `VariantDataManager.java:485` and the new oracle |
| `verify_variant_recalibrator.py:436-440` (`:439-440`) | `culprit=full-covariance-gmm` in `gmm_text` | `culprit=MQ` in `gmm_text` **and** the old string asserted absent | `--max-gaussians 2` case |
| `verify_variant_recalibrator.py:455-461` (`:460-461`) | `culprit=full-covariance-gmm` in `full_text` | `culprit=MQ` **and** the old string asserted absent | `--full-covariance` case |
| `verify_variant_recalibrator.py:480-498` (`:493-498`) | `AS_culprit=full-covariance-gmm` | `AS_culprit=full-covariance-gmm` asserted absent + exact per-record list `["QD","MQ","MQ","QD","QD","QD"]` + every entry ∈ `{"QD","MQ"}` (requested annotation names) | see the honesty note below |
| `verify_variant_recalibrator_gatk_model_oracle.py:180-189` (`:189`) | `native_info["culprit"] == "serialized-gmm"` | `native_info["culprit"] == gatk_info["culprit"] == "MQ"` | same defect class as `verify_indel.py`: a registered test pinning a native-only string where GATK has the annotation name. Not in the brief's list, but it *depends on the old behaviour*; the line above it (`:178`) already pins GATK's `MQ`, so the native side is now pinned to the same exact value |

Lines **not** touched, deliberately:

* `verify_variant_recalibrator.py:260/433/453` — `manifest_payload["model"] ==
  "full-covariance-gmm"`. The JSON manifest is native-only provenance metadata
  with no GATK counterpart (GATK writes no manifest), so this is not the culprit
  column and asserting it does not pin GATK-divergent behaviour.
* The C++ `model` local itself, used for the native manifest / VQSLOD header
  descriptions.

**Honesty note on the AS line.** The AS invocation in the registered test uses
training-only resources and the default `--bad-lod-score-cutoff (-5.0)`; pinned
GATK cannot build a negative model for exactly those arguments — measured: exit
2, `A USER ERROR has occurred: No data found.` (all data either pass the cutoff
or are excluded by `STD_THRESHOLD`). So there is no GATK byte-value to pin for
that configuration, and the assertion is: (a) the field no longer contains a
model name, (b) the exact native per-record values, (c) every entry is one of the
requested annotation names. Asserting `AS_culprit=MQ` there would have been
guesswork. Non-AS byte-identical GATK gating lives in the new oracle; the AS
column of a *runnable* GATK configuration was measured separately (with
`--bad-lod-score-cutoff 100.0`: GATK `QD,QD,MQ,MQ,MQ,MQ` vs native
`QD,MQ,MQ,QD,QD,QD`) and is **not** asserted anywhere as GATK parity — that is a
model-fidelity difference of the AS mode, reported here rather than hidden.

---

## 5. Gate results (mandatory order)

| step | command | exit |
|---|---|---|
| a | `python3 fastgatk-native/scripts/verify_variant_recalibrator_culprit_gatk_oracle.py` | **0** |
| a' | same, `--expect-divergence` (diagnostic mode) | **0** |
| b | `python3 fastgatk-native/scripts/verify_variant_recalibrator.py` | **0** |
| b' | `python3 fastgatk-native/scripts/verify_variant_recalibrator_gatk_model_oracle.py` (corrected) | **0** |
| b'' | `ctest -R 'variant-recalibrator' -V` (7 registered recalibrator tests, omp build) | **0** (7/7) |
| c | `ctest -R 'window-invariance-gatk-oracle\|alleles-overlap-gate-oracle\|span-del-qual-gatk-oracle\|gvcf-symbolic-prior-gatk-oracle\|af-zero-format-gatk-oracle\|arbitrary-ploidy-span-del-prior-oracle\|polyploid-gvcf-span-del-prior-oracle\|spanning-prior-genotype-gq-oracle\|multialt-owner-annotation-oracle\|mutect2-recheck\|gvcf-indel-end-gatk-oracle' -V` (in `fastgatk-native/build`) | **0** (13/13 passed, 834.6 s) |
| d | `fastgatk-native/scripts/run_regression.sh --label culprit-fix` | **0** (294/294 both backends) |
| a/b/b' re-confirmation | all three re-run last, on the final tree | **0 / 0 / 0** |

Rebuilds performed before (d): `cmake --build fastgatk-native/build --target
fastgatk-variant-recalibrator` and `.../build-serial --target
fastgatk-variant-recalibrator` (both exit 0). `fastgatk-hc-call` was **not**
modified: it does not link `variant_recalibrator_tool.cpp`, and a no-op
`cmake --build ... --target fastgatk-hc-call` in both trees reports
`[100%] Built target fastgatk-hc-call` (nothing to rebuild), so the HC binaries
already match their own sources. The runner's staleness *warning*
(`源码比二进制新 ... variant_recalibrator_tool.cpp`) is its conservative
any-source-newer-than-`fastgatk-hc-call` heuristic firing on a file that belongs
to a different tool; the binary that file belongs to was rebuilt in both trees.

### (d) `run_regression.sh --label culprit-fix`

```
# 双后端回归证据  (2026-09-11 15:29:32 CST, 标签 culprit-fix, filter <全量>, jobs 8)
| 后端   | 构建目录                     | 结果 | 通过/总数 | 耗时         |
| omp    | OpenMP (fastgatk-native/build)        | 通过 | 294/294 | 1217.43 sec |
| serial | Serial (fastgatk-native/build-serial) | 通过 | 294/294 | 1195.50 sec |
REGRESSION_EXIT=0
100% tests passed, 0 tests failed out of 294      (omp.log)
100% tests passed, 0 tests failed out of 294      (serial.log)
evidence: .diag/regression/20260911-150914/
```

`FASTGATK_REQUIRE_GATK_ORACLE=1` is the runner default, so every GATK oracle in
the 294 really ran (no silent skips). Note the run happened *before* the extra
serial-backend oracle invocation below only in wall-clock terms; no source or
binary changed between them.

### (e) Extra: the new oracle on the serial build

```
$ python3 fastgatk-native/scripts/verify_variant_recalibrator_culprit_gatk_oracle.py \
      --native fastgatk-native/build-serial/fastgatk-variant-recalibrator     # exit 0
NATIVE culprit multiset (tiny x2): {'MQ': 6} / {'MQ': 6}
NATIVE culprit multiset (rich):    {'QD': 25, 'MQ': 35}   (GATK {'QD': 26, 'MQ': 34})
rule fidelity: gatk 0 disagreements, native 0 disagreements on all three cases
cross_tool_differences (rich): POS 70/80/350 GATK=QD native=MQ, POS 120/250 GATK=MQ native=QD
```

Identical to the OpenMP build, including the same five differing records — so the
fix is backend-deterministic.

The new oracle is **not** registered in CTest: the brief forbids editing
`fastgatk-native/CMakeLists.txt` (gate registration is the orchestrator's job).
It is runnable standalone and would slot in next to
`fastgatk-variant-recalibrator-contract` as e.g.
`fastgatk-variant-recalibrator-culprit-gatk-oracle`.

---

## 6. State of the tree

**The tree CONTAINS the change** (no revert was needed — every criterion above
passed):

```
$ git status --porcelain
 M fastgatk-native/scripts/verify_variant_recalibrator.py
 M fastgatk-native/scripts/verify_variant_recalibrator_gatk_model_oracle.py
 M fastgatk-native/src/variant_recalibrator_tool.cpp
?? fastgatk-native/scripts/verify_variant_recalibrator_culprit_gatk_oracle.py
```

`git diff --stat`: 3 files changed, 183 insertions(+), 12 deletions(-), plus the
new oracle script. No commit, no branch, no `git checkout` of whole files; all
edits are targeted. Scratch (GATK fixtures, probe scripts) lived only in
`tempfile.TemporaryDirectory` / throwaway `/tmp` dirs.

The reported 294/294 regression therefore describes the tree *with* the change
(`git c492165` + 4 uncommitted paths, per the runner's own evidence block).

## 7. What remains unproven

1. **AS-mode GATK parity**: native's per-allele `AS_culprit`/`AS_VQSLOD` do not
   match pinned GATK for a configuration GATK can actually run (measured above);
   the AS column is therefore pinned to the GATK contract + native regression
   values, not to GATK bytes. Not investigated further (AS-mode model fidelity,
   not the culprit string).
2. **`rich-per-datum` byte-identity**: 5/60 records differ, for the measured
   model-numeric reason in §3.3. Gate 2 passes on both sides, so the *rule* is
   proven; only cross-tool equality of the underlying multi-gaussian fit is not.
3. **The `--input-model` path was verified only on the oracle's 6-record fixture**
   (`verify_variant_recalibrator_gatk_model_oracle.py`, `culprit == "MQ"` on both
   sides). A richer replayed-model fixture was not run.
4. **`NULL` fallback is implemented but not exercised**: no fixture with all-null
   annotation dimensions was constructed, so `culprit=NULL` is code-path
   reasoning from `VariantDataManager.java:485`, not a measurement.
5. *Speculation:* whether a real production callset would ever hit the tie-break
   (identical positive/negative models) as the 6-record fixture does — no
   evidence either way.
