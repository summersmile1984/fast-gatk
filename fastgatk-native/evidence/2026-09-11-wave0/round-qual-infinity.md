# Round: GATK's `Infinity` QUAL on the dense-mode REF-only row (`covered-star-only-record-dense`)

Scope: the REPORTED-ONLY case `covered-star-only-record-dense` in
`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` — the
REF-only no-call row that pinned GATK `GenotypeGVCFs --include-non-variant-sites`
materializes for a locus whose only surviving ALT was the symbolic spanning
deletion `*`.  Pinned GATK writes that row with QUAL rendered as the literal
token `Infinity`; native wrote `0`.

## 0. Verdict

* **GATK truth (measured, pinned 4.6.2.0):** the row is
  `chr1 3 . A . Infinity . DP=20;MLEAC=.;MLEAF=. GT ./.` — the token is exactly
  `Infinity`, no quotes, no `inf`, no `.`.
* **The token is the rendering of a genuine `Double.POSITIVE_INFINITY`**, not a
  sentinel and not a saturating conversion.  For a monomorphic site
  `GenotypingEngine.java:158-163` selects the *complementary* confidence
  `AFresult.log10ProbVariantPresent()`, `MathUtils.log10OneMinusPow10(0.0)` is
  `Double.NEGATIVE_INFINITY`, and `:183` assigns that straight into the record as
  `log10PError`; phred is `-10 * log10PError` = `+Infinity`.  htsjdk then writes
  it through plain decimal formatting, so the token is just Java's `%.2f`
  spelling of an infinite double (§1.3).
* **The native defect was a *numeric* guard, with a rendering half.**  Native's
  `log10_one_minus_pow10` already returns `-Infinity` for this input (the model
  agrees), but `materialize_gatk_monomorphic_ref_call` wrapped the assignment in
  `if (std::isfinite(complement))`, so the infinite QUAL was **never assigned**
  and the record kept the value left by the earlier AF pass (`0` here,
  `159.55` in the previous round's `implausible-upstream-plausible-star-dense`
  probe).  Additionally `gatk_qual_output()` maps a non-finite double to `0.0F`,
  and `format_gatk_qual_value()` had no `Infinity` spelling for HTSlib's `inf`.
  So: one *numeric* suppression plus two *rendering* gaps.
* **Minimal fix, kept in the tree:** three lines of behaviour in one translation
  unit (`fastgatk-native/src/genotype_gvcf_tool.cpp`) — publish the infinite
  QUAL on the `-Infinity` complement branch only, and translate HTSlib's
  `inf`/`-inf` token to `Infinity`/`-Infinity` in the GATK-compat QUAL
  formatter.  No restructuring, no kernel, no ABI, no `CMakeLists.txt`, no
  Mutect2, no other tool.
* **Shared code:** none.  `gatk_qual_output`, `format_gatk_qual_value` and
  `gatk_compatible_record_text` all live in the **anonymous namespace** of
  `genotype_gvcf_tool.cpp` (`namespace {` at `:48`) and are referenced nowhere
  else in the tree.  `gatk_qual_output()` itself was **not** modified.
  Reachability is also unchanged: `apply_gatk_output_allele_subset()` returns
  early unless `--gatk-compatible-annotations` is set (`:2900`), so the
  non-compat diagnostic profile is bit-for-bit untouched.
* Oracle: **exit 1 before** the fix (2 violations, literal rows in §3),
  **exit 0 after** (17 gated + 1 reported-only case, 0 violations).
  `covered-star-only-record-dense` is now **GATED**, plus one new GATED
  isolating case (§3).
* Step 5: (a) exit 0, (b) exit 0 / 15/15, (c) exit 0 / 19/19,
  (d) **300/300 on both backends**, `REG_EXIT=0`.  Change **kept in the tree**.

## 1. STEP 1 — measured GATK truth

### 1.1 Fixture and commands

The fixture is the oracle's own `STAR_ONLY_COVERED_RECORD`
(`verify_genotype_gvcf_spandel_gatk_oracle.py:243-248`), the same 100 bp `chr1`
of `A` with hand-written `.fai`/`.dict` the oracle already builds
(`write_reference`, `:712-721`), indexed with `gatk IndexFeatureFile` for both
tools:

```
chr1	2	.	AA	A,<NON_REF>	.	PASS	DP=20	GT:DP:AD:PL	0/1:20:0,20,0:100,0,100,100,100,100
chr1	3	.	A	*,<NON_REF>	.	PASS	DP=20	GT:DP:AD:PL	0/1:20:0,20,0:100,100,100,0,100,100
```

```
third_party/jdk17/bin/java -Xmx1g -jar \
  third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
  GenotypeGVCFs -R <chr1.fa> -V <fixture>.g.vcf -O out.vcf \
  --create-output-variant-index false [OPTIONS]
```

The probe is `.diag/qual_inf_probe.py`; its two runs are
`.diag/qual-inf-probe-before.log` and `.diag/qual-inf-probe-after.log`.  The
QUAL column is printed as **raw field bytes** (`split("\t")[5]`), not as a parsed
number, so the token itself is the evidence.

### 1.2 Literal rows and QUAL tokens, per mode (measured)

| fixture (probe case) | OPTIONS | GATK literal rows | GATK QUAL token | native QUAL token — before → after |
| --- | --- | --- | --- | --- |
| `2 AA A` + `3 A *` (**the reported case**) | `--include-non-variant-sites` | pos-2 row + `chr1 3 . A . Infinity . DP=20;MLEAC=.;MLEAF=. GT ./.` | `Infinity` | `0` → **`Infinity`** ✅ |
| same, record `FILTER=LowQual`, `DP=7` | `--include-non-variant-sites` | pos-2 row + `chr1 3 . A . Infinity . DP=7;MLEAC=.;MLEAF=. GT ./.` | `Infinity` | `0` → **`Infinity`** ✅ |
| same, two samples | `--include-non-variant-sites` | `chr1 2 . AA A 190.46 … 0/1 … 0/1` + `chr1 3 . A . Infinity . DP=20;MLEAC=.;MLEAF=. GT ./. ./.` | `Infinity` | `0` → **`Infinity`** ✅ |
| **default mode** (`2 AA A` + `3 A *`) | *(none)* | pos-2 row only — **no pos-3 row at all** | — (only `92.60`) | `92.60` → `92.60` (unchanged) ✅ |
| **star-free** `2 A <NON_REF>` hom-ref block | `--include-non-variant-sites` | `chr1 2 . A . . . DP=20 GT:AD ./.:20` | **`.`** (missing) | `.` → `.` (unchanged) ✅ |
| star-free, flat `PL=0,0,0` | `--include-non-variant-sites` | same row, QUAL `.` | `.` | `.` → `.` (unchanged) ✅ |
| star-free, default | *(none)* | no rows | — | — ✅ |
| finite sibling: `2 A *,G,<NON_REF>` (both ALTs pruned, **same branch**) | `--include-non-variant-sites` | `chr1 2 . A . 127.78 . DP=20;MLEAC=.;MLEAF=. GT ./.` | `127.78` | `127.78` → `127.78` (unchanged) ✅ |
| finite sibling: implausible concrete ALT `A G,<NON_REF>` | `--include-non-variant-sites` | `chr1 2 . A . 127.78 . …` | `127.78` | `127.78` → `127.78` (unchanged) ✅ |

Literal bytes, tab-unfolded:

```
# the reported case, dense mode — GATK and (after the fix) native, byte-identical
chr1	2	.	AA	A	92.60	.	AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63	GT:AD:DP:GQ:PL	0/1:0,20:20:99:100,0,100
chr1	3	.	A	.	Infinity	.	DP=20;MLEAC=.;MLEAF=.	GT	./.

# native BEFORE the fix — the divergence (row 1 was already byte-identical)
chr1	3	.	A	.	0	.	DP=20;MLEAC=.;MLEAF=.	GT	./.
```

Answers to the questions the brief asked:

* **Is it a genuine `+Infinity` double or a saturating conversion?** Genuine.
  The value is `-10 * log10ProbVariantPresent()` and
  `log10ProbVariantPresent()` is `-Infinity` (see §2), so nothing saturates.
* **Is the token a special VCF escape?** No.  It is Java's `%f` rendering of an
  infinite double (§1.3).  `NaN` and `-Infinity` would print `NaN` /
  `-Infinity` by the same mechanism.
* **Was this one defect or a real numeric difference?** The *branch* was
  correct in native; the *value* was being thrown away.  The previous round's
  "GATK `Infinity` vs native `159.55`" observation is, **by code inspection (I
  did not re-measure that specific fixture)**, the same defect seen from a
  fixture where the stale value happened to be non-zero: the `std::isfinite`
  guard skips the assignment entirely, so `record.value->qual` keeps whatever
  the only earlier QUAL assignment stored — `:2170`
  `record.value->qual = gatk_qual_output(cohort.qual)` — which is `0` when the
  AF pass wrote 0 and a finite AF-derived value (`159.55`) when it wrote one.
  So it is **one** defect, and it is *numeric* (a suppressed assignment) with a
  *rendering* half, not a rendering-only change.

### 1.3 Where the `Infinity` comes from, with `file:line`

**GATK (read from `gatk-source/`, then confirmed by measurement).**
`GenotypingEngine.calculateGenotypes()`:

```java
// GenotypingEngine.java:158-163
final double log10Confidence =
            !outputAlternativeAlleles.siteIsMonomorphic || configuration.annotateAllSitesWithPLs
                    ? AFresult.log10ProbOnlyRefAlleleExists() + 0.0 : AFresult.log10ProbVariantPresent() + 0.0;
// Add 0.0 removes -0.0 occurrences.
final double phredScaledConfidence = (-10.0 * log10Confidence) + 0.0;
...
// GenotypingEngine.java:183  -- this is the QUAL the writer sees
builder.log10PError(log10Confidence);
if ( ! passesCallThreshold(phredScaledConfidence) ) {   // :184, only the FILTER decision
    builder.filter(GATKVCFConstants.LOW_QUAL_FILTER_NAME);
}
```

* The site is monomorphic: `calculateOutputAlleleSubset()` sets
  `siteIsMonomorphic &= !(isPlausible && !isSpuriousSpanningDeletion)`
  (`:318`), and an empty output allele set implies every ALT was implausible or
  spurious, so the flag stays `true` and `annotateAllSitesWithPLs` is `false` for
  GenotypeGVCFs — the **second** arm is taken.
* `log10ProbVariantPresent()` is the complement of
  `log10ProbOnlyRefAlleleExists()`.  Here the posterior puts *all* its mass on
  "only the reference allele exists", so that value is `-0.0` and
  `MathUtils.log10OneMinusPow10(-0.0)` returns `Double.NEGATIVE_INFINITY`
  (the sibling native helper `log10_one_minus_pow10` implements exactly this
  rule, including the `log10_value == 0.0 -> -infinity` case).
* `log10PError = -Infinity`, so `getPhredScaledQual() = +Infinity`.  Note the
  distinction from the *finite* sibling: for the `STAR_RECORD` fixture the
  complement is a small negative number (`-12.778`), giving the finite
  `127.78`.

**htsjdk encoder (verified by `javap -c` on the pinned
`gatk-package-4.6.2.0-local.jar`, then confirmed on the pinned JDK 17).**

```
VCFEncoder.vcfFormat:  if (!vc.hasLog10PError()) append(".")
                       else append(formatQualValue(vc.getPhredScaledQual()))
VCFEncoder.formatQualValue(double):
     String.format(Locale.US, "%.2f", qual)      // then strip a trailing ".00"
```

and, measured with `third_party/jdk17/bin/java`:

```
String.format("%.2f", Double.POSITIVE_INFINITY) = "Infinity"
String.format("%.2f", Double.NEGATIVE_INFINITY) = "-Infinity"
String.format("%.2f", Double.NaN)               = "NaN"
```

so the literal `Infinity` token is simply Java's decimal formatting of the
infinite double.  There is nothing VCF-specific about it.

**Native, before the fix** (`genotype_gvcf_tool.cpp:2829-2833` pre-fix):

```cpp
if (record.cohort_quality_available) {
    const auto complement = gatk_log10_one_minus_pow10(record.cohort_log10_p_no_variant);
    if (std::isfinite(complement))
        record.value->qual = gatk_qual_output(-10.0 * complement);
}
```

`gatk_log10_one_minus_pow10` *does* return `-infinity` here — native models GATK
exactly — but the `std::isfinite` guard drops it, so `record.value->qual` keeps
the value written earlier at `:2170`.  Two further gaps stood behind it:
`gatk_qual_output()` (`:253-258`) begins `if (!std::isfinite(value)) return 0.0F;`
and `format_gatk_qual_value()` (`:4290`) returned any non-finite token verbatim,
i.e. HTSlib's `inf`, not `Infinity`.

## 2. Why the branch selection is provably the same

Native reaches `materialize_gatk_monomorphic_ref_call()` **only** from
`apply_gatk_output_allele_subset()` when the surviving output allele set is empty
(`:2941-2951`).  GATK's `siteIsMonomorphic` is `&=`-accumulated over every ALT as
`!(isPlausible && !isSpuriousSpanningDeletion)`, so an empty output set always
leaves it `true`; and `configuration.annotateAllSitesWithPLs` is
`args.annotateAllSitesWithPLs`, which `GenotypeGVCFsEngine.createMinimalArgs()`
never sets (`GenotypeGVCFsEngine.java:373-384`).  Therefore
`materialize_gatk_monomorphic_ref_call()` *is* the `log10ProbVariantPresent()`
arm — no extra condition is needed to decide the branch, and the fix only
completes the value that branch produces.

The cross-validation is the finite sibling: on the same code path with
`cohort_log10_p_no_variant = -12.778…`, native already produced GATK's `127.78`
byte-for-byte (measured before *and* after the fix), which shows the plumbing
from the AF posterior into this assignment is right; only the non-finite case
was discarded.

## 3. STEP 2 — the gate, written and run before the fix

`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`:

* `covered-star-only-record-dense` **promoted from REPORTED ONLY to GATED**
  (`gated: True`, `expect=[GATK_DEL_UPSTREAM_LOCUS_ROW, GATK_STAR_ONLY_DENSE_ROW]`),
  `why` rewritten to cite `GenotypingEngine.java:158-163` + `:183`,
  `MathUtils.log10OneMinusPow10`, and the htsjdk `%.2f` rendering.
* One new **GATED** isolating case, `covered-star-only-record-dense-two-samples`
  (new `HEADER_TWO_SAMPLES` constant, new
  `STAR_ONLY_COVERED_TWO_SAMPLES_RECORD` fixture, two new measured row
  constants): the same locus with two samples.  It separates "the token belongs
  to the monomorphic-confidence branch" from "the token is an artefact of a
  single-sample path"; both its rows are byte-identical after the fix.
* `run_case()` now honours an optional per-case `"header"` key
  (`case.get("header", HEADER)`), so a case can carry its own sample set.  No
  existing case is affected — none of them sets the key.
* The docstring gained a "The dense-mode QUAL token" section quoting
  `:158-163`, `:183`, the `log10OneMinusPow10(0.0)` evaluation and the htsjdk
  encoder disassembly, so the gate is self-describing.
* The finite sibling of the branch was **already gated** as
  `include-non-variant-sites` (`GATK_DEFAULT_ROW`, `127.78`), which is what makes
  "monomorphic ⟹ `Infinity`" falsifiable.
* No existing case's expectation changed.

### 3.1 Exit status and literal rows

**BEFORE the fix** — **exit 1**, `"status": "divergence"`, **2 violations**
(`.diag/qual-inf-gate-before.log`), both exactly the QUAL token and nothing else:

```
[covered-star-only-record-dense] row 1 is not byte-identical:
  GATK  ='chr1\t3\t.\tA\t.\tInfinity\t.\tDP=20;MLEAC=.;MLEAF=.\tGT\t./.'
  NATIVE='chr1\t3\t.\tA\t.\t0\t.\tDP=20;MLEAC=.;MLEAF=.\tGT\t./.'
[covered-star-only-record-dense-two-samples] row 1 is not byte-identical:
  GATK  ='chr1\t3\t.\tA\t.\tInfinity\t.\tDP=20;MLEAC=.;MLEAF=.\tGT\t./.\t./.'
  NATIVE='chr1\t3\t.\tA\t.\t0\t.\tDP=20;MLEAC=.;MLEAF=.\tGT\t./.\t./.'
```

The other 15 cases (including the finite `127.78` sibling and all eleven
pre-existing gates) were green **before** the fix, so the two violations isolate
the infinite value alone.  The GATK-truth assertion (`expect`) also passed before
the fix — the measured rows did not move.

**AFTER the fix: exit 0**, `"status": "pass"`, 0 violations
(`.diag/qual-inf-gate-after.log`).  Case census: **17 gated, 1 reported-only**.

| case | gated | GATK rows | native rows |
| --- | --- | --- | --- |
| the 15 previously gated cases | yes | — | byte-identical |
| **`covered-star-only-record-dense`** | **yes (was reported-only)** | 2 | 2, incl. `Infinity` |
| **`covered-star-only-record-dense-two-samples`** | **yes (new)** | 2 | 2, incl. `Infinity` |
| `unemitted-upstream-deletion-star-plus-concrete-alt` | no | 1 | 1 (ALT differs — unrelated, §6) |

## 4. STEP 3 — the fix

Three lines of behaviour, both in `fastgatk-native/src/genotype_gvcf_tool.cpp`;
full patch `.diag/qual-inf-fix.diff`, `git diff --stat` =
**2 files changed, 122 insertions(+), 8 deletions(-)** (oracle +106/-… ,
tool +24/-…):

```diff
@@ materialize_gatk_monomorphic_ref_call  (was :2829-2833)
         const auto complement = gatk_log10_one_minus_pow10(record.cohort_log10_p_no_variant);
         if (std::isfinite(complement))
             record.value->qual = gatk_qual_output(-10.0 * complement);
+        else if (std::isinf(complement) && complement < 0.0)
+            // GenotypingEngine.java:158-163 assigns that same
+            // log10ProbVariantPresent() straight into the record
+            // (builder.log10PError(log10Confidence) at :183), so when the
+            // posterior puts ALL of its mass on "no variant present", ...
+            record.value->qual = std::numeric_limits<float>::infinity();
     }

@@ format_gatk_qual_value  (was :4290-4291)
 std::string format_gatk_qual_value(const std::string& value) {
     if (value.empty() || value == ".") return value;
+    // htsjdk's VCFEncoder.formatQualValue() is Java's `%.2f`, which prints an
+    // infinite double as "Infinity"; HTSlib's kputd() renders the same float
+    // through C's "%.g", i.e. "inf"/"-inf" (htslib kstring.c:38, reached from
+    // vcf.c:4109), so translate the spelling.
+    if (value == "inf") return "Infinity";
+    if (value == "-inf") return "-Infinity";
     try {
```

Why this is the whole fix and nothing more:

* `else if (std::isinf(complement) && complement < 0.0)` is deliberately narrow:
  a **NaN** complement (`log10Confidence > 0`, a numerical artefact) keeps the
  previous leave-unchanged behaviour, exactly as the old `std::isfinite` guard
  did, so nothing about the NaN case changes.  Only the case that GATK itself
  turns into `+Infinity` is affected.
* `gatk_qual_output()` is **not** touched, so its other two callers
  (`:2035` HaplotypeCaller-posterior path, `:2170` cohort AF path) are
  unaffected.
* The infinite value is kept in the `bcf1_t` float QUAL and rendered at the
  compat text boundary, which is where the htsjdk spelling belongs; HTSlib's
  `kputd()` (third_party/htslib-build/htslib-src/kstring.c:38, called from
  `vcf.c:4109`) would otherwise emit `inf`.
* `std::numeric_limits` was already used in this file (`:257`), so no include
  changed.
* Both serialization paths share `apply_gatk_output_allele_subset()`, and both
  encode stages call `gatk_compatible_record_text()` (`:5711`, `:6420`).
  Measured directly: the aggregate and `--stream-by-locus` invocations now both
  emit `chr1 3 . A . Infinity . …`.

Not restructured: the genotyping engine, the AF/PL Kokkos kernels and their
ABIs, the union/merge stages, the deletion-span construction, the
output-allele-subset predicate fixed in the previous round, and every other
tool.

## 5. STEP 4 — stale assertions

**No registered test pinned the old token; no test line needed correcting.**
Verified rather than assumed:

| File / check | Why it could have been stale | Outcome |
| --- | --- | --- |
| whole `fastgatk-native/scripts/*.py` | a test could pin `Infinity` or the old `0` | `grep` for the REF-only row signature `MLEAC=.;MLEAF=.` returns exactly 5 hits: `verify_genotype_gvcf.py:417` (`127.78`, finite, unchanged) and four inside the oracle itself.  No hit carries a QUAL of `0`. |
| `verify_genotype_gvcf.py:403-419` (STAR_RECORD `*,G` dense, compat) | it asserts the literal row `… A . 127.78 . …` — the *finite* sibling of the same branch | the fix does not fire on a finite complement; unchanged and passing (`B1_EXIT=0`) |
| `verify_genotype_gvcf.py:422-435` (`--gp-qual` dense, compat) | same row | unchanged and passing |
| `verify_genotype_gvcf.py:831-870` (dense REF-only block, **non-compat** profile, asserts ALT/GT/PL) | `materialize_gatk_monomorphic_ref_call` is unreachable without `--gatk-compatible-annotations` | untouched by construction; unchanged and passing |
| `verify_genotype_gvcf.py:870-879` (dense aggregate vs `--stream-by-locus`, non-compat, byte-equality) | both sides could have diverged | untouched; unchanged and passing |
| the 15 `ctest -R genotype-gvcf` tests, incl. `include-non-variant-gatk-oracle` and `multisample-reference-confidence-gatk-oracle` | they exercise the compat dense path | 15/15 passed after the fix (§6.1) |
| C++ / ctest sources | a unit test could pin the QUAL assignment | `grep` for `materialize_gatk_monomorphic_ref_call`, `format_gatk_qual_value` and `gatk_qual_output` over `*.cpp/*.hpp/*.py/*.sh/*.txt/*.cmake` → only `genotype_gvcf_tool.cpp`; no unit test pins them |

The only test-script change in this round is the oracle itself (§3).

## 6. STEP 5 — gate results (mandatory order)

| Step | Command | Result |
| --- | --- | --- |
| a | `python3 fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` | **exit 0** (17 gated, 0 violations); **exit 1 before the fix** (2 violations, §3.1) |
| b1 | `python3 fastgatk-native/scripts/verify_genotype_gvcf.py` | `B1_EXIT=0` — `{"status": "pass", "output_records": 1}` |
| b2 | `ctest --test-dir fastgatk-native/build -R 'genotype-gvcf' -j 8` | `B2_EXIT=0` — **15/15 passed**, incl. `fastgatk-genotype-gvcf-spandel-gatk-oracle` (210.69 s) |
| c | the 17-name strict-gate filter from the task | `C_EXIT=0` — **19/19 matched tests passed**, 294.78 s |
| d | `fastgatk-native/scripts/run_regression.sh --label qual-inf` | `REG_EXIT=0` — **300/300 omp**, **300/300 serial** |

Binaries rebuilt in both trees **after** the last source edit and before every
gate above: `md5sum .diag/qual-inf-binaries.md5` →
`fbf2ae1aa76fc888058a06a42fa52f77` (omp),
`34b502faa51edbda69587d5aacb34031` (serial).

Logs: `.diag/qual-inf-{gate-before,gate-after,verify-gvcf,ctest-genotype-gvcf,strict-gates,probe-before,probe-after,regression}.log`,
build logs `.diag/qual-inf-build-{omp,serial}.log`, patch
`.diag/qual-inf-fix.diff`, probe script `.diag/qual_inf_probe.py`.

### 6.1 Strict gates (c) and the double-backend regression (d)

* **(c)** ran after the oracle revision was frozen: `C_EXIT=0`, **19/19 passed**
  (`.diag/qual-inf-strict-gates.log`).  The filter matched
  `fastgatk-hc-ploidy-window-invariance-gatk-oracle` and a second
  `mutect2-recheck` test in addition to the named 17; the
  `spandel-gatk-oracle` entry (`fastgatk-genotype-gvcf-spandel-gatk-oracle`,
  222.62 s) is the promoted gate itself and ran the final revision.
* **(d)** `fastgatk-native/scripts/run_regression.sh --label qual-inf`
  (omp + serial, `FASTGATK_REQUIRE_GATK_ORACLE=1` by default, ctest parallelism
  8), `REG_EXIT=0`:

  ```
  | 后端   | 构建目录                               | 结果 | 通过/总数 | 耗时  |
  | omp    | OpenMP (fastgatk-native/build)         | 通过 | 300/300   | 1238.94 sec |
  | serial | Serial (fastgatk-native/build-serial)  | 通过 | 300/300   | 1269.46 sec |
  ```

  Evidence block `.diag/regression/20260912-045714/{omp.log,serial.log}`; the
  runner's summary names the same two trees and reports the staleness check as
  clean.

**Nothing was edited after the suite run, so no rebuild or re-run was needed.**
`md5sum -c .diag/qual-inf-binaries.md5` re-confirms
`fbf2ae1aa76fc888058a06a42fa52f77` (omp) and
`34b502faa51edbda69587d5aacb34031` (serial).  Timestamps, all **before** the
suite window 04:57:14 → 05:18:24:
`genotype_gvcf_tool.cpp` 04:45:11 < omp binary 04:45:29 < serial binary
04:45:42, and the oracle `verify_genotype_gvcf_spandel_gatk_oracle.py` 04:42:05.
The only file touched afterwards is this report, which is not a test input.

## 7. What is shared, and what else it affects

* All three helpers are in the **anonymous namespace** of a single translation
  unit (`namespace {` at `genotype_gvcf_tool.cpp:48`); `grep` over the whole tree
  finds no other reference to `format_gatk_qual_value` or
  `gatk_compatible_record_text`.  The nearest look-alike,
  `variants_to_table_tool.cpp:681 format_qual()`, is a different function and was
  not touched.
* `format_gatk_qual_value()` is reached only from `gatk_compatible_record_text()`,
  which is reached only on the `--gatk-compatible-annotations` path: two call
  sites, both in the genotype-gvcf streaming encode stages.
* Consequently the changed surface is exactly "GenotypeGVCFs with
  `--gatk-compatible-annotations`, dense mode, a locus that turned monomorphic".
  The non-compat diagnostic profile, every other tool, Mutect2, the Kokkos
  kernels and the registered ctest list are unaffected — confirmed by (b), (c)
  and (d) rather than by inspection alone.

## 8. What remains unproven / open

* **New residual, measured, not fixed (unrelated to QUAL):** the dense-mode
  REF-only row for a **non-PASS** source record leaf diverges in FILTER.  GATK
  writes `chr1 3 . A . Infinity . DP=7;MLEAC=.;MLEAF=. GT ./.` while native
  writes `chr1 3 . A . Infinity RGQ DP=7;…` — the QUAL now matches, the FILTER
  does not (`LowQual` in the input, `RGQ` in native's output).  *Speculation:*
  native appears to carry the source record's FILTER **id** into a header whose
  FILTER dictionary has different indices.  I did not investigate it and did not
  gate it, because it is a different column; the DEFAULT-mode statement of the
  same fixture is already gated as `covered-star-only-record-non-pass` and
  passes.
* *Speculation:* the `-inf → "-Infinity"` branch of the formatter is not
  reachable in this tool (a `-Infinity` QUAL would need `log10PError = +Infinity`,
  which the code cannot produce).  It is included only so that Java's `%.2f`
  spelling is complete; the measured, reachable branch is `inf → "Infinity"`.
* The residual **input-records-vs-emitted-alleles** divergence is untouched and
  still open as
  `unemitted-upstream-deletion-star-plus-concrete-alt`; nothing in this round
  changes it.
* **Untested shapes** of the infinite QUAL: multi-shard input, `--max-alternate-alleles`
  reduction together with a monomorphic dense locus, and
  `--force-output-intervals`.  *Speculation:* the rule is per-locus and
  allele-count independent, and the two-sample case supports that, but these
  were not measured.
* Still-open pre-existing items from earlier rounds are untouched:
  `FILTER=LowQual` is not implemented in native GenotypeGVCFs, `-all-sites` is
  not accepted as a short alias, the `genotypeLikelihoods == null` arm of
  PREFER_PLS is not implemented, and a `./.` source call may still produce a
  record where GATK writes none.
* **Unproven by construction:** the oracle asserts byte-identical data rows on
  the fixtures it pins.  It is not a proof that native's `cohort_log10_p_no_variant`
  always equals GATK's `AFresult.log10ProbOnlyRefAlleleExists()`; the finite
  `127.78` sibling is the cross-check that it does on this branch, and the
  `Infinity` case now agrees too, but both are point measurements, not a proof.

## 9. Tree state

* Change **kept in the tree** (no commit, no branch, no branch switch).
* Modified — exactly two files, `git status --short`:

  ```
   M fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py
   M fastgatk-native/src/genotype_gvcf_tool.cpp
  ```

  `git diff --stat` = **2 files changed, 122 insertions(+), 8 deletions(-)**.
  Scratch lives under `.diag/` only, which is git-ignored, so it does not appear
  in `git status`.
* Not touched: `fastgatk-native/CMakeLists.txt`, every root `*.md`, Mutect2,
  every other tool, and every registered test script other than the oracle.
* Both build trees were rebuilt after the last **source** edit and before every
  gate in §6; the oracle's last edit preceded steps (a), (b), (c) and (d).
* Binaries are byte-identical to the ones the suite ran
  (`md5sum -c .diag/qual-inf-binaries.md5` → both `OK`).
