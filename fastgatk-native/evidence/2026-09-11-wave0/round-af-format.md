# Round: `AF`/`MLEAF` zero-value formatting (GATK vs native)

Date: 2026-09-11 (local). Repo `/home/turing-agents/Documents/fast-gatk`.
Observation recovered from `.diag/round-fixture-construction.md` section **D.1**
(“`INFO AF`/`MLEAF` zero formatting on multi-ALT records”).

Toolchain used: pinned GATK 4.6.2.0
(`third_party/jdk17/bin/java -Xmx1g -jar third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar HaplotypeCaller`),
native OpenMP `fastgatk-native/build/fastgatk-hc-call`, native serial
`fastgatk-native/build-serial/fastgatk-hc-call`.

---

## 0. Bottom line

| item | result |
|---|---|
| Divergence reproduced? | **YES**, exactly as reported: GATK `AF=0.500,0.00` vs native `AF=0.500,0.000` (same for `MLEAF`) |
| Root cause | `fastgatk-native/src/hc_call.cpp` `format_allele_frequency()` formatted *every* value with three fixed places (except ~1.0), so an exact zero became `0.000`; GATK formats each list element through htsjdk's magnitude-driven `VCFEncoder.formatVCFDouble`, whose zero is the literal `"0.00"` |
| GATK rule established from | **decompiled bytecode of the pinned jar** (`javap -p -c htsjdk.variant.vcf.VCFEncoder`), corroborated by **GATK source** (`GenotypingEngine.composeCallAttributes` puts a `List<Double>` into `MLEAF`) and by observed GATK rows |
| Fix | one function, one file: `hc_call.cpp` (`format_allele_frequency`) — 27 insertions / 4 deletions including the comment |
| Tree state | **CONTAINS** the change (`fastgatk-native/src/hc_call.cpp`, plus the new oracle `fastgatk-native/scripts/verify_hc_af_zero_format_gatk_oracle.py`) |
| Oracle | `fastgatk-native/scripts/verify_hc_af_zero_format_gatk_oracle.py`: **exit 1 before the fix**, exit 0 after (OpenMP and serial) |
| Gate chain | oracle red → fix → oracle green; `ctest -R 'gatk-oracle\|gate-oracle'` **159/159**; the 3 named gates that regex misses **3/3**; `verify_hc_alleles_gatk_oracle.py` pass; regression **289/289 on omp and 289/289 on serial** |
| Nothing else touched | no `CMakeLists.txt`, no root `*.md`, no existing oracle, no EventMap/PairHMM/AF/prior machinery, no Mutect2 |

---

## 1. The reproduced divergence (literal rows)

Fixture (built by the new oracle; identical synthetic reference to the fixture
library: 1500 bp `chr1`, `random.seed(11)`, 30 × 300M reference reads from 440
step 3, 12 × 300M reads carrying `chr1:698 C>A` from 440 step 3), plus a
bgzip/tabix’d single-record feature VCF `chr1 698 . C A,G,T . PASS .` fed to
both callers with `--alleles … --max-alternate-alleles 3 -L chr1:500-780
--min-pruning 1`.

`--case multialt-3alt-ploidy2` (`.diag/af_zero_prefix_strict.log`, pre-fix):

```
GATK   chr1 698 . C A,G,T 403.64 . AC=1,0,0;AF=0.500,0.00,0.00;AN=2;BaseQRankSum=0.000;DP=42;ExcessHet=0.0000;FS=0.000;MLEAC=1,0,0;MLEAF=0.500,0.00,0.00;MQ=60.00;MQRankSum=0.000;QD=9.61;ReadPosRankSum=-2.995;SOR=0.162  GT:AD:DP:GQ:PL  0/1:30,12,0,0:42:99:411,0,1217,501,1253,1754,501,1253,1754,1754
NATIVE chr1 698 . C A,G,T 403.64 . AC=1,0,0;AF=0.500,0.000,0.000;AN=2;BaseQRankSum=0.000;DP=42;ExcessHet=0.0000;FS=0.000;MLEAC=1,0,0;MLEAF=0.500,0.000,0.000;MQ=60.00;MQRankSum=0.000;QD=9.61;ReadPosRankSum=-2.995;SOR=0.162  GT:AD:DP:GQ:PL  0/1:30,12,0,0:42:99:411,0,1217,501,1253,1754,501,1253,1754,1754
```

Field-level diff of that pair — **only** the two zero entries differ:

```
INFO AF        '0.500,0.00,0.00'   -> '0.500,0.000,0.000'
INFO MLEAF     '0.500,0.00,0.00'   -> '0.500,0.000,0.000'
```

Every case of the pre-fix diagnostic run, with the native text in the last
column (measured unless marked otherwise):

| case | GATK | native (pre-fix) |
|---|---|---|
| `multialt-3alt-ploidy2` (A,G,T) | `AF=0.500,0.00,0.00` / `MLEAF=0.500,0.00,0.00` | `0.500,0.000,0.000` |
| `multialt-2alt-ploidy2` (A,G) | `AF=0.500,0.00` | `0.500,0.000` |
| `multialt-2alt-ploidy3` (A,G) | `AF=0.333,0.00` | `0.333,0.000` |
| `multialt-2alt-ploidy4` (A,G) | `AF=0.250,0.00` | `0.250,0.000` |
| `biallelic-forced-control` (A,T) | `AF=0.500,0.00` | `0.500,0.000` |
| `allzero-forced-control` (A,G,T, no variant reads) | `AF=0.00,0.00,0.00` | `0.00,0.00,0.00` (**already agreed**) |
| `noalleles-biallelic-snp-control` (no `--alleles`) | `AF=0.500` | `0.500` (**already agreed**) |
| `multialt-2alt-ploidy3-span-del` | `AF=0.333` (biallelic: the forced record has no zero entry) | `0.333` (**no zero to format**) |
| `gvcf-ploidy3-span-del-mleaf` (probe) | `MLEAF=0.333,0.00` / `MLEAF=0.00,0.333,0.00` | identical (**gVCF writer already conformant**) |
| `multialt-zero-first-ploidy2` (reads on G) — added after the fix | `AF=0.500,0.00,0.00` on ALT `G,A,T` | `0.000,0.500,0.000` on ALT `A,G,T` (**inferred, see §5**) |

The two already-agreeing zero rows are why the defect survived earlier rounds:
native's ordinary writer carried a *special case* for the all-zero forced-allele
record (`forced_hom_ref_feature && ac[i] == 0 ? "0.00" : format_allele_frequency(...)`,
`hc_call.cpp:3660/3683`) and the gVCF writers carry `count == 0 ? "0.00" : ...`
(`hc_call.cpp:4664/5030`).  The divergence is exactly the gap between those two
special cases: **a zero entry in a record that also has a non-zero entry**.

### 1.1 Adjacent observation (not this defect, not gated)

`multialt-zero-first-ploidy2` shows a *separately tracked* ordering difference
that the fixture library's canonicalizer absorbs: when the read support is on a
non-first forced allele, GATK emits the ALT list in **haplotype order**
(`G,A,T`) while native emits sorted concrete alleles (`A,G,T`).  The oracle
canonicalizes both sides identically (the library permutes every allele-indexed
field with the ALT list) and reports the normalization as a note; it does **not**
gate ALT order.  `--alleles`-driven ALT ordering is a distinct defect class (the
earlier round saw the same thing with insertion ALTs, its section D.4).

### 1.2 Where it is *not*

* The ploidy-3 `-ERC BP_RESOLUTION` spanning-deletion record, which the earlier
  round’s D.1 note also blamed, **already agrees** — its `MLEAF=0.333,0.00,0.00`
  comes from the gVCF writer’s `count == 0` special case.  Probe
  `gvcf-ploidy3-span-del-mleaf` (281 rows per side, ALT order canonicalized)
  is byte-identical pre- **and** post-fix, and it also shows the earlier round’s
  `*`-prior QUAL fix is still in place (`chr1:698 QUAL 368.90` on both sides,
  formerly `368.90 / 369.37`).
* `AF`/`MLEAF` are the only fields affected: `QUAL`, `QD`, every other INFO
  annotation, every FORMAT field and `GT` are byte-identical in all cases; the
  oracle gates **all** of them, not just the two fields.

---

## 2. GATK’s formatting rule, and how it was established

### 2.1 The rule (per value, magnitude-driven)

`htsjdk.variant.vcf.VCFEncoder.formatVCFDouble(double)`, decompiled from the
pinned `gatk-package-4.6.2.0-local.jar` with
`third_party/jdk17/bin/javap -p -c -cp <jar> htsjdk.variant.vcf.VCFEncoder`
(the JVM bytecode, reproduced verbatim):

```
 0: dload_0 / dconst_1 / dcmpg / ifge 43   ->  d >= 1.0       -> format = "%.2f"
 6: dload_0 / ldc 0.01  / dcmpg / ifge 36   ->  d >= 0.01      -> format = "%.3f"
14: dload_0 / Math.abs / ldc 1.0E-20 / dcmpl / iflt 32
                                            ->  abs(d) >= 1e-20 -> format = "%.3e"
32: ldc "0.00" / areturn                    ->  otherwise      -> return "0.00"
46: String.format(Locale.US, format, d)
```

So the serialization is **not** “one precision for the whole list” and **not**
“trim each value to a minimum of two decimals” (the earlier round’s guess in
`.diag/round-fixture-construction.md` D.1); it is a per-element decision on the
value’s own magnitude.  For a frequency list this means:

| value | text |
|---|---|
| `0.0` (and `abs(d) < 1e-20`) | `0.00` (literal) |
| `0.01 <= d < 1.0` | three fixed places, e.g. `0.500`, `0.333`, `0.250` |
| `d >= 1.0` | two fixed places, e.g. `1.00` |
| `1e-20 <= abs(d) < 0.01` | scientific, e.g. `5.000e-03` |

which is exactly why one record can legally carry `0.500,0.00`.

### 2.2 Why this function applies to `AF`/`MLEAF` and not to `FS`/`QD`

* `htsjdk`’s `VCFEncoder.formatVCFField(Object)` sends every `Double` (and every
  element of a `double[]`/`List<Double>`) through `formatVCFDouble`; `String`
  values are written verbatim (same bytecode dump).
* GATK hands `MLEAF` to the writer as a **`List<Double>`**:
  `GenotypingEngine.composeCallAttributes` →
  `calculateMLEAlleleFrequencies(alleleCountsofMLE, genotypes)` whose body is
  `alleleCountsofMLE.stream().map(AC -> Math.min(1.0, (double) AC / AN))`
  (`gatk-source/.../walkers/genotyper/GenotypingEngine.java:439-443, 471-474`);
  `AF` is filled the same way from the chromosome counts.
* GATK’s *other* HC INFO annotations are pre-formatted **`String`s**:
  `QualByDepth.annotate` returns `String.format("%.2f", QD)`
  (`gatk-source/.../annotator/QualByDepth.java:91`) and `FisherStrand` returns
  `String.format("%.3f", …)` (`FisherStrand.java:77-79`).  That is why GATK
  prints `FS=0.000`/`BaseQRankSum=0.000` (a literal zero as a *string*) **and**
  `AF=0.500,0.00` in the same row — a raw `Double` 0.0 would have gone through
  `formatVCFDouble` and printed `0.00`, which is precisely the divergence found
  on `AF`/`MLEAF`.

This is the difference between the two classes, and both halves are checkable:
the bytecode says the rule, the source says which fields are `Double`s, and the
observed GATK rows match both.

Evidence strength, stated precisely: `MLEAF` being a `List<Double>` is **read
directly from GATK source** (`GenotypingEngine.java:439-443, 471-474`); `AF` is
filled by htsjdk's `VariantContextUtils.calculateChromosomeCounts`, which the
GATK source calls (`GATKVariantContextUtils.java:1243, 1840`) but which lives in
the jar, so its Java type was **not** decompiled this round — the evidence for
`AF` is that its observed text (`0.500` next to the literal `0.00`) is exactly
`formatVCFDouble`’s output and is inconsistent with any fixed-precision
formatter.  Both fields go through the same htsjdk encoder either way.

---

## 3. The diff (minimal fix)

`fastgatk-native/src/hc_call.cpp`, function `format_allele_frequency` only:

```diff
+// ``AF``/``MLEAF`` are Number=A Float INFO values carried as raw Java doubles:
+// htsjdk's ``VariantContextUtils.calculateChromosomeCounts`` fills ``AF`` and
+// GATK's ``GenotypingEngine.composeCallAttributes`` fills ``MLEAF`` from
+// ``calculateMLEAlleleFrequencies``, whose body is
+// ``alleleCountsofMLE.stream().map(AC -> Math.min(1.0, (double) AC / AN))``
+// (a ``List<Double>``), so htsjdk's encoder sees a Double per element.  The
+// encoder is ``htsjdk.variant.vcf.VCFEncoder.formatVCFDouble(double)``, whose
+// rule (decompiled from the pinned gatk-package-4.6.2.0-local.jar:
+// ``d >= 1.0`` -> "%.2f"; ``d >= 0.01`` -> "%.3f"; ``abs(d) >= 1e-20`` ->
+// "%.3e"; otherwise the literal "0.00") is *per value and magnitude-driven*.
+// That is why one list legitimately mixes "0.500" (interior value) with the
+// literal "0.00" (an exact zero), which the previous value-independent
+// three-place formatter could not express.
 std::string format_allele_frequency(const double frequency) {
-    std::ostringstream value;
-    if (std::abs(frequency - 1.0) < 1.0e-12)
+    if (frequency >= 1.0) {
+        std::ostringstream value;
         value << std::fixed << std::setprecision(2) << frequency;
-    else
+        return value.str();
+    }
+    if (frequency >= 0.01) {
+        std::ostringstream value;
         value << std::fixed << std::setprecision(3) << frequency;
-    return value.str();
+        return value.str();
+    }
+    if (std::abs(frequency) >= 1.0e-20) {
+        std::ostringstream value;
+        value << std::scientific << std::setprecision(3) << frequency;
+        return value.str();
+    }
+    return "0.00";
 }
```

`git diff --stat`: `fastgatk-native/src/hc_call.cpp | 27 +++++++++++++++++++++++----`
(1 file changed, 27 insertions, 4 deletions).

Behaviour change, exhaustively (every call site is `AF`/`MLEAF` in the ordinary
writer and `MLEAF` in the two gVCF writers; the gVCF sites only reach the
function for a non-zero count, so they are untouched):

* `d == 0.0` → `0.00` instead of `0.000` — **the fix**;
* `1e-20 <= d < 0.01` → `5.000e-03` instead of `0.005` (see §6, unexercised);
* `d >= 1.0` → `1.00` (unchanged; `AF = AC/AN` is exactly `1.0` iff `AC == AN`);
* `0.01 <= d < 1.0` → unchanged (`0.500`, `0.333`, `0.250`);
* everything else unchanged.

The redundant `count == 0 ? "0.00"` / `forced_hom_ref_feature && ac[i] == 0`
special cases were deliberately **left in place**: they now produce exactly the
same text as the function, and removing them would widen the diff into three
other writer sites for no behavioural gain.

---

## 4. The oracle

`fastgatk-native/scripts/verify_hc_af_zero_format_gatk_oracle.py` (new,
strict, `main() -> int`, `argparse --native`/`FASTGATK_HC_BINARY`,
`tempfile.TemporaryDirectory`, `--expect-divergence`, final JSON status line —
the style of `verify_hc_span_del_qual_oracle.py`).  It builds the reference,
BAMs and `--alleles` feature VCFs itself and compares **every data row field by
field** (REF/ALT/QUAL/FILTER, each INFO sub-field, FORMAT keys and values) after
the fixture library’s documented ALT-order canonicalization.

Cases and gates:

| case | gated | purpose |
|---|---|---|
| `multialt-3alt-ploidy2` | yes | 3 forced ALTs, `AF=0.500,0.00,0.00` — the reported row |
| `multialt-2alt-ploidy2` | yes | the reported fixture’s exact shape (2 ALTs) |
| `multialt-2alt-ploidy3` | yes | `0.333,0.00` at AN=3 |
| `multialt-2alt-ploidy4` | yes | `0.250,0.00` at AN=4 |
| `multialt-zero-first-ploidy2` | yes | reads on the **second** forced allele: the zero is the first list entry (and the ALT order differs between callers — see §1.1) |
| `multialt-2alt-ploidy3-span-del` | yes | forced alleles on the span-del pileup (comes out biallelic: no zero entry) |
| `biallelic-forced-control` | yes | forced A,T with reads on A — a zero in the list, no third non-zero entry |
| `allzero-forced-control` | yes | `0.00,0.00,0.00` — the shape that already agreed; a blanket `%.2f` fix would also pass this |
| `noalleles-biallelic-snp-control` | yes | the canonical single-ALT record (`0.500`), i.e. the main serialization path |
| `gvcf-ploidy3-span-del-mleaf` | probe | gVCF `MLEAF` zero path, 281 rows; reported, not gated (its record needs ALT-order canonicalization) |

Ten cases total: nine gated, one probe.  Logs:
`.diag/af_zero_prefix.log` (pre-fix diagnostic, 9 cases),
`.diag/af_zero_prefix_strict.log` (pre-fix strict, exit 1),
`.diag/af_zero_final.log` (post-fix strict, 10 cases, exit 0),
`.diag/af_zero_serial.log` (post-fix against the serial binary, exit 0),
`.diag/af_zero_zerofirst.log` (the zero-first case alone).

Fixture-validity guard: a gated case that claims a zero entry must actually
carry, on the GATK side, a numerically-zero `Number=A` entry spelled `0.00`
next to a non-zero entry **in a byte-identical row**; otherwise the oracle fails
even if the rows happen to match, so the gate cannot silently stop exercising
the rule.

---

## 5. Gate results

Sequence as prescribed: oracle red first, then fix, then the full gate chain.

| # | command | result | exit | log |
|---|---|---|---|---|
| 1 | `python3 fastgatk-native/scripts/verify_hc_af_zero_format_gatk_oracle.py --case multialt-3alt-ploidy2 --case multialt-2alt-ploidy2` **before** the fix | 4 violations, `"status": "fail"` | **1 (red, as required)** | `.diag/af_zero_prefix_strict.log` |
| 1b | full oracle `--expect-divergence` **before** the fix (9 cases) | 5 zero-entry cases red, 4 green | 0 (diagnostic) | `.diag/af_zero_prefix.log` |
| 2 | rebuild OpenMP: `cmake --build fastgatk-native/build --target fastgatk-hc-call -j 16` | built | **0** | — |
| 3 | full oracle, OpenMP, **after** the fix (10 cases) | `"status": "pass"`, every row byte-identical | **0** | `.diag/af_zero_final.log` |
| 4 | rebuild serial: `cmake --build fastgatk-native/build-serial --target fastgatk-hc-call -j 8` | built | **0** | — |
| 5 | oracle against `fastgatk-native/build-serial/fastgatk-hc-call` (3 cases) | pass | **0** | `.diag/af_zero_serial.log` |
| 6 | `ctest --test-dir fastgatk-native/build -R 'gatk-oracle\|gate-oracle' -V` | **100% tests passed, 0 tests failed out of 159** (3135 s) | **0** | `.diag/af_zero_ctest.log` |
| 7 | `python3 fastgatk-native/scripts/verify_hc_alleles_gatk_oracle.py` | `"status": "pass"` | **0** | `.diag/af_zero_alleles_oracle.log` |
| 8 | `fastgatk-native/scripts/run_regression.sh --label 'AF/MLEAF 零值格式化修复后'` | **omp 289/289, serial 289/289** (1135 s / 1136 s) | **0** | `.diag/af_zero_regression.out`, `.diag/regression/20260911-100546/{omp,serial}.log` |
| 9 | full oracle re-run with the final case metadata (10 cases) | `"status": "pass"`, 10/10 rows identical | **0** | `.diag/af_zero_final.log` |
| 10 | `ctest -R '^fastgatk-hc-arbitrary-ploidy-span-del-prior-oracle$'`, `…-polyploid-gvcf-span-del-prior-oracle$`, `…-spanning-prior-genotype-gq-oracle$` (run one by one, with `FASTGATK_HC_BINARY`/`FASTGATK_REQUIRE_GATK_ORACLE`) | 3 × **100% tests passed** (52.5 s / 59.0 s / 56.6 s) | **0** | job `bash-101` output |

The regression evidence block (`.diag/af_zero_regression.out`) records
`git 951d265` with 2 uncommitted changes (the modified `hc_call.cpp` and the new
oracle script) — i.e. the run measured exactly this working tree, and **no test
failed on either backend**, so the fix does not regress the 289-test suite.

**Caveat found while gating (worth fixing in a future round, not fixable here
because `fastgatk-native/CMakeLists.txt` is off-limits):** only **5** of the
**8** named strict gates are actually selected by
`ctest -R 'gatk-oracle|gate-oracle'`.  Tests #88/#89/#90
(`fastgatk-hc-arbitrary-ploidy-span-del-prior-oracle`,
`fastgatk-hc-polyploid-gvcf-span-del-prior-oracle`,
`fastgatk-hc-spanning-prior-genotype-gq-oracle`) do **not** contain the literal
substring `gatk-oracle` (their names end in `-prior-oracle` / `-gq-oracle`), so
the prescribed regex silently skips the three fixtures landed by the previous
round.  They were therefore executed explicitly (row 10) and all three passed,
so the earlier fix is **not** regressed by this round either.  The remaining
five named gates are inside the 159-test selection.

**What was measured red vs inferred red.**  The 5 cases `multialt-3alt-ploidy2`,
`multialt-2alt-ploidy2`, `multialt-2alt-ploidy3`, `multialt-2alt-ploidy4` and
`biallelic-forced-control` were observed red on the pre-fix binary (table in
§1).  `multialt-zero-first-ploidy2` was added to the oracle *after* the fix
landed (to cover the zero-in-first-position shape), so its pre-fix redness is
**inferred, not measured**: the pre-fix `format_allele_frequency` printed three
fixed places for every value not within `1e-12` of `1.0`, hence `0.0 -> "0.000"`,
and GATK's row for that fixture is `0.500,0.00,0.00`; the case is measured green
post-fix.

---

## 6. What remains unproven / speculation

* **Speculation-free but unexercised**: the `1e-20 <= abs(d) < 0.01` → `%.3e`
  branch is copied from the decompiled bytecode and is **not** exercised by any
  fixture in this round: reaching it needs an `AF` below `0.01`, i.e. a single
  sample with ploidy > 100 (`AF = AC/AN`), which no oracle builds.  If a future
  round shows GATK printing `0.005` there instead of scientific notation, this
  branch is the line to change; the bytecode says otherwise, but the bytecode is
  a source-level fact, not an observed GATK row.
* **Not proven**: that htsjdk’s `formatVCFDouble` is the *only* formatter for
  `AF`/`MLEAF` in every GATK output path.  It is proven for HaplotypeCaller
  ordinary and `-ERC BP_RESOLUTION` gVCF output (observed rows + bytecode +
  the `List<Double>` assignment).  Other tools that write `AF` were not touched
  or tested this round.
* **Not addressed** (adjacent, reported only): `fastgatk-native/src/genotype_gvcf_tool.cpp`
  and the other tools with their own float serializers
  (`filter_mutect_tool.cpp`, `mutect2_tool.cpp`, `genotype_gvcf_tool.cpp`) keep
  their own local rules; in particular `genotype_gvcf_tool.cpp` implements the
  same htsjdk ordering by hand (`token_precision = (number <= 0.0 || number >= 1.0) ? 2 : 3`
  for `AF`/`MLEAF`), which coincides with `formatVCFDouble` for `d == 0`,
  `0.01 <= d < 1` and `d >= 1` but differs for `0 < d < 0.01`.  Not changed
  here (out of scope; Mutect2 paths are off-limits this round).
* **Not proven**: that no *other* native writer emits a zero `AF`/`MLEAF`
  elsewhere.  The four `format_allele_frequency` call sites are the only ones in
  `hc_call.cpp`; the gVCF sites were shown identical by probe.
* The **serial** backend was exercised on 3 of the 10 cases by this oracle plus
  the full regression suite (289/289); the remaining cases were run on OpenMP
  only.
* The **new oracle is not registered in CTest** (editing
  `fastgatk-native/CMakeLists.txt` is forbidden this round), so `ctest` alone
  does not protect this rule yet — same caveat as the previous round’s oracles.
  Registering it needs one `add_test` block; the test itself is self-contained
  and runs in ~7 minutes for all ten cases.
* **Target 2 (missing `GP`/`PG` FORMAT fields under
  `--genotype-assignment-method USE_POSTERIOR_PROBABILITIES`, the earlier
  round’s §D.2) was deliberately not started**: this round’s time-box says it
  may only begin once target 1 is fixed *and gated*, and gating target 1
  (159-test ctest selection + the three unselected fixture gates + the
  289-test double-backend regression) consumed the budget.  Nothing in this
  round changes or measures that divergence; it remains exactly as reported,
  and the new oracle’s `GT:AD:DP:GQ:PL` rows are unaffected by it because no
  case uses the posterior assignment method.
