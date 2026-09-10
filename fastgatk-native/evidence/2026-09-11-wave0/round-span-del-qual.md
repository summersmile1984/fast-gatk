# Round: HaplotypeCaller spanning-deletion `QUAL`/`QD` divergence (no `--alleles`)

Date: 2026-09-11. Repo: `/home/turing-agents/Documents/fast-gatk`, git `0c06575`
(clean working tree at start). **No `git commit`, no branch, no branch switch, no
edit to `fastgatk-native/CMakeLists.txt`, no edit to any root `*.md`, no Mutect2
change, no attempt at the separate symbolic-ALT emission gap.** Binaries:
`fastgatk-native/build/fastgatk-hc-call` (OpenMP/Kokkos),
`fastgatk-native/build-serial/fastgatk-hc-call` (Serial). Oracle: pinned GATK
4.6.2.0 (`third_party/jdk17/bin/java -Xmx1g -jar
third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar
HaplotypeCaller`). All scratch fixtures live under `tempfile.TemporaryDirectory`;
nothing was written outside `.diag/` (ignored) and the two files listed in §8.

**Bottom line: the observation is a GENUINE native-vs-GATK divergence, it
reproduces exactly (including the reported literals `282.04` / `291.07`), it is
NOT an artifact of the original probe's fixture, and it is localized to one
initializer. A minimal, source-grounded fix was landed and fully gated. The tree
DOES CONTAIN the change.**

The probe's own hypothesis is **partially refuted**: native does *not* exclude
`*` from the allele-frequency calculation. The real cause is the `*` allele's
Dirichlet **prior pseudocount**.

---

## 1. The fixture (reproducible, not anecdotal)

New reusable strict gate:
`fastgatk-native/scripts/verify_hc_span_del_qual_oracle.py` (modelled on
`verify_hc_alleles_overlap_gate_oracle.py`: `main() -> int`, `--native` /
`FASTGATK_HC_BINARY`, `TemporaryDirectory` scratch, `--case`, `--threads`,
`--expect-divergence`, final JSON status line, non-zero exit on mismatch).

Reference and reads are synthesised inside a `TemporaryDirectory`:

* 1500 bp synthetic `chr1`, `random.seed(11)`, 10×`A` homopolymer at 600, 8×`CAG`
  tandem repeat at 900. Base `chr1:698` is `C` and bases `698-700` are the `CCC`
  run the 1 bp deletion slides in.
* 30 reference-only reads: 300M from 440, step 3 (cover 440–826).
* 12 reads carrying a real `chr1:698 C>A` SNP: 300M from 440, step 3, base 698
  replaced by `A`.
* 6 reads carrying the 1 bp deletion: from 560, step 4, two encodings
  (`frameshift` = the encoding the reported probe reused — sequence removal with
  a 299M CIGAR; `cigar` = a real `138M1D161M` deleting base 698).
* interval `chr1:500-780`, `--min-pruning 1`,
  `--create-output-variant-index false`, `--add-output-vcf-command-line false`,
  **no `--alleles`**; native additionally `--threads 2`.

Exact command (paths are absolute inside the script; the SAM/FASTA/BAM are
written to a temp dir and both callers get identical `-R/-I/-L` arguments):

```
timeout 2400 python3 fastgatk-native/scripts/verify_hc_span_del_qual_oracle.py            # strict, exit 0 expected
timeout 2400 python3 fastgatk-native/scripts/verify_hc_span_del_qual_oracle.py --expect-divergence   # diagnostic, exit 0
```

Inside the script, per case:

```
java -Xmx1g -jar <pinned gatk 4.6.2.0> CreateSequenceDictionary -R ref.fa -O ref.dict --TRUNCATE_NAMES_AT_WHITESPACE true
java -Xmx1g -jar <pinned gatk 4.6.2.0> SortSam -I case.sam -O case.bam -SO coordinate --CREATE_INDEX true
java -Xmx1g -jar <pinned gatk 4.6.2.0> HaplotypeCaller -R ref.fa -I case.bam -L chr1:500-780 --min-pruning 1 \
     --create-output-variant-index false --add-output-vcf-command-line false -O gatk.vcf
fastgatk-native/build/fastgatk-hc-call  -R ref.fa -I case.bam -L chr1:500-780 --min-pruning 1 \
     --create-output-variant-index false --add-output-vcf-command-line false --threads 2 -O native.vcf
```

Comparison is field-by-field over every data row (REF/ALT/QUAL/FILTER, each INFO
sub-field, FORMAT keys, each FORMAT value); `#` provenance lines and headers are
ignored by construction. `QUAL`/`QD` are printed explicitly for every row.

---

## 2. Literal rows

### 2.1 Pre-fix (native binary built from `0c06575`), case `reported-frameshift-span-del`

Raw log: `.diag/span_del_prefix_final.log` (strict run, **exit 1**).

```
GATK    chr1  697  .  TC  T  62.60  .  AC=1;AF=0.500;AN=2;BaseQRankSum=0.000;DP=48;ExcessHet=0.0000;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=1.30;ReadPosRankSum=3.633;SOR=0.026   GT:AD:DP:GQ:PL  0/1:42,6:48:70:70,0,1410
GATK    chr1  698  .  C   A  282.04  .  AC=1;AF=0.500;AN=2;BaseQRankSum=0.000;DP=48;ExcessHet=0.0000;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=5.88;ReadPosRankSum=-2.995;SOR=0.162  GT:AD:DP:GQ:PL  0/1:36,12:48:99:400,0,1235
NATIVE  chr1  697  .  TC  T  62.60  .  AC=1;AF=0.500;AN=2;BaseQRankSum=0.000;DP=48;ExcessHet=0.0000;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=1.30;ReadPosRankSum=3.633;SOR=0.026   GT:AD:DP:GQ:PL  0/1:42,6:48:70:70,0,1410
NATIVE  chr1  698  .  C   A  291.07  .  AC=1;AF=0.500;AN=2;BaseQRankSum=0.000;DP=48;ExcessHet=0.0000;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=6.06;ReadPosRankSum=-2.995;SOR=0.162  GT:AD:DP:GQ:PL  0/1:36,12:48:99:400,0,1235
```

Field diff of the divergent row (GATK -> NATIVE), produced by the gate:

```
QUAL '282.04' -> '291.07'; INFO QD '5.88' -> '6.06'
```

That is the **entire** difference: ALT, FILTER, every other INFO sub-field
(`AC/AF/AN/BaseQRankSum/DP/ExcessHet/FS/MLEAC/MLEAF/MQ/MQRankSum/ReadPosRankSum/SOR`)
and every FORMAT field (`GT:AD:DP:GQ:PL = 0/1:36,12:48:99:400,0,1235`) are
identical, and the `chr1:697` row is byte-identical. These are exactly the
literals the prior probe reported (GATK 282.04 / QD 5.88, native 291.07 / QD
6.06), reproduced on the first run of the rebuilt fixture.

### 2.2 Post-fix (final binary)

Raw log: `.diag/span_del_postfix.log` (strict run, **exit 0**).

```
GATK    chr1  698  .  C   A  282.04  .  ... QD=5.88 ...  GT:AD:DP:GQ:PL  0/1:36,12:48:99:400,0,1235
NATIVE  chr1  698  .  C   A  282.04  .  ... QD=5.88 ...  GT:AD:DP:GQ:PL  0/1:36,12:48:99:400,0,1235
data_rows_byte_identical=True
```

---

## 3. Verdict: **genuine divergence** (not an artifact)

Attribution evidence (raw log `.diag/span_del_full_diag.log` for the pre-fix
sweep, `.diag/span_del_prefix_final.log` for the final pre-fix run):

| case | reads / flags | pre-fix GATK QUAL | pre-fix native QUAL | verdict |
|---|---|---|---|---|
| `reported-frameshift-span-del` | 30 ref + 12 SNP + 6 del (probe's encoding) | 282.04 | 291.07 | diverges |
| `clean-cigar-del-6` | same, del encoded as real `138M1D161M` | 282.04 | 291.07 | diverges |
| `reported-deep-ref-60` | 60 ref + 12 SNP + 6 del | 282.01 | 290.78 | diverges |
| `reported-narrow-window` | `-L chr1:600-760` | 282.04 | 291.07 | diverges |
| `reported-max-mnp-distance-0` | `--max-mnp-distance 0` | 282.04 | 291.07 | diverges |
| `reported-max-mnp-distance-2` | `--max-mnp-distance 2` | 282.04 | 291.07 | diverges |
| `no-span-del-control` | same pileup, **no deletion reads** | 403.64 | 403.64 | parity |
| `clean-cigar-del-18` / `frameshift-del-18` | 18 deletion reads | 454.60 | 454.60 | parity |
| `reported-ploidy-3` (out of scope) | `--sample-ploidy 3` | 358.64 (ALT `*,A`) | 358.64 (ALT `A`) | QUAL parity; ALT gap is the *other* defect |

Why this is not a fixture artifact:

* the divergence is **not** tied to the probe's odd read encoding — a real
  `138M1D161M` CIGAR deletion reproduces it identically;
* it is not tied to the deletion depth, window, ploidy-2 genotyping options or
  reference depth within the probed range;
* the control without deletion reads is byte-identical, so the trigger is the
  deletion (the `*`/spanning-deletion path), not the SNP;
* this is **not** the separately-tracked symbolic-ALT emission gap: at `chr1:698`
  both callers emit the same biallelic `C A` row and the `chr1:697 TC>T` deletion
  row is byte-identical on both sides — the only difference is `QUAL`/`QD`.

Refutation attempts that **failed** to refute it: original read encoding vs real
CIGAR (both diverge), 6 vs 18 deletion reads (6 diverges), 30 vs 60 reference
reads (both diverge), `-L chr1:500-780` vs `chr1:600-760` (both), `--max-mnp-distance`
0/1(default)/2 (all), ploidy 2 vs 3 (ploidy 3 has no QUAL divergence at all).
The only probe that produced a matching `chr1:698` row trivially was
`*del-18`: with 18 deletion reads neither caller emits a `chr1:698` row at all
(only `chr1:697`), so there is no confidence to disagree about — the divergence
is removed by removing the SNP row, not by removing the deletion's effect.

---

## 4. Localization (native) + GATK citation

### 4.1 The hypothesis under test is refuted in its specifics

`fastgatk::calling::Result run(...)`, `fastgatk-native/src/calling_pipeline.cpp`,
already has a dedicated spanning-deletion confidence policy
(`calling_pipeline.cpp:17297-17360`, post-fix numbering). Instrumented with a
temporary `FASTGATK_DEBUG_SPAN_DEL_AF` print (since removed, see §6) on the
reported fixture, native's OpenMP binary reported:

```
[FASTGATK_SPAN_DEL_AF] cand=0 pos=696 ref=TC alt=T derived=0 biallelic_qual=62.6018 spanning_pl=
[FASTGATK_SPAN_DEL_AF] cand=1 pos=697 ref=C alt=A derived=1 biallelic_qual=392.64 spanning_pl=400,0,1235,262,848,1204,
[FASTGATK_SPAN_DEL_AF] cand=1 pos=697 af_qual=291.074 log10_p_no_variant=-29.1074 effective_counts=1,1,7.80888e-30, samples=1
```

So native **does** build the three-allele `[REF, concrete ALT, *]` genotype
likelihood row (`derive_pairhmm_spanning_deletion_pl`, `calling_pipeline.cpp:5959-5989`)
and **does** feed it to the AF kernel
(`calculate_allele_frequency_kokkos(spanning_pl, 1, 3, 2, ...)`). The `*` allele
is *not* dropped before the confidence computation, so the probe's
"ordinary diploid path forces a concrete-only allele set" explanation is not the
cause of the QUAL gap. (Those concrete-only flags — `hc_call.cpp:3238`, `3248` —
govern the *emitted* allele list/AD/PL, not the confidence.)

### 4.2 The actual defect: the `*` allele's Dirichlet prior pseudocount

`fastgatk-native/src/calling_pipeline.cpp`, function
`fastgatk::calling::Result run(...)` (`calling_pipeline.cpp:13975`-ish start;
the defect at **`calling_pipeline.cpp:17311-17324`** post-fix / `17308-17312`
pre-fix), inside the `if (!options.somatic_mode && !result.spanning_deletion_read_likelihoods.empty())`
block. Pre-fix:

```cpp
const auto spanning_pseudocount = options.indel_heterozygosity * ref_pseudocount;
```

Matching GATK source: `gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/genotyper/afcalc/AlleleFrequencyCalculator.java`

```java
172:    private AFCalculationResult calculate(final int numAlleles,
173:                                          final List<Allele> alleles,
...
175:        final double[] priorPseudocounts = alleles.stream()
176:                .mapToDouble(a -> a.isReference() ? refPseudocount : (a.length() == refLength ? snpPseudocount : indelPseudocount)).toArray();
```

with `refLength` bound at the public entry point
(`AlleleFrequencyCalculator.java:139`: `calculate(numAlleles, alleles,
vc.getGenotypes(), defaultPloidy, vc.getReference().length())`).

The symbolic spanning-deletion allele is a **1 bp allele**, verified directly
against the pinned jar (not against prose):

```
$ third_party/jdk17/bin/javac -cp <pinned gatk 4.6.2.0 jar> -d .diag/jprobe .diag/jprobe/SpanLen.java
$ third_party/jdk17/bin/java -cp <pinned gatk 4.6.2.0 jar>:.diag/jprobe SpanLen
SPAN_DEL display=* length=1 isRef=false isSymbolic=false
NON_REF length=0
ref C length=1
```

So on a 1 bp REF record GATK gives `*` the **SNP** pseudocount
(`a.length() == refLength` -> 1 == 1), and only on longer REF records the indel
one. Native always used the indel prior. The prior feeds the Dirichlet EM
(`AlleleFrequencyCalculator.java:184-193`), which moves
`log10PNoVariant = sum over samples of log10SumLog10(posteriors of REF-only and
REF/* genotypes)` (`AlleleFrequencyCalculator.java:198-222`, the
`spanningDeletionPresent` branch using
`genotypeIndicesWithOnlyRefAndSpanDel`, lines 101-113), i.e. exactly the
`AFresult.log10ProbOnlyRefAlleleExists()` value that becomes `QUAL`
(`GenotypingEngine.java:158-165`; `usePosteriorProbabilitiesToCalculateQual` is
`false` by default — `GenotypeCalculationArgumentCollection.java:33` — so the AF
result is the final QUAL).

Same locus, same PL row, only the prior changed:

| | `*` prior | `log10_p_no_variant` | QUAL |
|---|---|---|---|
| pre-fix native | indel | -29.1074 | 291.074 -> `291.07` |
| post-fix native | SNP (REF length 1) | -28.2043 | 282.043 -> `282.04` |
| GATK 4.6.2.0 | SNP (`a.length()==refLength`) | (not observable) | `282.04` |

---

## 5. The change (minimal, one initializer)

`git diff` (also saved as `.diag/span_del_fix.diff`):

```diff
diff --git a/fastgatk-native/src/calling_pipeline.cpp b/fastgatk-native/src/calling_pipeline.cpp
index e276f09..823b4bb 100644
--- a/fastgatk-native/src/calling_pipeline.cpp
+++ b/fastgatk-native/src/calling_pipeline.cpp
@@ -17308,7 +17308,19 @@ Result run(const io::ReadBatch& reads,
                 const auto concrete_pseudocount =
                     (concrete_indel ? options.indel_heterozygosity : options.heterozygosity) *
                     ref_pseudocount;
-                const auto spanning_pseudocount = options.indel_heterozygosity *
+                // AlleleFrequencyCalculator derives every allele's Dirichlet
+                // pseudocount from its length: `a.length() == refLength ?
+                // snpPseudocount : indelPseudocount`
+                // (AlleleFrequencyCalculator.java:175-176, refLength =
+                // vc.getReference().length()).  The symbolic spanning deletion
+                // is a 1 bp allele (`Allele.SPAN_DEL.length() == 1`), so it
+                // takes the SNP prior on a 1 bp REF record and the indel prior
+                // otherwise.  Using the indel prior unconditionally changed
+                // P(no variant) for REF/* genotypes and therefore QUAL/QD.
+                const auto spanning_pseudocount =
+                    (candidate_reference(result.candidates[candidate]).size() == 1
+                         ? options.heterozygosity
+                         : options.indel_heterozygosity) *
                     ref_pseudocount;
                 if (!(concrete_pseudocount > 0.0) || !std::isfinite(concrete_pseudocount) ||
                     !(spanning_pseudocount > 0.0) || !std::isfinite(spanning_pseudocount))
```

Why this is minimal **and** safe: for any candidate whose REF allele is not 1 bp
the ternary is the old expression verbatim (`indel_heterozygosity`), so the
behaviour changes **only** for 1 bp REF records that have a spanning-deletion PL
row — precisely the records where GATK's own rule changes the prior. The
sibling concrete-ALT prior is already length-driven
(`concrete_indel = ref.size() != alt.size()`), which is the same rule, so this
change makes the two priors consistent with the single GATK expression. No
temporary instrumentation remains in the tree: the `FASTGATK_DEBUG_SPAN_DEL_AF`
prints quoted in §4.1 were added to `calling_pipeline.cpp` while localizing and
removed again before this diff was taken (the only surviving change is the one
shown above; §7 records the final `git status`/`git diff --stat`).

---

## 6. Gate ladder (in the mandated order), final source

| step | command | result |
|---|---|---|
| a | `timeout 1800 cmake --build fastgatk-native/build --target fastgatk-hc-call -j 16` then `timeout 2400 python3 fastgatk-native/scripts/verify_hc_span_del_qual_oracle.py` | **PASS, exit 0**, `"status": "pass"`, `violations: []`, 9/9 gated cases byte-identical (log `.diag/ladder_a_gate.log`) |
| b | `ctest --test-dir fastgatk-native/build -R 'window-invariance-gatk-oracle\|alleles-overlap-gate-oracle\|alleles-gatk' -V` | **PASS, exit 0**, `100% tests passed, 0 tests failed out of 5` (log `.diag/ladder_b_ctest.log`) |
| c | `python3 fastgatk-native/scripts/verify_hc_alleles_gatk_oracle.py` | **PASS, exit 0**, `"status": "pass"` (log `.diag/ladder_c_alleles.log`) |
| d | serial rebuild + `fastgatk-native/scripts/run_regression.sh --label 'span-del QUAL'` | **PASS, exit 0**, `284/284` on **both** backends (omp and serial; `100% tests passed, 0 tests failed out of 284` in each log) |

Pre-fix control for step (a) (the gate must fail without the fix): the fix was
temporarily reverted, the OpenMP target rebuilt, and the final gate re-run —
**exit 1**, 6 of 9 gated cases divergent (`reported-frameshift-span-del`,
`clean-cigar-del-6`, `reported-narrow-window`, `reported-max-mnp-distance-0`,
`reported-max-mnp-distance-2`, `reported-deep-ref-60`), raw log
`.diag/span_del_prefix_final.log`. The fix was then restored byte-for-byte from
`.diag/calling_pipeline.fixed.cpp` and the target rebuilt; step (a) above is the
post-restore run of that same final source.

### 6.1 Step (d)

```
$ third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/cmake --build fastgatk-native/build-serial --target fastgatk-hc-call -j 16
[100%] Built target fastgatk-hc-call
$ fastgatk-native/scripts/run_regression.sh --label 'span-del QUAL'
# 双后端回归证据
- 标签：span-del QUAL   git：`0c06575`（未提交变更 2 项）   过滤：`<全量>`   并行度：8
| 后端   | 构建目录                              | 结果 | 通过/总数 | 耗时       |
| omp    | OpenMP (fastgatk-native/build)        | 通过 | 284/284   | 1080.31 sec |
| serial | Serial (fastgatk-native/build-serial) | 通过 | 284/284   | 1015.24 sec |
$ cat .diag/regression/20260911-053816/{omp,serial}.status
PASS
PASS
$ grep -E "tests passed" .diag/regression/20260911-053816/{omp,serial}.log
100% tests passed, 0 tests failed out of 284     <- omp
100% tests passed, 0 tests failed out of 284     <- serial
```

Evidence directory: `.diag/regression/20260911-053816/`
(`evidence.md`, `summary.txt`, `omp.log`, `omp.status`, `serial.log`, `serial.status`);
script exit status **0**. No test regressed: the suite is still 284/284 on both
backends.

Extra cross-backend check (not required by the ladder, run for completeness):
the new gate against the **Serial** binary also passes —
`python3 fastgatk-native/scripts/verify_hc_span_del_qual_oracle.py --native
fastgatk-native/build-serial/fastgatk-hc-call` -> **exit 0**, `"status": "pass"`,
9/9 gated cases byte-identical, the only non-identical probe being the
out-of-scope ploidy-3 one (log `.diag/span_del_postfix_serial.log`).

Diagnostic mode also works on the final script:
`python3 fastgatk-native/scripts/verify_hc_span_del_qual_oracle.py
--expect-divergence --case reported-frameshift-span-del --case reported-ploidy-3`
-> **exit 0**, `"status": "diagnostic"`, `violations: []`, 1 out-of-scope note
(log `.diag/span_del_diag_final.log`).

---

## 7. Does the tree contain the change?

**YES — the tree CONTAINS one production change.**

```
$ git status --short
 M fastgatk-native/src/calling_pipeline.cpp
?? fastgatk-native/scripts/verify_hc_span_del_qual_oracle.py

$ git diff --stat
 fastgatk-native/src/calling_pipeline.cpp | 14 +++++++++++++-
 1 file changed, 13 insertions(+), 1 deletion(-)
```

* Change kept in the working tree (no commit, no branch).
* `fastgatk-native/CMakeLists.txt` and all root `*.md` untouched; the new gate is
  **not** registered in CTest (that file is off-limits this round) and must be
  run directly, or registered by the orchestrator.
* Mutect2 untouched. The separately-tracked symbolic-ALT emission gap was not
  attempted (the ploidy-3 probe that exposes it is reported, not gated — §4 of
  the script docstring).
* Scratch logs under `.diag/` (git-ignored): `span_del_prefix_final.log`,
  `span_del_postfix.log`, `span_del_full_diag.log`, `ladder_a_gate.log`,
  `ladder_b_ctest.log`, `ladder_c_alleles.log`, `span_del_fix.diff`.

---

## 8. What remains unproven / open

1. **The prior probe's exact fixture is still unknown** (it was lost). My
   reconstruction reproduces the reported GATK `282.04`/QD `5.88` and native
   `291.07`/QD `6.06` exactly, plus the reported "everything else identical"
   property including the `chr1:697` row, so it is almost certainly the same
   fixture — but that is inference from the literals, not a byte-level identity
   proof. (The probe recorded the `chr1:697` row only as "identical"; in my
   rebuild it is `chr1 697 . TC T 62.60 ... QD=1.30 ... 0/1:42,6:48:70:70,0,1410`
   on both sides.)
2. **`*` on a multi-base REF record is untested here.** My change preserves the
   indel prior for `REF.length() != 1`, mirroring
   `AlleleFrequencyCalculator.java:176`, but I did not build a fixture in which a
   spanning deletion coexists with a longer REF allele at the same locus, so the
   `indel` branch of the new ternary is source-derived, not gate-proven. (It is
   bit-identical to the previous behaviour, so it cannot regress that path.)
3. **The `*` allele's prior in the gVCF/`<NON_REF>` path was deliberately not
   touched.** `hc_call.cpp:682-687` selects a spanning-deletion prior from
   `genotype_indel_heterozygosity`; that is a different quantity for a different
   (genotype-prior, not AF-Dirichlet) calculation, and `NON_REF` has
   `length == 0`, so its GATK prior is always the indel one. Not investigated.
4. **The ploidy-3 ALT divergence is real and still open** (GATK `*,A` vs native
   `A` at `chr1:698`, QUAL/QD already equal at 358.64/7.47). It is exactly the
   out-of-scope symbolic-ALT emission gap from
   `round-forced-alleles-emission.md`; this round only measured it.
5. **Not proven**: that `282.043` and GATK's `282.04` are the same double rather
   than two doubles that round to the same 2-decimal string. The gate compares
   serialized `VCF` fields, so the acceptance criterion is the serialization, not
   the internal double.
6. **Speculation (labelled as such)**: that this same prior mismatch is the
   reason the `*,A` row in the `--alleles` fixture (`round-forced-alleles-emission.md`
   §2.3, where GATK's `*,A` row has `QD=0.00` and native never emits the row) also
   disagrees internally. Not measured this round — no `--alleles` fixture was run.
