# Round: `--alleles` forced-allele **emission bypass** (LowQual + symbolic `*` ALT)

Date: 2026-09-11. Repo: `/home/turing-agents/Documents/fast-gatk`, git `f4b46a7`
(clean working tree at start; see §1). **No `git commit`, no branch, no branch switch,
no edit to `fastgatk-native/CMakeLists.txt`, no edit to any root `*.md`.**
Binaries: `fastgatk-native/build/fastgatk-hc-call` (OpenMP/Kokkos),
`fastgatk-native/build-serial/fastgatk-hc-call` (Serial). Oracle: pinned GATK 4.6.2.0
(`third_party/jdk17/bin/java -Xmx1g -jar
third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar HaplotypeCaller`).
All scratch fixtures live under `tempfile.TemporaryDirectory`; nothing was written into
the source tree.

**Bottom line: STEP 1 and STEP 2 are complete and the divergence is reproduced and
localized to a single statement; STEP 3 was STOPPED deliberately because the change is
neither minimal nor safe (§4), and the task explicitly makes that a sanctioned outcome.
The tree DOES NOT CONTAIN any production change: `git diff` is empty, so nothing had to
be reverted.** The only tree addition is the mandated strict gate script.

---

## 1. Changes in the working tree

```
$ git status --short
?? fastgatk-native/scripts/verify_hc_forced_alleles_emission_gate_oracle.py

$ git diff --stat
(empty — no tracked file was modified)
```

**Production change: NONE.** `fastgatk-native/src/calling_pipeline.cpp`,
`fastgatk-native/src/hc_call.cpp`, `fastgatk-native/include/fastgatk/calling/pipeline.hpp`
and every other tracked file are byte-identical to `f4b46a7`. No temporary
instrumentation was ever added, so there is nothing to remove (no `FASTGATK_TEMP_*`
probe exists in the tree). Both HC binaries on disk were built from the same source
content that `f4b46a7` records (`calling_pipeline.cpp` mtime `04:23:59`, binaries
`04:24:44` (OpenMP) and `04:26:04` (Serial), commit `05:05:34`); neither was rebuilt this
round, so the 284-test result for `f4b46a7` applies unchanged.

The single new file is the STEP 1 strict gate
`fastgatk-native/scripts/verify_hc_forced_alleles_emission_gate_oracle.py` (untracked,
not registered in CTest — the orchestrator registers gates).

---

## 2. STEP 1 — strict gate (built first, before touching production code)

### 2.1 The re-established triggering fixture

The prior round's report (`fastgatk-native/evidence/2026-09-11-wave0/round-trackb-fix-report.md`
§5.2) only described the fixture in prose. Re-established exactly as follows, in a
`TemporaryDirectory`, by the new oracle:

* 1500 bp synthetic `chr1` (`random.seed(11)`, a 10×`A` homopolymer at 600, an 8×`CAG`
  tandem repeat at 900);
* 30× 300M reference-only reads from 440 (covers 440–739);
* **6 reads carrying an assembled 1 bp deletion at the reference base 700**, which
  left-aligns into the local `C` run and is emitted as `chr1:697 TC>T`;
* `--alleles` feature VCF (bgzipped + tabix-indexed, ascending order):

```
chr1  698  .  C   A   .  PASS  .
chr1  700  .  CA  C   .  PASS  .
```

* interval `chr1:500-780`, `--min-pruning 1`,
  `--create-output-variant-index false`, `--add-output-vcf-command-line false`,
  native additionally `--threads 2`.

The forced SNP `chr1:698 C>A` sits on the base that the emitted deletion `chr1:697 TC>T`
removes, i.e. inside a deletion that has real read support — this is the "forced allele
retained but implausible" trigger named in the task.

### 2.2 The gate

`fastgatk-native/scripts/verify_hc_forced_alleles_emission_gate_oracle.py` — modelled on
`verify_hc_alleles_overlap_gate_oracle.py` (`main() -> int`, `argparse` with
`--native` / `FASTGATK_HC_BINARY`, `TemporaryDirectory` scratch, final JSON status line,
`--expect-divergence` inversion, `--case`, `--threads`).

* **Cases:** `forced-snp-inside-supported-deletion` (both feature records) and
  `forced-deletion-only` (same pileup, the SNP record removed) — the second case is the
  trigger control: if the `*,…` row survives without the forced SNP, the fixture no
  longer isolates the trigger, and the gate reports a violation for that too.
* **Modes:** `alleles` and `drop-alleles` (the `--alleles` flag removed). The gate
  additionally requires that the `--drop-alleles` GATK row count be **strictly fewer**
  than the with-`--alleles` row count, so the divergent row is proven `--alleles`-driven
  rather than a pileup/discovery difference.
* **Assertion:** native **data rows must be byte-identical to GATK's**, field by field
  (REF/ALT/QUAL/FILTER/INFO sub-fields/FORMAT keys and values); header and `#` lines are
  ignored by construction. Non-zero exit on any mismatch. Field labels are correct here
  (the sibling oracle mislabels VCF column 6 as `QUAL`; column 6 is `FILTER`, and this
  script labels 3/4/5/6 = REF/ALT/QUAL/FILTER).

### 2.3 Pre-fix run — **it FAILS as required** (exit 1, `status: fail`)

Command:

```
timeout 900 python3 fastgatk-native/scripts/verify_hc_forced_alleles_emission_gate_oracle.py
```

Literal output (raw log `.diag/forced_alleles_gate_prefix.log`):

```
[forced-snp-inside-supported-deletion / alleles] gatk_rows=3 native_rows=2 data_rows_byte_identical=False divergent_positions=['698/700', 'missing-in-native@2:POS=700']
    feature VCF rows: chr1	698	.	C	A	.	PASS	.; chr1	700	.	CA	C	.	PASS	.
    GATK   chr1	697	.	TC	T	98.60	.	AC=1;AF=0.500;AN=2;BaseQRankSum=0.000;DP=36;ExcessHet=0.0000;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=2.74;ReadPosRankSum=3.417;SOR=0.050	GT:AD:DP:GQ:PL	0/1:30,6:36:99:106,0,961
    GATK   chr1	698	.	C	*,A	0	LowQual	AC=1,0;AF=0.500,0.00;AN=2;DP=36;ExcessHet=0.0000;FS=0.000;MLEAC=1,0;MLEAF=0.500,0.00;MQ=60.00;MQRankSum=0.000;QD=0.00;ReadPosRankSum=3.523;SOR=0.050	GT:AD:DP:GQ:PL	0/1:30,6,0:36:99:106,0,961,207,1055,1478
    GATK   chr1	700	.	CA	C	144.77	.	AC=0;AF=0.00;AN=2;DP=36;ExcessHet=0.0000;FS=0.000;MLEAC=0;MLEAF=0.00;MQ=60.00;SOR=0.001	GT:AD:DP:GQ:PL	0/0:36,0:36:99:0,108,1481
    NATIVE chr1	697	.	TC	T	98.60	.	AC=1;AF=0.500;AN=2;BaseQRankSum=0.000;DP=36;ExcessHet=0.0000;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=2.74;ReadPosRankSum=3.417;SOR=0.050	GT:AD:DP:GQ:PL	0/1:30,6:36:99:106,0,961
    NATIVE chr1	700	.	CA	C	144.77	.	AC=0;AF=0.00;AN=2;DP=36;ExcessHet=0.0000;FS=0.000;MLEAC=0;MLEAF=0.00;MQ=60.00;SOR=0.001	GT:AD:DP:GQ:PL	0/0:36,0:36:99:0,108,1481
    first divergent row field diff (GATK -> NATIVE): REF 'C' -> 'CA'; ALT '*,A' -> 'C'; QUAL '0' -> '144.77'; FILTER 'LowQual' -> '.'; INFO AC '1,0' -> '0'; INFO AF '0.500,0.00' -> '0.00'; INFO MLEAC '1,0' -> '0'; INFO MLEAF '0.500,0.00' -> '0.00'; INFO MQRankSum '0.000' -> None; INFO QD '0.00' -> None; INFO ReadPosRankSum '3.523' -> None; INFO SOR '0.050' -> '0.001'; FORMAT GT '0/1' -> '0/0'; FORMAT AD '30,6,0' -> '36,0'; FORMAT PL '106,0,961,207,1055,1478' -> '0,108,1481'
[forced-snp-inside-supported-deletion / drop-alleles] gatk_rows=1 native_rows=1 data_rows_byte_identical=True divergent_positions=[]
[forced-deletion-only / alleles] gatk_rows=2 native_rows=2 data_rows_byte_identical=True divergent_positions=[]
[forced-deletion-only / drop-alleles] gatk_rows=1 native_rows=1 data_rows_byte_identical=True divergent_positions=[]
# 1 violation(s):
#   - forced-snp-inside-supported-deletion [alleles]: data rows differ at ['698/700', 'missing-in-native@2:POS=700'] (gatk_rows=3 native_rows=2)
{"status": "fail", ... "strict_mode": true, "release": "GATK 4.6.2.0", ...}
EXIT=1
```

**Key numbers.** GATK `3` data rows vs native `2`; the missing row is the forced SNP's
own locus `chr1:698` with `ALT=*,A`, `QUAL=0`, `FILTER=LowQual`. All three controls pass
byte-identically, so the fixture's attribution is proven:

| case | mode | GATK rows | native rows | byte-identical |
|---|---|---|---|---|
| forced-snp-inside-supported-deletion | `--alleles` | 3 | 2 | **False** (1 violation) |
| forced-snp-inside-supported-deletion | `--drop-alleles` | 1 | 1 | True |
| forced-deletion-only | `--alleles` | 2 | 2 | True |
| forced-deletion-only | `--drop-alleles` | 1 | 1 | True |

`--drop-alleles` yields exactly one identical row on both sides, so the `698` row exists
**only** because of `--alleles`; removing only the SNP record makes the `*,…` row vanish
on both sides, which pins the forced SNP record as the trigger.

### 2.4 Diagnostic mode works

```
$ timeout 900 python3 fastgatk-native/scripts/verify_hc_forced_alleles_emission_gate_oracle.py --expect-divergence
... "status": "diagnostic" ...
DIAG EXIT=0
```

Raw log: `.diag/forced_alleles_gate_diag.log`. `python3 -m py_compile` on the script: OK.

---

## 3. Confirming the two reported GATK semantics against the pinned source

All quotes are from `gatk-source/` at this checkout.

**(a) `passesEmitThreshold` is bypassed when the forced-allele collection is
non-empty — CONFIRMED.**

`gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/genotyper/GenotypingEngine.java:153-170`:

```java
153:        final AFCalculationResult AFresult = alleleFrequencyCalculator.calculate(reducedVC, defaultPloidy);
154:        final Set<Allele> forcedAlleles = AssemblyBasedCallerUtils.allelesConsistentWithGivenAlleles(givenAlleles, vc);
155:        final OutputAlleleSubset outputAlternativeAlleles = calculateOutputAlleleSubset(AFresult, vc, forcedAlleles);
...
163:        final double phredScaledConfidence = (-10.0 * log10Confidence) + 0.0;
...
167:        if ( !passesEmitThreshold(phredScaledConfidence, outputAlternativeAlleles.siteIsMonomorphic) && !emitAllActiveSites()
168:                && noAllelesOrFirstAlleleIsNotNonRef(outputAlternativeAlleles.alleles) && forcedAlleles.isEmpty()) {
169:            return null;
170:        }
```

`forcedAlleles.isEmpty()` at line 168 is the bypass; `passesEmitThreshold`
(`GenotypingEngine.java:425-432`) is `(outputMode == EMIT_ALL_CONFIDENT_SITES || !bestGuessIsRef) && passesCallThreshold(conf)`,
and `passesCallThreshold(conf)` is `conf >= configuration.genotypeArgs.standardConfidenceForCalling`
(line 430-432). Because `passesEmitThreshold(...)` is `false` for our row
(`phredScaledConfidence == 0 < 30`) but `forcedAlleles` is non-empty, the null return is
skipped.

**Instead the record is merely `LowQual`-filtered — CONFIRMED.**
`GenotypingEngine.java:183-186`:

```java
183:        builder.log10PError(log10Confidence);
184:        if ( ! passesCallThreshold(phredScaledConfidence) ) {
185:            builder.filter(GATKVCFConstants.LOW_QUAL_FILTER_NAME);
186:        }
```

With `log10Confidence == 0` this serializes as `QUAL=0` + `FILTER=LowQual`, exactly the
observed row. (`forcedAlleles.isEmpty()` is the *only* forced-allele exemption in the
whole file — `grep -n forcedAlleles GenotypingEngine.java` returns lines 154, 155, 168,
296, 316.)

**(b) GATK emits the symbolic `*` ALT in that situation — CONFIRMED.**
`GenotypingEngine.java:304-327` (`calculateOutputAlleleSubset`):

```java
311:                final boolean isPlausible = afCalculationResult.passesThreshold(allele, configuration.genotypeArgs.standardConfidenceForCalling);
313:                //it's possible that the upstream deletion that spanned this site was not emitted, mooting the symbolic spanning deletion allele
314:                final boolean isSpuriousSpanningDeletion = GATKVCFConstants.isSpanningDeletion(allele) && !isVcCoveredByDeletion(vc);
316:                final boolean toOutput = (isPlausible || forceKeepAllele(allele) || isNonRefWhichIsLoneAltAllele || forcedAlleles.contains(allele) ) && !isSpuriousSpanningDeletion;
318:                siteIsMonomorphic &= !(isPlausible && !isSpuriousSpanningDeletion);
```

At `chr1:698` the merged VC's alleles are `[C, *, A]`: `*` is present because an
upstream deletion event was rewritten to it —
`HaplotypeCallerGenotypingEngine.java:380-397` `replaceWithSpanDels`/`replaceWithSpanDelVC`
(`.alleles(Arrays.asList(refAllele, Allele.SPAN_DEL))`, line 392) — and
`AssemblyBasedCallerUtils.createAlleleMapper` (`AssemblyBasedCallerUtils.java:696-707`)
maps the deletion-carrying haplotypes onto `Allele.SPAN_DEL`. `*` is *retained* here
because it **is** plausible (6 spanning reads) and **not** spurious:
`isVcCoveredByDeletion` (`GenotypingEngine.java:365-371`) is true, since the deletion
`chr1:697 TC>T` was emitted first and recorded by `recordDeletions`
(`GenotypingEngine.java:342-357`, `deletionSize = vc.getReference().length() - allele.length()`
→ interval `697..698`). `A` is retained only through `forcedAlleles.contains(allele)`
(line 316), and `forcedAlleles` comes from
`AssemblyBasedCallerUtils.allelesConsistentWithGivenAlleles` (`AssemblyBasedCallerUtils.java:1005-1009`).
Hence the emitted ALT list is `*,A` — the symbolic allele is emitted *because* a forced
allele is present at a site whose only other ALT would have been a lone `*`, which is
separately suppressed by `GenotypingEngine.java:172-175`. The `*,A` order is the merged
VC's allele order.

**Native's corresponding numbers are already right where it counts.** The forced
candidate's QUAL computed by native is **0.0** — identical to GATK's — so the divergence
is *purely* an emission/allele-set gap, not a confidence gap.

---

## 4. STEP 2 — localization (with evidence)

### 4.1 The native drop site

Instrumentation: pre-existing debug gates only, no code change
(`FASTGATK_DEBUG_CALL_GATES=1`, `FASTGATK_DEBUG_EVENTMAP=1`,
`FASTGATK_DEBUG_EVENTMAP_LIFECYCLE=1`, `FASTGATK_DEBUG_ANNOTATION_POSITION=<0-based pos>`).
Fixture built inside one command; nothing written outside the `TemporaryDirectory`.

Trace for the divergence fixture (OpenMP build; 0-based positions in the trace,
1-based in VCF):

```
[FASTGATK_EVENTMAP_RAW_SNP] tid=0 pos=699 ref=C alt=A support=6
[FASTGATK_EVENTMAP_RAW_SNP] tid=0 pos=700 ref=A alt=T support=6
[FASTGATK_EVENTMAP_MATERIALIZED] path=2 position=697 ref=C alt=A
[FASTGATK_EVENTMAP_MATERIALIZED] path=3 position=699 ref=CA alt=C
[FASTGATK_EVENTMAP_MATERIALIZED] path=4 position=698 ref=CCA alt=C
[FASTGATK_EVENTMAP_INJECT] path=1 tid=0 position=696 ref=TC alt=T support=36 appended=1
[FASTGATK_EVENTMAP_INJECT] path=2 tid=0 position=697 ref=C alt=A support=36 appended=1
[FASTGATK_EVENTMAP_INJECT] path=3 tid=0 position=699 ref=CA alt=C support=36 appended=1
[FASTGATK_EVENTMAP_INJECT] path=4 tid=0 position=698 ref=CCA alt=C support=36 appended=1
[FASTGATK_EVENTMAP_REPLAY] tid=0 position=697 ref=C alt=A inserted=0 reason=already-present
[FASTGATK_EVENTMAP_REPLAY] tid=0 position=699 ref=CA alt=C inserted=0 reason=already-present
[FASTGATK_EVENTMAP_REGION] tid=0 position=696 ref=TC alt=T graph_derived=1 raw_class=1 calling_class=1
[FASTGATK_EVENTMAP_REGION] tid=0 position=697 ref=C alt=A graph_derived=1 raw_class=1 calling_class=1
[FASTGATK_EVENTMAP_REGION] tid=0 position=698 ref=CCA alt=C graph_derived=1 raw_class=1 calling_class=1
[FASTGATK_EVENTMAP_REGION] tid=0 position=699 ref=CA alt=C graph_derived=1 raw_class=1 calling_class=1
[FASTGATK_CALL_GATE] pos=696 ref=TC alt=T support=6  somatic_pair_evidence=0 genotype=1 qual=98.6018 decision=emit
[FASTGATK_EVENTMAP_CALL] tid=0 position=696 ref=TC alt=T graph_derived=1 qual=98.6018
[FASTGATK_CALL_GATE] pos=697 ref=C alt=A support=0 somatic_pair_evidence=0 genotype=0 qual=0 threshold=30 decision=confidence_suppressed
[FASTGATK_CALL_GATE] pos=698 ref=CCA alt=C support=0 qual=1.44649e-14 decision=hom_ref
[FASTGATK_CALL_GATE] pos=699 ref=CA alt=C support=0 genotype=0 qual=144.775 decision=emit
[FASTGATK_EVENTMAP_CALL] tid=0 position=699 ref=CA alt=C graph_derived=1 qual=144.775
```

**Function and line.** The record is lost in
`fastgatk::calling::Result run(...)` (`fastgatk-native/src/calling_pipeline.cpp:13975`,
the per-region call loop at the end of the function) at

```
fastgatk-native/src/calling_pipeline.cpp:17521-17534
17521:        if (qual + 1.0e-10 < options.standard_confidence_for_calling) {
...
17531:                          << " decision=confidence_suppressed\n";
17532:            }
17533:            continue;            // <- the forced-but-implausible row is dropped here
17534:        }
```

together with the candidate's pre-computed confidence at

```
fastgatk-native/src/calling_pipeline.cpp:17289-17296   // forced-allele QUAL policy (gives 0.0, = GATK's QUAL)
fastgatk-native/src/calling_pipeline.cpp:17297-17361   // spanning-deletion AF policy (REF, concrete ALT, `*`)
```

The gate at 17521 has **no forced-allele exemption**: the only `forced_by_alleles_feature`
uses in the loop are the *negative* filters at 17437 (`low_support_singleton`) and 17462
(`genotype_host(i) == 0` early-out) plus the gVCF branch at 17485-17510. Native's
behaviour therefore corresponds to GATK's `return null` branch at
`GenotypingEngine.java:167-169` with the `&& forcedAlleles.isEmpty()` clause missing.

**Matching GATK citation.** `GenotypingEngine.java:167-170` (the bypass) and
`GenotypingEngine.java:184-186` (LowQual instead of dropping); the `*,A` ALT set comes
from `GenotypingEngine.java:304-327` with `isVcCoveredByDeletion`
(`GenotypingEngine.java:365-371`) and `recordDeletions` (`GenotypingEngine.java:342-357`).

### 4.2 Why the missing row is not merely "one `continue`"

The `chi` trace above shows native never had a `*` allele for this locus at all: native's
candidate set at 1-based 698 is the *biallelic* forced SNP `ref=C alt=A` (there is also a
separate graph-derived `pos=698 ref=CCA alt=C` virtual deletion, which is folded back as
`hom_ref`). GATK's row needs five things native does not produce at that locus:

| # | GATK field | Native today | Where it would have to come from |
|---|---|---|---|
| 1 | emitted at all, `QUAL=0` | dropped at `calling_pipeline.cpp:17521` | forced-allele exemption on the gate (QUAL is already `0.0`, i.e. already correct) |
| 2 | `FILTER=LowQual` | **no FILTER emission path exists**: every writer site hardcodes `.` (`hc_call.cpp:3615`, `3624`, `3862`, `4587`, `4598`, `4917`, `4945`); `LowQual` appears only in the header (`hc_call.cpp:3037`, `3746`); `GenotypeCall` (`pipeline.hpp:648-656`) has no filter field | new field + all 7 emission sites |
| 3 | `ALT=*,A` (`*` **before** the concrete ALT) | ordinary-VCF diploid output has no `*`: `hc_call.cpp:3238` (`if (include_spanning_deletion) return false;`, guarded by the ploidy-2 comment at 3234-3236) and `hc_call.cpp:3248` (`has_spanning_deletion = sample_ploidy != 2 && …`); the only `*` emission in the whole writer is the **gVCF** path, `hc_call.cpp:4922` | enable + re-order for the ordinary diploid writer |
| 4 | `GT:AD:DP:GQ:PL 0/1:30,6,0:36:99:106,0,961,207,1055,1478` | only biallelic PL/AD available in that path | Number=G and Number=R re-permutation (native's helpers order alleles `[REF, concrete ALT, *]`, GATK needs `[REF, *, concrete ALT]` → new index map `[0,3,5,1,4,2]`) |
| 5 | `QD=0.00;ReadPosRankSum=3.523;SOR=0.050;FS=0.000;MQ=60.00;MQRankSum=0.000;BaseQRankSum=0.000;DP=36;ExcessHet=0.0000;AC=1,0;MLEAC=1,0;AF=0.500,0.00;MLEAF=0.500,0.00` | **the annotations for this candidate are never computed at all** — measured: `FASTGATK_DEBUG_ANNOTATION_POSITION=697` on this fixture prints **no** `[FASTGATK_ANNOTATION_SUMMARY]` line (only `696` and `699` do), because `calculate_variant_annotations` is invoked inside the emit branch that `continue` at 17533 skips | annotations are computed per single concrete candidate and reused verbatim by the writer (`calls[best]->annotations`, `hc_call.cpp:3207`; the later published-allele variant is `output_calls[best]->annotations`, `hc_call.cpp:3494`); a multi-allelic record whose alt-supporting partition is the `*` allele (`ReadPosRankSum` 3.523 differs from the sibling row's 3.417 precisely because the ALT set differs) needs a new annotation pass over the published allele partition |

Item 5 is the one that makes the change unsafe to attempt blind: the annotation path is
owned by the per-candidate Host boundary (`calculate_variant_annotations`, called at
`calling_pipeline.cpp:17550`), and the writer never recomputes rank sums.

**Conclusion of STEP 2: the correct behaviour is not ambiguous — it is fully pinned by
the source — but the change is neither minimal nor safe inside this round.** It would
touch `calling_pipeline.cpp` (call gate + annotation boundary),
`include/fastgatk/calling/pipeline.hpp` (`GenotypeCall` gains a filter field) and the
shared diploid writer in `hc_call.cpp` (allele list, FILTER column, Number=G/Number=R
re-permutation, AF priors for the `*` allele, annotation recomputation) — code paths used
by **all 284 registered tests**, for one record shape. Per the task's STEP 2 clause I
therefore STOP here and report, rather than land a half-working change.

### 4.3 Recommendation (for the orchestrator)

Open a dedicated track "ordinary-VCF symbolic spanning-deletion ALT + forced-allele
emission bypass", scoped as follows:

1. **Emission bypass (small, safe on its own).** In the call-gate loop
   (`calling_pipeline.cpp:17521`), do not `continue` for a
   `forced_by_alleles_feature` candidate; instead mark the call `LowQual` and keep it,
   mirroring `GenotypingEngine.java:167-170` + `184-186`. *On its own this is not
   GATK-parity* and must not be landed as a "fix" — it would emit a biallelic `C A` row
   where GATK emits `C *,A`, i.e. it would replace a missing row with a wrong row.
2. **FILTER column.** Add a filter field to `GenotypeCall` and thread it through the
   seven emission sites listed above; keep everything else `.`.
3. **Symbolic `*` in the ordinary diploid writer.** The disabled branch at
   `hc_call.cpp:3238` / `3248` is the intended home; the gVCF path
   (`hc_call.cpp:4693-4922`) is the working reference implementation, including the
   "covered by an upstream emitted deletion" test that mirrors GATK's
   `isVcCoveredByDeletion` + `recordDeletions` in traversal order.
4. **Number=G/Number=R re-ordering** so `*` is emitted before the concrete forced ALT.
5. **Annotations** must be computed against the *published* allele partition, not reused
   from the single concrete candidate; `ReadPosRankSum=3.523` is the sentinel field that
   will fail if this is skipped.
6. Register `verify_hc_forced_alleles_emission_gate_oracle.py` in CTest at that point
   (it is deliberately left unregistered here).

---

## 5. STEP 4 status (gate ladder) on the unchanged tree

Because no production change was made, the ladder below is the *baseline* control; it is
reported for completeness, not as post-fix evidence.

| step | command | result |
|---|---|---|
| a | rebuild OpenMP HC then new strict gate | **not applicable** — no source change; new gate fails on the unchanged binary (exit 1, §2.3) |
| b | `ctest --test-dir fastgatk-native/build -R 'window-invariance-gatk-oracle\|alleles-overlap-gate-oracle' -V` | **PASS, exit 0**, `100% tests passed, 0 tests failed out of 3` — `1/3 Test #83: fastgatk-hc-window-invariance-gatk-oracle Passed 60.22 sec`, `2/3 Test #84: fastgatk-hc-ploidy-window-invariance-gatk-oracle Passed 29.44 sec`, `3/3 Test #85: fastgatk-hc-alleles-overlap-gate-oracle Passed 37.35 sec`, `Total Test time (real) = 127.01 sec`. Raw log `.diag/ctest_strict_gates_nofix.log` |
| c | `python3 fastgatk-native/scripts/verify_hc_alleles_gatk_oracle.py` | **PASS, exit 0**, `"status": "pass"` (raw log `.diag/alleles_oracle_baseline.log`) |
| d | `run_regression.sh` (284/284 both backends) | **not run** — it would re-test byte-identical binaries; `git diff` is empty and neither HC binary was rebuilt, so the 284-test result for `f4b46a7` carries over unchanged |

No revert was required: `git diff` is empty and both HC binaries were never rebuilt from
modified sources (§1).

---

## 6. Does the tree contain the change?

**NO — the tree DOES NOT CONTAIN any production change.** `git status --short` shows only
`?? fastgatk-native/scripts/verify_hc_forced_alleles_emission_gate_oracle.py` and
`git diff --stat` is empty. Nothing was added to, removed from or reverted in
`fastgatk-native/src` or `fastgatk-native/include`;
`fastgatk-native/CMakeLists.txt` and all root `*.md` files are untouched. No temporary
instrumentation exists (none was needed: the pre-existing `FASTGATK_DEBUG_*` gates were
sufficient).

---

## 7. What remains unproven / open

1. **The divergence itself is proven**, including its attribution to `--alleles`
   (`§2.3` table) and its native drop site (`§4.1`). What is **not** proven is whether a
   fix can be made byte-exact; §4.2 enumerates the five required deltas and item 5
   (annotation recomputation for the published allele partition) is the one I did not
   validate numerically. If a future round lands this, the expected per-field values are
   recorded literally in §2.3 and are the acceptance target.
2. **The `*,A` allele order** is asserted from GATK's observed output plus the merged-VC
   allele order (`AssemblyBasedCallerUtils.java:310-318` → `simpleMerge`,
   `HaplotypeCallerGenotypingEngine.java:174-184`). I did not trace
   `GATKVariantContextUtils.simpleMerge` to prove *why* `*` precedes `A`. **Speculation**:
   it follows from `replaceWithSpanDelVC` producing `[refBase, *]` for the event that
   starts before the locus and `simpleMerge` preserving first-seen allele order. Treat the
   ordering as an empirical requirement, not a derived one.
3. **Side observation, out of scope, measured (not fixed).** A pileup with 30 reference
   reads, 12 reads carrying a real `chr1:698 C>A` SNP and the same 6 deletion reads,
   **without** `--alleles`, gives at `chr1:698`:
   `GATK chr1 698 . C A 282.04 ... 0/1:36,12:48:99:400,0,1235` vs
   `NATIVE chr1 698 . C A 291.07 ... 0/1:36,12:48:99:400,0,1235` — identical ALT/INFO
   except `QUAL`/`QD` (282.04 / 5.88 vs 291.07 / 6.06), with the `chr1:697` row
   byte-identical. Both sides emitted a **biallelic** row (no `*`) here, so this is *not*
   the symbolic-ALT defect. **Speculation**: it is the same root area — GATK computes that
   record's AF over `REF, concrete ALT, *` and then subsets `*` out, while native's
   ordinary diploid path is forced concrete-only (`hc_call.cpp:3238`, `3248`). This is a
   *separate, non-`--alleles`* divergence and is **not** covered by the new gate; it would
   need its own fixture and its own strict gate before anyone fixes it.
4. **Still untested** (unchanged from the prior round): base-haplotype score-ranking ties
   (`NUM_HAPLOTYPES_TO_INJECT_FORCE_CALLING_ALLELES_INTO`), `--max-genotype-count` × NaN
   injected-haplotype scores, multi-sample/joint genotyping, symbolic feature ALTs, and
   the out-of-scope D3 exit-code divergence (`Null alleles are not supported`).
5. **The new oracle is not registered in CTest** (CMakeLists.txt is off-limits this
   round), so it is not part of the 284-test count; it must be run directly, or registered
   by the orchestrator when the fix lands. It is currently expected to **fail**; running
   it in a green suite requires `--expect-divergence`.

---

## 8. Artifacts

* `fastgatk-native/scripts/verify_hc_forced_alleles_emission_gate_oracle.py` — the STEP 1
  strict gate (new, untracked).
* `.diag/forced_alleles_gate_prefix.log` — pre-fix strict run, exit 1.
* `.diag/forced_alleles_gate_diag.log` — `--expect-divergence` run, exit 0.
* `.diag/ctest_strict_gates_nofix.log` — the three existing strict gates, exit 0.
* `.diag/alleles_oracle_baseline.log` — existing `--alleles` oracle, exit 0.
