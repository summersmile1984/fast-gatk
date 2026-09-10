# Round: "prior selected from allele length/properties" sweep (gVCF `*` / `<NON_REF>`)

Date: 2026-09-11 (local) / 2026-09-10T22:45Z. Repo `/home/turing-agents/Documents/fast-gatk`,
git `4d66251` (clean tree at start; the previous round's spanning-deletion fix is HEAD).
**No `git commit`, no branch, no switch, no edit to `fastgatk-native/CMakeLists.txt`, no edit to any
root `*.md`, no Mutect2 change, no attempt at the separately tracked symbolic-ALT emission gap.**
All fixtures are built inside `tempfile.TemporaryDirectory`; scratch logs live under `.diag/`.

Binaries: `fastgatk-native/build/fastgatk-hc-call` (OpenMP/Kokkos) and
`fastgatk-native/build-serial/fastgatk-hc-call` (Serial). Oracle: pinned GATK 4.6.2.0
(`third_party/jdk17/bin/java -Xmx1g -jar third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar HaplotypeCaller`).

**Bottom line: a second instance of the same root class was reproduced, localized and fixed — the
gVCF reference-confidence record's `*` pseudocount (`hc_call.cpp`, the gVCF candidate-site AF block),
which made every reference-confidence record carrying the spanning deletion print the wrong QUAL
(GATK 282.04 / native 291.07, and independently GATK 144.04 / native 153.07). Six further HC sites
(plus one in the out-of-scope GenotypeGVCFs tool) of the same class were localized by source citation
but could NOT be made observable; they are reported with their exact one-line fixes and a
recommendation. The tree CONTAINS one change
(`fastgatk-native/src/hc_call.cpp`, +15/-1) plus the new gate script.**

---

## 0. The rule being swept

`AlleleFrequencyCalculator` (AF/Dirichlet path) selects every allele's pseudocount from its **length**:

```java
// gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/genotyper/afcalc/AlleleFrequencyCalculator.java:175-176
final double[] priorPseudocounts = alleles.stream()
        .mapToDouble(a -> a.isReference() ? refPseudocount
                        : (a.length() == refLength ? snpPseudocount : indelPseudocount))
        .toArray();
```
`refLength = vc.getReference().length()` (`AlleleFrequencyCalculator.java:165`).

`GenotypePriorCalculator` (genotype-prior/GT-GQ path) uses **symbolic-ness then length**, and adds a
SNP normalization constant:

```java
// gatk-source/.../utils/genotyper/GenotypePriorCalculator.java:141-153
} else if (allele.isCalled() && !allele.isSymbolic()) {
    return allele.length() == referenceLength ? AlleleType.SNP : AlleleType.INDEL;
} else if (allele.equals(Allele.SV_SIMPLE_INS) || allele.equals(Allele.SV_SIMPLE_DEL)) {
    throw new IllegalArgumentException("cannot handle symbolic indels: " + allele);
} else { return AlleleType.OTHER; }
```
with `hetValues[SNP] = snpHet - log10(3)`, `hetValues[INDEL] = indelHet`,
`hetValues[OTHER] = max(snpHet, indelHet)` (`GenotypePriorCalculator.java:62-72, 118-121`).

**Measured htsjdk semantics** (probe program compiled against the pinned GATK jar, htsjdk 4.2.0;
`/tmp` scratch, output below). This is the fact the whole sweep hinges on: symbolic alleles are
**length 0**, `*` is a *called, non-symbolic* 1 bp allele:

```
SPAN_DEL   bases=*            length=1 isRef=false isCalled=true  isSymbolic=false
NON_REF    bases=<NON_REF>    length=0 isRef=false isCalled=true  isSymbolic=true
SV_DEL     bases=<DEL>        length=0 isRef=false isCalled=true  isSymbolic=true
A          bases=A            length=1 isRef=true  isCalled=true  isSymbolic=false
```

So on a 1 bp REF record: `*` → **SNP** pseudocount/type in both GATK paths, `<NON_REF>`/`<DEL>` →
**indel** pseudocount (AF path) / **OTHER** = `max(snpHet, indelHet)` (genotype-prior path).

---

## A. Audit table (deliverable A)

Every native site that turns an allele's length / indel-ness / symbolicity into a numeric prior,
heterozygosity, pseudocount or branch, checked against its GATK counterpart. `snp`/`indel` below
mean `options.heterozygosity` (1e-3) and `options.indel_heterozygosity` (1/8000).

| # | Native site | Native rule | GATK rule (citation) | Agree? | Reproduced? |
|---|---|---|---|---|---|
| 1 | `calling_pipeline.cpp:17320-17324` `spanning_pseudocount` (biallelic candidate span-del AF, ordinary VCF) | REF len `== 1` → snp, else indel | `AlleleFrequencyCalculator.java:175-176`; `Allele.SPAN_DEL.length()==1` | **agree** (fixed in HEAD `4d66251`) | yes, previous round |
| 2 | `calling_pipeline.cpp:17261-17277` `confidence_priors` per concrete candidate | `ref.size() != alt.size()` → indel | same rule; for *concrete* alleles `size()` inequality ⇔ `length() == refLength` | agree | n/a |
| 3 | **`hc_call.cpp:4820-4836`** gVCF candidate-site AF priors, `include_spanning_deletion` branch | `indel_heterozygosity` **unconditionally** for `*` | `*` is a 1 bp allele → **SNP** prior on a 1 bp REF | **DISAGREE** | **YES — fixed here** |
| 4 | `hc_call.cpp:4838-4841` gVCF writer `<NON_REF>` prior | `indel_heterozygosity` | htsjdk `NON_REF.length()==0 != refLength` → indel prior | agree | ✓ (probe: no QUAL mismatch on `<NON_REF>`-only rows) |
| 5 | `hc_call.cpp:3358-3360` arbitrary-ploidy **ordinary-VCF** writer, `has_spanning_deletion` | `indel_heterozygosity` unconditionally for `*` | same as #3 | **DISAGREE** (left unfixed, see §C.2) | not observable: triploid probe QUAL `358.64` on both sides before **and** after applying the fix |
| 6 | `hc_call.cpp:4485-4488` non-diploid / `hidden_spanning_deletion` **gVCF** writer | `indel_heterozygosity` unconditionally for `*` | same as #3 | **DISAGREE** (left unfixed, §C.2) | not reproduced: no haploid/triploid `-ERC` fixture emits a `*` record on either side |
| 7 | `hc_call.cpp:4489-4491` same writer, `<NON_REF>` | `indel_heterozygosity` | indel prior | agree | ✓ |
| 8 | `hc_call.cpp:682-694` `genotype_priors_for_group(..., include_spanning_deletion=true)` | `log10(indelHet)` for `*` (no `log10(3)`) | `GenotypePriorCalculator.java:141-153`: `*` is called & non-symbolic → **SNP type** → `log10(snpHet) - log10(3)` | **DISAGREE** (left unfixed, §C.2) | not observable: `GT`/`GQ` byte-identical in all six gated fixtures |
| 9 | `calling_pipeline.cpp:405-424` `prior_allele_type` | `value == "*"` counts as **symbolic** → `OTHER` → `max(snp,indel)` **without** `log10(3)` | `*` → **SNP** type **with** `log10(3)` | **DISAGREE** (left unfixed, §C.2) | not observable (same GT/GQ evidence as #8) |
| 10 | `calling_pipeline.cpp:6059-6074` `gvcf_candidate_is_non_monomorphic` | 3rd allele (`<NON_REF>`, see `derive_pairhmm_nonref_pl` `calling_pipeline.cpp:6002-6028`) gets `other_heterozygosity = max(snp, indel)` = **snp** | AF calc: `NON_REF.length()==0` → **indel** | **DISAGREE** (left unfixed, §C.2) | not reproduced (only shifts the site-level retention decision; all gated fixtures agree) |
| 11 | `hc_call.cpp:2915-2926` `estimate_mle_allele_counts(..., include_non_ref=true)` | `prior.back() = max(snp, indel)` for `<NON_REF>` | indel prior | **DISAGREE** (left unfixed, §C.2) | not reproduced (fallback path only) |
| 12 | `genotype_gvcf_tool.cpp:2047-2062` `build_cohort_prior_pseudocounts` (GenotypeGVCFs, out of HC scope) | `alleles[i].size() == reference_length ? snp : indel` | AF calc length rule | agree (string size vs htsjdk `length()` happen to fall on the same side for `*` and `<NON_REF>`) | n/a |
| 13 | `genotype_gvcf_tool.cpp:1585-1615` `build_log10_genotype_priors` (GenotypeGVCFs) | `allele == "*" ‖ front=='<'` → `Other` | `*` → SNP type | **DISAGREE** (out of HC scope, not touched) | n/a |
| 14 | kernels: `fastgatk-kernels/src/genotype.cpp:1215/1562` AF kernel, `.../genotype.cpp` priors kernel | priors are *inputs*; no allele-length branch inside (`grep heterozygosity` → 0 hits) | — | agree by construction | n/a |

Latent edge (speculation, not observed): native compares **string sizes** (`candidate_reference(c).size() !=
candidate_alternate(c).size()`, `group.reference.size() == 1`) where GATK compares htsjdk `length()`.
For symbolic ALT strings (`"<NON_REF>"` = 7 chars, `"<DEL>"` = 5) these differ from htsjdk's `0`
and land on the same side of the comparison *unless* the record's REF happens to have length 7 or 5.
HC's EventMap only produces a concrete or `*` ALT, so this looks unreachable in HC (speculation); it
is reachable in principle in the GenotypeGVCFs builder (#12).

Gate-script convention note: the previous round registered its oracle in
`fastgatk-native/CMakeLists.txt` (+10 lines). This round is forbidden to edit that file, so the new
gate is **not** a CTest entry: `ctest -N` still prints `Total Tests: 285`, and the oracle must be run
directly (§D.2). This is the one deliverable I could not satisfy in the previous round's style.

---

## B. Reproduced divergence: gVCF reference-confidence `*` AF pseudocount

### B.1 Fixture (deliverable B)

New strict gate: **`fastgatk-native/scripts/verify_hc_gvcf_symbolic_prior_gatk_oracle.py`**
(modelled on `verify_hc_span_del_qual_oracle.py`: `main() -> int`, `--native` /
`FASTGATK_HC_BINARY`, `TemporaryDirectory` scratch, `--case`, `--threads`, `--expect-divergence`,
final JSON status line, non-zero exit on mismatch).

Reference/reads are the *same* synthetic fixture the previous round used (so the two rounds are
comparable): 1500 bp `chr1`, `random.seed(11)`, 10×`A` homopolymer at 600, 8×`CAG` at 900, base
`chr1:698` = `C` inside the `698-700 CCC` run; 30 reference reads (300M from 440, step 3), 12 reads
carrying `chr1:698 C>A`, 6–10 reads carrying a 1 bp deletion in the run (two encodings: sequence
removal with a 299M CIGAR, and a real `138M1D161M`).

Difference from the previous round: `-ERC BP_RESOLUTION` (reference confidence) so that the
`<NON_REF>` allele — and, at a site covered by the assembled deletion, the symbolic `*` allele — enter
the AF pseudocount vector. Interval `chr1:500-780`, `--min-pruning 1`, no `--alleles`; native
additionally `--threads 2`.

```
timeout 1800 python3 fastgatk-native/scripts/verify_hc_gvcf_symbolic_prior_gatk_oracle.py
timeout 1800 python3 fastgatk-native/scripts/verify_hc_gvcf_symbolic_prior_gatk_oracle.py --expect-divergence
```

Gated cases: `bp-span-del-frameshift`, `bp-span-del-cigar`, `bp-span-del-deep` (QUAL divergence
expected < / = / cleared) and three controls that must be byte-identical with no normalization at
all (`bp-snp-only-control`, `bp-insertion-only-control`, `bp-ref-only-control`). Two plain
`-ERC GVCF` probes are run and reported as `out_of_scope_notes` because native additionally writes a
redundant `END=` there (§F(a)).

### B.2 Literal rows

**Pre-fix** (binary from `4d66251`, `.diag/gvcf_prior_prefix.log`, strict run **exit 1**):

```
GATK   chr1 698 . C *,A,<NON_REF> 282.04 . BaseQRankSum=0.000;DP=48;ExcessHet=0.0000;MLEAC=0,1,0;MLEAF=0.00,0.500,0.00;MQRankSum=0.000;RAW_MQandDP=172800,48;ReadPosRankSum=-0.533  GT:AD:DP:GQ:PGT:PID:PL:PS:SB  2|0:30,6,12,0:48:99:1|0:697_TC_T:400,262,1204,0,848,1235,482,1284,1177,1595:697:30,0,18,0
NATIVE chr1 698 . C A,*,<NON_REF> 291.07 . BaseQRankSum=0.000;DP=48;ExcessHet=0.0000;MLEAC=1,0,0;MLEAF=0.500,0.00,0.00;MQRankSum=0.000;RAW_MQandDP=172800,48;ReadPosRankSum=-0.533  GT:AD:DP:GQ:PGT:PID:PL:PS:SB  0|1:30,12,6,0:48:99:1|0:697_TC_T:400,0,1235,262,848,1204,482,1177,1284,1595:697:30,0,18,0
```

Deep-deletion case (10 deletion reads, independent magnitude):

```
GATK   chr1 698 . C *,A,<NON_REF> 144.04 . BaseQRankSum=0.000;DP=52;ExcessHet=0.0000;MLEAC=0,1,0;MLEAF=0.00,0.500,0.00;MQRankSum=0.000;RAW_MQandDP=187200,52;ReadPosRankSum=0.120 ...
NATIVE chr1 698 . C A,*,<NON_REF> 153.07 . BaseQRankSum=0.000;DP=52;ExcessHet=0.0000;MLEAC=1,0,0;MLEAF=0.500,0.00,0.00;MQRankSum=0.000;RAW_MQandDP=187200,52;ReadPosRankSum=0.120 ...
```

After canonicalizing the symbolic ALT order (§B.4) the **only** remaining field difference in all
three span-del cases was `QUAL` (and nothing else):

```
#   - bp-span-del-frameshift: record 198 (chr1:698): QUAL '282.04' -> '291.07'
#   - bp-span-del-cigar:      record 198 (chr1:698): QUAL '282.04' -> '291.07'
#   - bp-span-del-deep:       record 198 (chr1:698): QUAL '144.04' -> '153.07'
first divergent row field diff (GATK -> NATIVE): QUAL '282.04' -> '291.07'
```

The `chr1:697 TC>T,<NON_REF>` record and all 279 reference rows were already byte-identical.

**Post-fix** (rebuilt binary, `.diag/gvcf_prior_final_strict.log`, strict run **exit 0**,
`"status": "pass"`):

```
GATK   chr1 698 . C *,A,<NON_REF> 282.04 ...
NATIVE chr1 698 . C A,*,<NON_REF> 282.04 ...
data_rows_byte_identical=True  raw_qual_mismatches=0
```

and the deep case `NATIVE ... 144.04` likewise. The plain-`-ERC GVCF` probes lost their QUAL
mismatch too (`raw QUAL mismatches []`), leaving only the unrelated `END=` note.

### B.3 Root cause and localization

`hc_call.cpp`, gVCF candidate-site AF block (`std::string gvcf(...)`), pre-fix:

```cpp
if (include_spanning_deletion) {
    const auto spanning = group.candidates.size() + 1U;
    prior_pseudocounts[spanning] =
        result.genotype_indel_heterozygosity * ref_pseudocount;   // <- fixed constant
}
```

GATK rule: `AlleleFrequencyCalculator.java:175-176` (`a.length() == refLength ? snpPseudocount :
indelPseudocount`), with `Allele.SPAN_DEL.length() == 1` (measured, §0). The `*` allele therefore
takes the **SNP** pseudocount on a 1 bp REF record. Passing the indel pseudocount changes the
Dirichlet posterior for every genotype containing `*` and hence
`P(no variant) = AFresult.log10ProbOnlyRefAlleleExists()`, which is exactly the gVCF record's QUAL
(`GenotypingEngine.java:158-161`). The sibling initializer in the ordinary-VCF path
(`calling_pipeline.cpp:17320`) had already been fixed in `4d66251`; the gVCF writer constructs its
own pseudocount vector and was missed. `<NON_REF>` on the next lines is **correct** in native
(htsjdk length 0 → indel prior).

### B.4 Gate criterion (what is and is not gated)

* Every data row is compared field by field after a canonicalization applied **identically to both
  files**: concrete ALTs sorted, `<NON_REF>` pinned last, `AC/AF/MLEAC/MLEAF` (Number=A), `AD`
  (Number=R), `PL` (Number=G) permuted with it, and phased genotype strings written in ascending
  allele order. `QUAL` is additionally re-checked on the **raw, un-normalized** rows of every record
  of every case — that check alone fails pre-fix.
* Records whose ALT **set** differs are violations, never normalized; any comma-valued INFO/FORMAT
  key outside the documented whitelist makes a record non-normalizable (so it must be identical).
  Only one record in the whole suite needs normalization (`chr1:698`, 1 normalization per span-del
  case).
* The `PGT`/`PID`/`PS` triple is *reported but not gated* on a record whose symbolic ALT order
  differs: native's `PGT` there encodes the pre-reorder allele index of `*` (GATK `0|1` vs native
  `0|2` after canonicalization), which is the separately tracked symbolic-ALT ordering gap (§F(b)),
  not a prior-selection defect. Every other field, and the phasing triple on every other record, is
  gated strictly.
* Two other divergences found while probing are deliberately **not** in the gate (they are unrelated
  to this class and are reported here instead): §F(c) zero-QUAL serialization, §F(d) hom-ref block
  AD/PL for deletion-carrying reads, §F(e) haploid `-ERC` concrete-vs-block emission. Gated cases use
  `-ERC BP_RESOLUTION`, which keeps them out of the comparison.

### B.5 Diff (the change the tree contains)

```diff
diff --git a/fastgatk-native/src/hc_call.cpp b/fastgatk-native/src/hc_call.cpp
@@ -4819,8 +4819,22 @@ std::string gvcf(const fastgatk::io::HtsReader& reader,
                 if (include_spanning_deletion) {
                     const auto spanning = group.candidates.size() + 1U;
+                    // AlleleFrequencyCalculator selects every allele's
+                    // pseudocount from its length:
+                    // `a.length() == refLength ? snpPseudocount : indelPseudocount`
+                    // (AlleleFrequencyCalculator.java:175-176, refLength =
+                    // vc.getReference().length()).  htsjdk's symbolic spanning
+                    // deletion is a 1 bp allele (`Allele.SPAN_DEL.length() ==
+                    // 1`), so on a 1 bp REF record it takes the SNP prior and
+                    // only a longer REF gives it the indel prior.  The
+                    // ordinary-VCF writer already derives it this way; the
+                    // unconditional indel prior here changed P(no variant) for
+                    // REF/* genotypes and therefore the reference-confidence
+                    // QUAL of any gVCF record that carries `*`.
                     prior_pseudocounts[spanning] =
-                        result.genotype_indel_heterozygosity * ref_pseudocount;
+                        (group.reference.size() == 1
+                             ? result.genotype_snp_heterozygosity
+                             : result.genotype_indel_heterozygosity) * ref_pseudocount;
                 }
```

`git diff --stat`: `fastgatk-native/src/hc_call.cpp | 16 +++++++++++++++-` (1 file, +15/-1);
plus the new file `fastgatk-native/scripts/verify_hc_gvcf_symbolic_prior_gatk_oracle.py` (untracked).

---

## C. Fixes

### C.1 Landed (reproduced divergence #3)

One ternary in `hc_call.cpp` (§B.5). Minimal: no EventMap/PairHMM restructure, no CMakeLists edit, no
Mutect2 change, no root `*.md` edit.

### C.2 Localized, NOT fixed (unreproduced) — recommendation per item

Each of these is the same class (GATK selects from allele length/symbolic-ness; native uses a fixed
or differently-derived value) but I could not build an input where the divergence becomes an emitted
byte difference, so per the brief they were stopped at localization.

1. **#5 `hc_call.cpp:3358-3360`** — arbitrary-ploidy ordinary-VCF writer. Ready-to-apply fix:
   ```cpp
   if (has_spanning_deletion)
       prior_pseudocounts.back() =
           (candidate_reference(calls.front()->candidate).size() == 1
                ? result.genotype_snp_heterozygosity
                : result.genotype_indel_heterozygosity) * ref_pseudocount;
   ```
   I *did* build and test this: the only observable probe (`--sample-ploidy 3`, the previous round's
   `reported-ploidy-3` case) prints QUAL `358.64` / QD `7.47` on both sides **before and after** the
   change (`.diag/ploidy3_after_gvcfix.log` vs `.diag/ploidy3_after_5.log`), i.e. the path is
   provably insensitive in the only reachable configuration, so the change could not be validated —
   it was reverted to keep this round's tree minimal and fully gated. Recommendation: land it together
   with a fixture that exercises an arbitrary-ploidy record where `*` is the only plausible extra
   allele, or accept it as source-derived only.
2. **#6 `hc_call.cpp:4485-4488`** — non-diploid / `hidden_spanning_deletion` gVCF writer, same fix
   shape. A haploid `-ERC BP_RESOLUTION` fixture was built and neither caller emits a `*` record there
   (native and GATK agree on the `*`-free rows; the only difference is §F(e), an emission/threshold
   branch). Recommendation: same as (1).
3. **#8 `hc_call.cpp:682-687`** and **#9 `calling_pipeline.cpp:412-414`** — `*` classified as
   symbolic/`OTHER` instead of GATK's called-non-symbolic 1 bp allele. Fix: in
   `prior_allele_type`, treat `"*"` as non-symbolic (then `reference.size() == alternate.size()`
   → SNP for a 1 bp REF), and in `genotype_priors_for_group` use `log10(snpHet) - log10(3)` for the
   spanning allele on a 1 bp REF. Both feed `GT`/`GQ` selection; in all six gated fixtures (including
   12-read and 52-read depths the `GT:AD:DP:GQ:PL` block was byte-identical after the allele-order
   canonicalization). Recommendation: needs a fixture where the `*` prior actually decides GT/GQ
   (e.g. a tie between the concrete ALT and `*`, `GQ` below the 99 cap) before changing it.
4. **#10 `calling_pipeline.cpp:6065-6066`** — `other_heterozygosity = max(snp, indel)` for the
   `<NON_REF>` allele of the 3-allele `[REF, concrete, <NON_REF>]` matrix used by
   `gvcf_candidate_is_non_monomorphic`; GATK's AF calculator gives `<NON_REF>` the **indel**
   pseudocount. Fix: replace `max(...)` with the indel heterozygosity. Recommendation: the retention
   flag only flips when the `<NON_REF>` absent posterior crosses `-0.1*stand_call_conf`; build a
   fixture at that boundary (very low depth on a gVCF ref-block/site boundary) first.
5. **#11 `hc_call.cpp:2924-2926`** — `estimate_mle_allele_counts(..., include_non_ref=true)` uses
   `max(snp, indel)` for the trailing `<NON_REF>`; fix is to use the indel heterozygosity (MLEAC/MLEAF
   can differ). Only reachable when the full symbolic AF matrix is unusable (the `else if
   (!include_spanning_deletion)` fallback), which none of my fixtures hit. Recommendation: localize
   with a unit-level call rather than a whole-pipeline fixture.
6. **#13 `genotype_gvcf_tool.cpp:1585-1615`** (GenotypeGVCFs) — same `"*"`→`Other` misclassification;
   out of this round's HC scope, not touched.

---

## D. Gate commands and exit statuses (deliverable D)

| # | Command | Result |
|---|---|---|
| D.1 | `cmake --build fastgatk-native/build --target fastgatk-hc-call -j 16` | exit 0 (`[100%] Built target`, 18.8 s incremental) |
| D.2 | `python3 fastgatk-native/scripts/verify_hc_gvcf_symbolic_prior_gatk_oracle.py` (new strict gate) | **exit 0**, `"status": "pass"`, 6/6 gated cases byte-identical, 0 raw-QUAL mismatches (0 `RAW QUAL MISMATCH` lines in `.diag/gvcf_prior_final_strict.log`) |
| D.2b | same gate, pre-fix binary | **exit 1** (recorded in `.diag/gvcf_prior_prefix.log`; the divergence literals in §B.2) |
| D.3 | `ctest --test-dir fastgatk-native/build -R 'window-invariance-gatk-oracle\|alleles-overlap-gate-oracle\|span-del-qual-gatk-oracle' -V` | **exit 0**, `100% tests passed, 0 tests failed out of 4` (`window-invariance`, `ploidy-window-invariance`, `alleles-overlap-gate`, `span-del-qual`, 238.4 s) — `.diag/ctest_gates_final.log` |
| D.4 | `python3 fastgatk-native/scripts/verify_hc_alleles_gatk_oracle.py` | **exit 0**, `"status": "pass"` — `.diag/alleles_oracle_final.log` |
| D.5 | `cmake --build fastgatk-native/build-serial --target fastgatk-hc-call -j 16` | exit 0 |
| D.6 | `fastgatk-native/scripts/run_regression.sh --label 'length-prior sweep: gVCF symbolic `*` AF prior fix'` | **exit 0** |
| D.7 | `FASTGATK_HC_BINARY=fastgatk-native/build-serial/fastgatk-hc-call python3 fastgatk-native/scripts/verify_hc_gvcf_symbolic_prior_gatk_oracle.py` | **exit 0**, `"status": "pass"`, 0 raw-QUAL mismatches — the new gate is backend-independent (`.diag/gvcf_prior_final_serial.log`) |

### D.6.1 Full dual-backend regression (outdir `.diag/regression/20260911-064512`)

```
- 时间：2026-09-11 07:02:39 CST
- 标签：length-prior sweep: gVCF symbolic `*` AF prior fix
- git：`4d66251`（未提交变更 2 项）
- 过滤：`<全量>`　并行度：8

| 后端   | 构建目录                             | 结果 | 通过/总数 | 耗时          |
| ------ | ------------------------------------ | ---- | --------- | ------------- |
| omp    | `OpenMP (fastgatk-native/build)`     | 通过 | 285/285   | 1040.37 sec   |
| serial | `Serial (fastgatk-native/build-serial)` | 通过 | 285/285 | 1047.45 sec   |
```
`omp.status` = `PASS`, `serial.status` = `PASS`, script exit 0. So **285/285 on both backends**,
unchanged suite size (285), no test weakened or skipped.

Comparator self-test (sanity of the new gate's normalizer, run as a throwaway `python3` snippet):
a hand-permuted `C A,*` / `C *,A` record with correctly permuted `PL`/`AD`/`MLEAC`/`MLEAF`/`GT`
passes with exactly one reported normalization, while a single wrong `PL` entry, an unpermuted
`MLEAC`, an ALT **subset**, and an unknown comma-valued INFO key are each reported as a violation.


---

## E. Does the tree contain a change?

**CONTAINS.** Two paths:

```
$ git status --short
 M fastgatk-native/src/hc_call.cpp
?? fastgatk-native/scripts/verify_hc_gvcf_symbolic_prior_gatk_oracle.py
```

No file was reverted: every acceptance criterion in §D passed with the change in place. The change is
exactly one ternary plus its comment in the gVCF writer; nothing else under
`fastgatk-native/src|include` was modified (verified with `git diff --stat`). `git diff` is reproduced
verbatim in §B.5.

---

## F. Adjacent divergences found while probing (NOT this class, NOT fixed, NOT gated)

Reported because they were measured on the same fixtures and a future round should not rediscover
them from scratch. All are pre-existing (none is affected by the §B.5 change) and none is in the gate.

**(a) native-only redundant `END=` on concrete-indel records in plain `-ERC GVCF`** — measured on
`gvcf-span-del-frameshift` (`chr1:697`):

```
GATK   chr1 697 . TC T,<NON_REF> 62.60 . BaseQRankSum=0.000;DP=48;...
NATIVE chr1 697 . TC T,<NON_REF> 62.60 . END=698;BaseQRankSum=0.000;DP=48;...
```
`END=698` equals the value implied by `POS+len(REF)-1`, so htsjdk omits it; native writes it whenever
`!bp_resolution && has_indel` (`hc_call.cpp:4961-4964`). This is why the gate's gated cases use
`-ERC BP_RESOLUTION` and the two plain-GVCF probes are reported as `out_of_scope_notes`.

**(b) symbolic `*` ALT ordering / `PGT` encoding on the same record** — native always renders `*` as
the last concrete ALT, GATK places it by haplotype order:

```
GATK   ... C *,A,<NON_REF> ... GT:...:PGT:...:PL:...  2|0:30,6,12,0:48:99:1|0:...
NATIVE ... C A,*,<NON_REF> ... GT:...:PGT:...:PL:...  0|1:30,12,6,0:48:99:1|0:...
```
After canonicalization `GT`, `AD`, `MLEAC`, `MLEAF`, `PL` all agree; native's `PGT` does **not**
(GATK `0|1` = the phased haplotype carries `*`; native `0|2` = it carries the SNV), i.e. native's
`PGT` string encodes the pre-reorder allele index of `*`. This is the same family as the already
tracked symbolic-ALT emission gap (the previous round's ploidy-3 note), so it is reported, not gated.

**(c) zero-QUAL serialization in `-ERC` records** — probe `shallow-ref8-snp2-del4` (8 ref + 2 SNV +
4 CIGAR-deletion reads), `chr1:698`: GATK `... 0 . ...`, native `... 0.00 . ...`. htsjdk trims a
zero QUAL; native prints fixed two decimals. Unrelated to priors (the underlying posterior is ~0 on
both sides).

**(d) hom-ref block AD/PL for reads carrying a deletion inside a repeat** — probe
`shallow-ref5-del5-nosnp` (5 ref + 5 CIGAR-deletion reads, `-ERC BP_RESOLUTION`):

```
GATK   chr1 698 . C <NON_REF> . . .  GT:AD:DP:GQ:PL  0/0:5,5:10:0:0,0,50
NATIVE chr1 698 . C <NON_REF> . . .  GT:AD:DP:GQ:PL  0/0:10,0:10:30:0,30,448
```
Both callers agree at `chr1:700` (`0/0:5,5:10:0:0,0,50`), so the disagreement is the classification
of deletion-carrying reads at a repeat-ambiguous position, not a prior. Out of class, sizeable, not
touched.

**(e) haploid `-ERC BP_RESOLUTION` concrete-vs-block emission** — probe `haploid-bp-span-del`:
GATK emits a `<NON_REF>` HomRefBlock at `chr1:697`, native emits a concrete record
`chr1 697 . TC T,<NON_REF> 0.00 ... GT:AD:DP:GQ:PL:SB 0:42,6,0:48:99:0,1341,1555:42,0,6,0`. This is
an emission/threshold branch (native's biallelic `confidence_log10_p_alt_absent` test at
`calling_pipeline.cpp:17497-17501` vs GATK's `siteIsMonomorphic` computed over the `[REF, ALT,
<NON_REF>]` VC) rather than an allele-length prior, and fixing it means computing the absent
posterior over the symbolic matrix — not minimal, not attempted.

---

## G. What remains unproven

* **The `*` prior fix is proven only for the `-ERC` reference-confidence path with a concrete
  ALT present at the same site** (three gated cases, two independent depths, two deletion encodings)
  and only to the printed precision of the compared rows (`QUAL` 2 decimals; the PL/AD rows are
  byte-identical, so the underlying integer matrices agree exactly). The internal `double` posterior
  is not proven bit-identical.
* Sites #5, #6, #8, #9, #10, #11, #13 are **source-level divergences without a fixture**. #5 was
  empirically shown to be insensitive in the only reachable configuration (triploid), which is *not*
  proof that it is harmless elsewhere; #8/#9/#10/#11 are pure source reasoning.
* The `size()` vs htsjdk `length()` equivalence for symbolic alleles (§A, latent edge) is reasoning,
  not a measured divergence; it is labelled speculation.
* The plain-`-ERC GVCF` mode is only *probed*, not gated, because of the native-only `END=`
  divergence (§F(a)); its `chr1:698` QUAL does now match GATK.
* The gate is not registered in CTest (CMakeLists editing is forbidden this round), so `ctest` alone
  does not protect this rule; anyone running only `ctest` would miss a regression of it.
