# Track B — `--alleles` / GenotypeGivenAlleles injection boundary (findings)

Date: 2026-09-11. **Audit only** — no production source or CMake file was modified;
no `git commit`. Two new oracle scripts were created (both listed in (c)).
Binaries: `fastgatk-native/build/fastgatk-hc-call` (Kokkos/OpenMP) and
`fastgatk-native/build-serial/fastgatk-hc-call` (Serial). Oracle: pinned GATK 4.6.2.0
(`third_party/jdk17/bin/java -Xmx1g -jar
third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar HaplotypeCaller`).
All scratch under `TemporaryDirectory`; nothing written into the repo tree.

**Bottom line: the given-alleles injection boundary is NOT 1:1 with GATK 4.6.2.0.**
Smallest counterexample observed (D2): *reference-only reads + two `--alleles` records in one
region — a deletion and a SNP whose position lies inside that deletion's REF span.* GATK emits
two VCF rows; native emits one.

---

## (a) Exact commands and pass/fail matrix

Each script builds one shared synthetic fixture per run (1500 bp `chr1`, 300 bp 60MQ reads,
a 10×`A` homopolymer at 600, an 8×`CAG` tandem repeat at 900, and 6 reads carrying an assembled
1 bp deletion near 700) and then, per case, runs GATK and the native binary with **identical**
arguments, comparing VCF data rows field-by-field:

```
java  -Xmx1g -jar <gatk.jar> HaplotypeCaller -R ref.fa -I input.bam -L <interval> \
      [--alleles feature.vcf.gz] --min-pruning 1 --create-output-variant-index false \
      --add-output-vcf-command-line false [<case extra args>] -O gatk.<case>.vcf
<native> -R ref.fa -I input.bam -L <interval> [--alleles feature.vcf.gz] --min-pruning 1 \
      --create-output-variant-index false --add-output-vcf-command-line false \
      [<case extra args>] --threads 2 -O native.<case>.vcf
```

### Pass/fail matrix (data-row exactness, GATK 4.6.2.0 vs native)

| # | Case | Command (both backends unless noted) | OMP `build/` | Serial `build-serial/` |
|---|------|--------------------------------------|--------------|------------------------|
| 0 | baseline SNP/multi-allelic/indel/gVCF | `python3 fastgatk-native/scripts/verify_hc_alleles_gatk_oracle.py` | **PASS** (exit 0) | not run (pre-existing test; backend-independent result) |
| 1 | `homopolymer-insertion-left` | `verify_hc_alleles_deep_boundary.py` | PASS | PASS |
| 2 | `homopolymer-insertion-right` | " | PASS | PASS |
| 3 | `homopolymer-deletion` | " | PASS | PASS |
| 4 | `tandem-repeat-deletion-anchored` | " | PASS (0 rows on both sides) | PASS (0 rows) |
| 5 | `tandem-repeat-one-unit-literal` | " | **DIVERGE (D3)** java exit 3 / native exit 0 | **DIVERGE (D3)** |
| 6 | `assembled-event-overlap` | " | **DIVERGE (D1)** | **DIVERGE (D1)** |
| 7 | `adjacent-snps` | " | PASS | PASS |
| 8 | `overlapping-features` | " | **DIVERGE (D2)** | **DIVERGE (D2)** |
| 9 | `region-boundary-span` | " | PASS | PASS |
| 10 | controls: cases 5, 6, 8 with `--drop-alleles` | `--drop-alleles` | PASS (case 6: 1 identical row both sides; cases 5 and 8: 0 rows both sides, java exit 0) | not run (attribution control, backend-independent) |
| 11 | `max-alt-alleles-1` | `verify_hc_alleles_deep_limits.py` | PASS | PASS |
| 12 | `max-alt-alleles-2` | " | PASS | PASS |
| 13 | `max-genotype-count-3` | " | PASS | PASS |
| 14 | `max-genotype-count-2` | " | **DIVERGE (O2, not alleles-specific)** | **DIVERGE (O2)** |
| 15 | control: case 14 with `--drop-alleles` | `--drop-alleles` | **still diverges** → O2 is *not* a `--alleles` divergence | not run |
| 16 | `gvcf-homopolymer-insertion` | " | **DIVERGE (O1, not alleles-specific)** | **DIVERGE (O1)** |
| 17 | `gvcf-max-alt-alleles-1` | " | **DIVERGE (O1 only; rows 0-4 exact)** | **DIVERGE (O1)** |
| 18 | control: case 16 with `--drop-alleles` | `--drop-alleles` | **still diverges at the same INFO field** → O1 is *not* a `--alleles` divergence | not run |

Exact invocation used for the two backends (`FASTGATK_REQUIRE_GATK_ORACLE` unset; both scripts skip
cleanly if the jar/binary is absent):

```
timeout 900 python3 fastgatk-native/scripts/verify_hc_alleles_deep_boundary.py
timeout 900 python3 fastgatk-native/scripts/verify_hc_alleles_deep_boundary.py \
    --native fastgatk-native/build-serial/fastgatk-hc-call
timeout 900 python3 fastgatk-native/scripts/verify_hc_alleles_deep_limits.py
timeout 900 python3 fastgatk-native/scripts/verify_hc_alleles_deep_limits.py \
    --native fastgatk-native/build-serial/fastgatk-hc-call
timeout 600 python3 fastgatk-native/scripts/verify_hc_alleles_deep_boundary.py --drop-alleles \
    --case assembled-event-overlap --case overlapping-features --case homopolymer-insertion-left
timeout 600 python3 fastgatk-native/scripts/verify_hc_alleles_deep_limits.py --drop-alleles \
    --case gvcf-homopolymer-insertion
timeout 400 python3 fastgatk-native/scripts/verify_hc_alleles_deep_limits.py --drop-alleles \
    --case max-genotype-count-2
timeout 400 python3 fastgatk-native/scripts/verify_hc_alleles_deep_boundary.py --drop-alleles \
    --case tandem-repeat-one-unit-literal
```

Serial and OpenMP produced **identical** divergence sets in every case
(proof that none of these are a parallel-reduction artifact); the two backends were
not byte-compared against each other, only each against GATK.

---

## (b) Divergences

### D1 — a forced indel that is equivalent-but-not-identical to an assembled event is injected, then never emitted

**Repro:** `verify_hc_alleles_deep_boundary.py --case assembled-event-overlap`
(read set includes 6 reads carrying an assembled 1 bp deletion; feature VCF has the single
record `chr1 700 . CA C . PASS .`, interval `chr1:500-780`).

* first differing row: index 1 (native has no second row). **FIRST DIFFERING FIELD: ROW_COUNT**
  (row 0 is byte-identical, including `PL=106,0,961`).

```
GATK   chr1 697 . TC T 98.60 . AC=1;AF=0.500;AN=2;...;MLEAC=1;MLEAF=0.500;...  GT:AD:DP:GQ:PL 0/1:30,6:36:99:106,0,961
GATK   chr1 700 . CA C 144.77 . AC=0;AF=0.00;AN=2;DP=36;ExcessHet=0.0000;FS=0.000;MLEAC=0;MLEAF=0.00;MQ=60.00;SOR=0.001  GT:AD:DP:GQ:PL 0/0:36,0:36:99:0,108,1481
NATIVE chr1 697 . TC T 98.60 . AC=1;AF=0.500;AN=2;...;MLEAC=1;MLEAF=0.500;...  GT:AD:DP:GQ:PL 0/1:30,6:36:99:106,0,961
                            <no chr1:700 row>
```

**Control (proven cause):** with `--drop-alleles` both sides emit exactly the one `chr1:697` row,
so the missing `chr1:700` row is created by the `--alleles` injection path, not by general
candidate discovery.

**Native-side trace** (`FASTGATK_DEBUG_EVENTMAP=1 FASTGATK_DEBUG_EVENTMAP_LIFECYCLE=1`, OMP build):

```
[FASTGATK_EVENTMAP_INJECT]  path=1 tid=0 position=696 ref=TC alt=T  appended=1
[FASTGATK_EVENTMAP_INJECT]  path=2 tid=0 position=699 ref=CA alt=C  appended=1
[FASTGATK_EVENTMAP_INJECT]  path=3 tid=0 position=698 ref=CCA alt=C appended=1
[FASTGATK_EVENTMAP_REGION]  ... position=696 ... calling_class=1 ; 698 ... calling_class=1 ; 699 ... calling_class=1
[FASTGATK_EVENTMAP_CALL]    tid=0 position=696 ref=TC alt=T qual=98.6018      <-- only this one
```

So the injection itself works (three artificial paths were appended and all three regenerated
EventMap events were classified `calling_class=1`), and the loss happens **between candidate
classification and `EVENTMAP_CALL`/VCF emission**.

**Hypothesis (source-level).** GATK `HaplotypeCallerGenotypingEngine.genotypeVariants()`
(`HaplotypeCallerGenotypingEngine.java:190`) builds one merged `VariantContext` per locus from
`readLikelihoods.alleles()` and emits *every* locus, so the injected `CA>C` haplotype always
contributes its own site row even when `AC=0` (`GenotypeGivenAlleles` + `calculateOutputAlleleSubset`,
`GenotypingEngine.java:144-176`). Native reduces the region's `result.candidates` in
`calling_pipeline.cpp` (the EventMap/`calling_class` reduction feeding the per-region call loop,
`[FASTGATK_EVENTMAP_CALL]` vs `[FASTGATK_EVENTMAP_REGION]` sites in
`fastgatk-native/src/calling_pipeline.cpp`) with a span-overlap suppression that keeps only one
of several mutually overlapping candidate events in the same region
(cf. `multiallelic_suppressed` grouping near `calling_pipeline.cpp:16470-16485` and the
`restrict_to_forced_alleles` / HC ownership gate near `calling_pipeline.cpp:16163-16170`).
The injected alleles are present at the `AssemblyResult` boundary — the divergence is downstream,
in the *emission* of injected events, not in `inject_hc_given_alleles_into_haplotype_paths()`
(`calling_pipeline.cpp:12797`).

---

### D2 — **smallest counterexample**: a forced SNP inside the REF span of a forced deletion is injected and classified, but never emitted

**Repro:** `verify_hc_alleles_deep_boundary.py --case overlapping-features`.
Reference-only coverage over an 8×`CAG` tandem repeat; feature VCF (one file, two records):

```
chr1 900 . CAGCAG C . PASS .
chr1 904 . A      T . PASS .
```

* first differing row: index 1. **FIRST DIFFERING FIELD: ROW_COUNT** (row 0 identical).

```
GATK   chr1 900 . CAGCAG C 126.81 . AC=0;AF=0.00;AN=2;DP=30;...;MLEAC=0;MLEAF=0.00;MQ=60.00;SOR=0.001  GT:AD:DP:GQ:PL 0/0:30,0:30:90:0,90,1350
GATK   chr1 904 . A      T 117.78 . AC=0;AF=0.00;AN=2;DP=30;...;MLEAC=0;MLEAF=0.00;MQ=60.00;SOR=0.001  GT:AD:DP:GQ:PL 0/0:30,0:30:90:0,90,1343
NATIVE chr1 900 . CAGCAG C 126.81 . AC=0;AF=0.00;AN=2;DP=30;...;MLEAC=0;MLEAF=0.00;MQ=60.00;SOR=0.001  GT:AD:DP:GQ:PL 0/0:30,0:30:90:0,90,1350
                            <no chr1:904 row>
```

**Control (proven cause):** `--drop-alleles` → both sides emit **zero** rows. The two rows are
purely `--alleles`-driven, so the lost `chr1:904` row is an injection-boundary divergence.

**Native-side trace** (OMP build):

```
[FASTGATK_EVENTMAP_MATERIALIZED] path=1 position=899 ref=CAGCAG alt=C
[FASTGATK_EVENTMAP_MATERIALIZED] path=2 position=903 ref=A      alt=T
[FASTGATK_EVENTMAP_INJECT]  path=1 position=899 ref=CAGCAG alt=C appended=1
[FASTGATK_EVENTMAP_INJECT]  path=2 position=903 ref=A      alt=T appended=1
[FASTGATK_EVENTMAP_REPLAY]  position=899 inserted=0 reason=already-present
[FASTGATK_EVENTMAP_REPLAY]  position=903 inserted=0 reason=already-present
[FASTGATK_EVENTMAP_REGION]  position=899 ... calling_class=1
[FASTGATK_EVENTMAP_REGION]  position=903 ... calling_class=1
[FASTGATK_EVENTMAP_CALL]    position=899 qual=126.814                       <-- 903 lost here
```

Both forced events were injected into the haplotype population and both survived to
`calling_class=1`; only the one that is not inside another candidate's REF span is emitted.

**Hypothesis (source-level).** GATK's `AssemblyResultSet.addGivenAlleles()`
(`AssemblyResultSet.java:736`) injects each given `Event` independently into every base haplotype
whose EventMap has no *overlapping* event (`makeHaplotypeWithInsertedEvent`,
`AssemblyResultSet.java:774-790`), so the `904 A>T` haplotype exists beside the `900 CAGCAG>C`
haplotype, and `HaplotypeCallerGenotypingEngine` emits both loci (they are different
`VariantContext`s). Native injects into the same bounded base-path population
(`inject_hc_given_alleles_into_haplotype_paths`, `calling_pipeline.cpp:12797`, whose `base_paths`
snapshot and per-path `events_overlap()` test mirror GATK's rule) — and the trace shows the
injection succeeding — but the two resulting loci are then reduced to one site by native's
region-span/output reduction before `[FASTGATK_EVENTMAP_CALL]`.
**Candidate native location:** the candidate/region reduction between EventMap classification and
per-region emission in `calling_pipeline.cpp` (same code region as D1; the
`calling_class`/`multiallelic_suppressed` and HC ownership gates around lines 16163-16170 and
16470-16485). This is reported as a *location*, not a proven single line — see (d).

---

### D3 — a `--alleles` record whose `Event` minimal representation has an empty ALT: GATK aborts, native calls it

**Repro:** `verify_hc_alleles_deep_boundary.py --case tandem-repeat-one-unit-literal`, feature row
`chr1 903 . CAGCAG CAG . PASS .`

```
GATK   exit 3   java.lang.IllegalArgumentException: Null alleles are not supported   (no VCF written)
NATIVE exit 0   VCF written with zero data rows
```

**Control:** the same region with `--drop-alleles` → GATK exit 0, native exit 0, zero rows on both
sides. The abort is caused by this `--alleles` record alone (it is not a general region/argument
problem).

**Hypothesis (source-level, proven by reading).** `Event`'s constructor calls
`makeMinimalRepresentation` (`gatk-source/.../utils/haplotype/Event.java:33-68`); for
`ref="CAGCAG", alt="CAG"` the 3-base suffix overlap is removed, leaving `alt = ""`, and
`Allele.create(empty, false)` throws `IllegalArgumentException: Null alleles are not supported`
(this is exactly the exception observed, and
`AssemblyResultSet.isEventPresentInAssembly` relies on that minimalization: *"note that Events are
forced to have a minimal representation"*). Native's `load_forced_alleles()`
(`fastgatk-native/src/hc_call.cpp:2213`) stores the raw `ref`/`alt` strings and never applies
`makeMinimalRepresentation`, so `inject_hc_given_alleles_into_haplotype_paths()` deletes the
3 bases literally instead of rejecting the record. Consequence: on this input class native
**accepts and exits 0** where GATK aborts the whole tool run; native itself emits no data rows
here, so this is a *failure-mode/exit-status* divergence rather than an extra-call divergence —
but a 1:1 boundary would require native to reject the record the way GATK does.

---

### O1 / O2 — out-of-track observations that the controls prove are NOT `--alleles` divergence

* **O1 (gVCF indel `END=`).** Case `gvcf-homopolymer-insertion` (and the gVCF case of the limits
  script) diverges only on row index 6, first differing field = INFO:
  GATK `BaseQRankSum=0.000;DP=36;...` vs native `END=698;BaseQRankSum=0.000;DP=36;...`.
  All GT/AD/DP/GQ/PL values and every `<NON_REF>` block boundary are identical, and the sequence
  of rows is identical (37 vs 37 / 42 vs 42). **The identical divergence occurs with
  `--drop-alleles`**, so this is a plain gVCF indel-symbolic-END annotation difference, unrelated
  to `--alleles`. Note also that a forced *homopolymer insertion* produces **no** forced ALT row in
  gVCF on either side (both sides agree) — the baseline `verify_hc_alleles_gatk_oracle.py` gVCF case
  (forced SNP) still passes.
* **O2 (`--max-genotype-count 2`).** GATK 4.6.2.0 aborts (exit 3) with
  `java.lang.IllegalArgumentException: VariantContext at chr1:620 has only a single reference allele,
  but getLog10PNonRef requires at least alternate allele`; native returns a call. Reproduced with
  `--drop-alleles` too, so it is a general GATK `--max-genotype-count 2` defect (ref-only VC reaching
  `getLog10PNonRef` via `removeAltAllelesIfTooManyGenotypes`/`removeExcessAltAllelesFromVC`,
  `HaplotypeCallerGenotypingEngine.java:413-437`), **not** an injection-boundary divergence.
  For the Java-accepted `--max-genotype-count 3` with three forced SNPs the result is exact.

---

## (c) Scripts and their CMake registration (snippet only — NOT applied)

### `fastgatk-native/scripts/verify_hc_alleles_deep_boundary.py`
Pins the `AssemblyResultSet.addGivenAlleles()` structural boundary: homopolymer indel placement
(leftmost / rightmost / deletion), one-CAG-unit deletion in an 8×`CAG` tandem repeat (anchored form
and the literal `ref=2 units/alt=1 unit` form), a `--alleles` REF span crossing the `-L` interval
end, an already-assembled-equivalent forced indel (D1), two adjacent forced SNPs, two mutually
overlapping forced events (D2), plus the region-wide baseline SNP path.
`--case NAME` selects cases; `--drop-alleles` is the attribution control; `--native PATH`
selects the backend; non-zero exit on divergence.

### `fastgatk-native/scripts/verify_hc_alleles_deep_limits.py`
Pins `--alleles` crossed with the allele budget: `--max-alternate-alleles 1/2` with three forced
SNPs, `--max-genotype-count 3` and `2`, and two gVCF (`-ERC GVCF`) combinations. Imports the
fixture helpers from the boundary script so both share one reference/BAM construction.

```cmake
    add_test(NAME fastgatk-hc-alleles-deep-boundary-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_alleles_deep_boundary.py")
    set_tests_properties(fastgatk-hc-alleles-deep-boundary-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    add_test(NAME fastgatk-hc-alleles-deep-limits-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_alleles_deep_limits.py")
    set_tests_properties(fastgatk-hc-alleles-deep-limits-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
```

(Place next to the existing `fastgatk-hc-complex-multiallelic-gatk-oracle` block at
`fastgatk-native/CMakeLists.txt:1384-1388`; the snippet above was **not** inserted.)
Both scripts support `--native`, so the serial tree can be exercised with
`FASTGATK_HC_BINARY=…/build-serial/fastgatk-hc-call` or `--native`, and exit non-zero while
D1/D2/D3 are open — they are intended to become **regression tests**, not to be registered green.

---

## (d) Proven vs not proven

**Proven (observed directly, with a control):**

1. D2 is a genuine `--alleles` divergence: GATK emits rows at `chr1:900` **and** `chr1:904`;
   native emits only `chr1:900`; `--drop-alleles` makes both sides emit zero rows
   (so the rows exist only because of `--alleles`). First differing field = ROW_COUNT.
2. D1 is a genuine `--alleles` divergence of the same shape: GATK emits `chr1:697` **and**
   `chr1:700` (`CA>C`, `AC=0`, `0/0`); native emits only `chr1:697`; `--drop-alleles`
   gives one identical row on both sides. First differing field = ROW_COUNT.
3. In both D1 and D2 the native *injection* step itself ran and produced `calling_class=1`
   EventMap events for every forced allele (native debug trace, quoted verbatim above);
   the divergence is downstream of the injection, between candidate classification and call/emission.
4. D3 is a genuine boundary divergence in the other direction: GATK aborts (exit 3) with
   `IllegalArgumentException: Null alleles are not supported` on this `--alleles` record while
   native exits 0 with an empty (header-only) VCF; the same command without `--alleles` gives
   GATK exit 0, so the abort is caused by that record. Observed on both backends.
5. Homopolymer placement (leftmost vs rightmost vs deletion), anchored tandem-repeat deletion,
   region-boundary-spanning forced deletion, adjacent forced SNPs, `--max-alternate-alleles 1/2`
   with three forced SNPs, and `--max-genotype-count 3` are **exact** — no divergence observed in
   these cases on either backend.
6. O1 and O2 reproduce with `--alleles` removed, so neither is attributable to the
   given-alleles boundary.

**Not proven / not exercised (be explicit):**

* **Base-haplotype score-ranking ties.** Not exercised. Every fixture here has at most one
  non-reference base haplotype, so the `NUM_HAPLOTYPES_TO_INJECT_FORCE_CALLING_ALLELES_INTO = 5`
  cut (`AssemblyBasedCallerUtils.java:58`, consumed at `AssemblyResultSet.java:745-750`) never
  discards anything and GATK's `.sorted(...)` stability vs native's `std::stable_sort`
  (`calling_pipeline.cpp:12855-12870`) is untested. A tie-breaking divergence remains **unproven,
  not refuted**.
* **`--max-genotype-count` × NaN injected-haplotype scores.** The suspicion (GATK's
  `AlleleScoredByHaplotypeScores.compareTo`, `HaplotypeCallerGenotypingEngine.java:487-502`, uses
  `>`/`<` on `Double` so a NaN `bestHaplotypeScore` from `Haplotype.insertAllele()`
  (`Haplotype.java:34` `score = Double.NaN`) loses every comparison, matching native's
  `haplotype_score > score.highest` handling in `subset_hc_alleles_by_haplotype_scores`,
  `calling_pipeline.cpp:3647-3740`) produced **no observable divergence** at
  `--max-genotype-count 3`: the output was byte-identical. The ranking rule is therefore
  *untested at the point where it would matter* (needs ≥ `computeMaxAcceptableAlleleCount` forced
  alleles plus competing non-reference haplotypes) — treated as unproven, not as verified.
* **Exact native code line for D1/D2.** The losing step is bracketed to the code between
  `[FASTGATK_EVENTMAP_REGION]` (classify) and `[FASTGATK_EVENTMAP_CALL]` (emit) in
  `calling_pipeline.cpp`; I did **not** instrument further, so the precise predicate/branch that
  drops the second site is a hypothesis, not a proven line.
* **Multi-sample / joint genotyping, `--force-call-filtered-alleles` with repeat indels,
  symbolic (`<DEL>`, `*`, `<NON_REF>`) feature ALTs, and Mutect2's separate
  `restrict_to_forced_alleles` replay** were not re-audited here; the baseline oracle's
  filtered-allele and multi-allelic cases still pass, and Track C owns the Mutect2 side.
* The serial column of rows 5, 13 in matrix (a) was re-run after the case rename/split and
  matched the OpenMP column; the serial `tandem-repeat-deletion` row predates that split and is
  the same input.

---

## Answer to the track question

**Is the given-alleles injection boundary 1:1 with GATK 4.6.2.0?** No.

**Smallest counterexample (D2).** Reference-only coverage over a tandem repeat, one `--alleles`
VCF containing two PASS records, interval covering both:

```
chr1 900 . CAGCAG C . PASS .
chr1 904 . A      T . PASS .
```

GATK 4.6.2.0 emits two rows (`chr1:900 CAGCAG>C 0/0` and `chr1:904 A>T 0/0`); native emits only
`chr1:900`. Removing `--alleles` makes both sides emit nothing, so the missing row is a
given-alleles-injection-boundary divergence, not general candidate discovery. Both forced events
are injected and classified in native (`[FASTGATK_EVENTMAP_INJECT]`, `calling_class=1`); the
second is lost before `[FASTGATK_EVENTMAP_CALL]`.
