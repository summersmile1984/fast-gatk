# Round 49 — how far does the GenotypeGVCFs reverse trim actually reach, and could the pre-existing corpus have caught it?

Date: 2026-09-13. Workspace: `/home/turing-agents/Documents/fast-gatk`.
Subject of the measurement: `apply_gatk_reverse_trim()` in
`fastgatk-native/src/genotype_gvcf_tool.cpp` (working-tree change, **not committed**:
`git diff --stat HEAD -- fastgatk-native/src/genotype_gvcf_tool.cpp` = `111 insertions`,
and `git show HEAD:...genotype_gvcf_tool.cpp | grep -c apply_gatk_reverse_trim` = `0`),
gated by `fastgatk-native/scripts/verify_genotype_gvcf_reverse_trim_gatk_oracle.py`.

Nothing under `fastgatk-native/src`, `fastgatk-native/include`, `fastgatk-kernels`,
any `CMakeLists.txt` or `third_party/` was modified. No `git add` / `git commit`.
No binary was rebuilt. Scratch lives in `.diag/round49/` (git-ignored).

## 0. Answer up front

* **Reach on the project's real corpora: zero records.**
  On **142** real (non-LFS) gVCFs — **165,527 records** — GATK's own
  `GATKVariantContextUtils.reverseTrimAlleles()` rewrites **0** records
  (`changed=0`, `errors=0`, measured by calling the pinned jar's method directly).
  The python predicate agrees: 0 fires.
* **Why: the `:1458` guard, and specifically the one-base REF.**
  165,524 of 165,527 records carry an allele of length 1 that is not `*`
  (**162,712** of them a one-base **REF** — every SNP and every one-base-REF indel).
  Only **3** records in the
  entire corpus even reach `normalizeAlleles()`, and all three share no trailing
  base. Even a counterfactual run with the `:1458` guard deleted fires **0** records.
* **Variant records that fire: 0** (0 of the 25,305 records that carry a concrete ALT).
* **Pre-existing registered corpus: it could not have caught this.**
  No registered gate other than the new one feeds GenotypeGVCFs a record of the
  trimmable shape. The two pre-existing registered fixtures that do contain one
  (`verify_left_align.py:52`, `verify_reblock_gvcf.py:422`) belong to
  LeftAlignAndTrimVariants and ReblockGVCF, neither of which calls
  `reverseTrimAlleles()`. The one GenotypeGVCFs fixture with the shape outside the
  new gate is the **unregistered** diagnostic gate added in the previous round
  (`verify_genotype_gvcf_dense_materialize_gatk_oracle.py:152`), whose `why` string
  names the divergence.
* **Practical significance: no row of the main real-data parity fixtures changes.**
  Measured: every emitted row of both real-data fixtures has a one-base REF
  (3/3 default rows at chr17, 1001/1001 chr17 dense rows, 810/810 chr20 dense rows),
  so the `:1458` guard short-circuits at every locus, and the project's own
  pre-fix log shows the real-data parity gate **passing before the fix**
  (`.diag/emitted-ownership-ctest-genotype-gvcf.log`, 2026-09-12T20:45,
  `fastgatk-genotype-gvcf-gatk-oracle ... Passed`, `100% tests passed, 0 tests failed out of 18`).

| question | measured answer |
| --- | --- |
| records the trim rewrites, pinned GATK gVCF corpora | 0 of 165,527 (142 corpora) |
| records that pass the `:1458` guard at all | 3 of 165,527 |
| of those, records that fire | 0 |
| records that fire and are variant records | 0 |
| predicates fire without the guard (counterfactual) | 0 |
| rows changed in the chr17 real parity fixture | 0 (3 default rows, 1001 dense rows) |
| rows changed in the chr20 dense real fixture | 0 (810 rows) |
| pre-existing registered tests exercising a firing record through GenotypeGVCFs | 0 |

## 1. Corpora located (task item 1)

| location | result |
| --- | --- |
| `gatk-source/src/test/resources/**/*.g.vcf` + `*.g.vcf.gz` | 284 matches, of which **142 are real files** and **20 are git-lfs pointers** (130–133 bytes, first bytes `version https://git-lfs...`; e.g. `large/NA12878.prod.chr20snippet.g.vcf.gz`, `large/gvcfs/HG00268.g.vcf.gz`, `large/testProductionGVCF.expected.g.vcf`, `large/gvcfs/CEUTrio.20.21.gatk3.4.g.vcf`) |
| the 1291-record chr20 HaplotypeCaller expected gVCF | `gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/haplotypecaller/expected.testGVCFMode.gatk4.g.vcf` — **measured 1291 data lines**, real file (153,088 bytes). A byte-identical copy exists at `tools/GenomicsDBImport/expected.testGVCFMode.gatk4.g.vcf` (1290 data lines, one fewer record). |
| `testdata/` | **no gVCFs at all** (`find testdata -name '*.g.vcf*'` → empty). It has BAMs (`testdata/real/ceutrio/CEUTrio.chr20.bam`, `CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam`, `testdata/real/dream_synthetic/*`), VCFs and references. |
| `fixtures/` | **no gVCFs**; `fixtures/chr20/` holds VCF outputs (`gatk.vcf.gz`, `native.vcf.gz`) for the HaplotypeCaller parity fixture. |
| runtime-generated real gVCFs (the actual GenotypeGVCFs gate inputs) | generated and probed in this round: chr17 `17:69000-70000` from `gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam` (**126 records**), chr20 dense region `20:10019901-10020710` from `fixtures/chr20/mnp.bam` (**49 records**), and a 100 kb real chr20 window from `testdata/real/ceutrio/CEUTrio.chr20.bam` (**1 record** — that window has no calls, which is the point of the registered `fastgatk-hc-chr20-100k-nocall-gatk-oracle` gate). |

Which corpora were scanned: all 142 real gVCFs listed in `.diag/round49/corpora.txt`
plus the three generated ones above. The 20 LFS pointers were **not** scanned (no data
in the checkout) — stated as a limit in §9.

`grep -rn 'testGVCFMode' --include=*.py --include=CMakeLists.txt fastgatk-native` finds
**no** registered gate that names these pinned resources: the project's GenotypeGVCFs
gates either build inline synthetic loci or generate their input gVCF at runtime. So the
"pinned GATK gVCF corpora" and "the project's real-data parity fixtures" are two distinct
sets here, and both are measured below.

## 2. The predicate, read from the source (task item 2)

Read directly (not from a paraphrase):

* `GATKVariantContextUtils.java:1443-1445` — `reverseTrimAlleles(vc)` = `trimAlleles(vc, false, true)`.
* `GATKVariantContextUtils.java:1454-1481` — the wrapper: guard at **`:1458`**
  (`getNAlleles() <= 1`, or any allele with `length() == 1 && !equals(Allele.SPAN_DEL)` → return input),
  candidate list = alleles that are neither symbolic nor `*` (`:1462-1463`),
  `normalizeAlleles(sequences, ranges, 0, true)` (`:1465`), `endTrim = shifts.right`,
  `startTrim = -shifts.left` (`:1466-1467`), `emptyAllele = any range.size() == 0` (`:1469`),
  `restoreOneBaseAtEnd = emptyAllele && startTrim == 0` (`:1470`),
  `endBasesToClip = restoreOneBaseAtEnd ? endTrim - 1 : endTrim` (`:1475`), and the recursion
  `trimAlleles(vc, (trimForward ? startBasesToClip : 0) - 1, trimReverse ? endBasesToClip : 0)` (`:1477`).
* `GATKVariantContextUtils.java:1489-1521` — the worker: `fwdTrimEnd == -1 && revTrim == 0`
  returns the input untouched (`:1492-1493`); symbolic and `*` are copied untouched (`:1497-1501`);
  otherwise `copyOfRange(a.getBases(), fwdTrimEnd+1, a.length()-revTrim)` (`:1504`);
  `start = inputVC.getStart() + (fwdTrimEnd + 1)` (`:1515`) and
  `stop = start + alleles.get(0).length() - 1` (`:1516`) — **POS never moves for a reverse-only trim**.
* `AlignmentUtils.java:818-845` — `normalizeAlleles`: the second and third loops cannot change the
  answer for `trimForward=false, maxShift=0` (the third needs `startShift < 0`, impossible); the
  first loop trims while `minSize > 0 && lastBaseOnRightIsSame(...)` (`:826-830`), so the clip is the
  common trailing run capped by the shortest candidate allele: `endShift == minSize` ⇒ clip
  `minSize - 1` (restore-one-base), otherwise clip `endShift`. Helpers at `:848-857`, `:859-868`, `:870-879`.

Net result (and the exact statement of the predicate): **a record is rewritten iff it has ≥2 alleles,
no allele of length 1 other than `*`, and its non-symbolic non-`*` alleles share ≥1 trailing base;
the clip is the shared trailing run capped at the shortest candidate, minus one when the cap was hit.**
The clip also has to be ≥ 1, otherwise `:1492` returns the input untouched.

`htsjdk` note: `Allele.isSymbolic()` and `length()` were checked against the class files in the pinned
fat jar (`unzip -l ... | grep htsjdk/variant/variantcontext/VariantContext.class` → present); its
`VariantContext.class` md5 (`86b0ef23513489a8781734f19dcbf0b1`) differs from the gradle-cached
`htsjdk-4.2.0.jar` (`700baf08b1fbbcb2a010fd8afe442381`), so I do not name an htsjdk version and treat
the jar as the oracle.

### The predicate used (verbatim, `.diag/round49/reverse_trim.py`)

```python
#!/usr/bin/env python3
"""Literal port of GATK 4.6.2.0's reverse allele trim, used for round-49 reach
measurement.

Sources read (not paraphrased):

* ``GATKVariantContextUtils.trimAlleles(inputVC, trimForward, trimReverse)``
  -- ``gatk-source/src/main/java/org/broadinstitute/hellbender/utils/variant/GATKVariantContextUtils.java:1454-1481``
  with ``reverseTrimAlleles`` at ``:1443-1445`` = ``trimAlleles(vc, false, true)``.
* ``GATKVariantContextUtils.trimAlleles(inputVC, fwdTrimEnd, revTrim)``
  -- same file ``:1489-1521`` (the ``copyOfRange(a.getBases(), fwdTrimEnd+1, a.length()-revTrim)``
  rewrite and the ``stop = start + alleles.get(0).length() - 1`` rebuild).
* ``AlignmentUtils.normalizeAlleles(sequences, bounds, maxShift, trim)``
  -- ``gatk-source/src/main/java/org/broadinstitute/hellbender/utils/read/AlignmentUtils.java:818-845``
  plus its helpers ``lastBaseOnRightIsSame`` ``:848-857``,
  ``firstBaseOnLeftIsSame`` ``:859-868``, ``nextBaseOnLeftIsSame`` ``:870-879``.

The port is deliberately *literal* (including the two loops that cannot run for
the reverse-only case, ``trimForward=false`` / ``maxShift=0``, and the
``Utils.nonEmpty`` precondition) so that it can be differentially tested against
the real method in the pinned jar rather than trusted.
"""
from __future__ import annotations

SPAN_DEL = "*"


def is_symbolic(allele: str) -> bool:
    """htsjdk ``Allele.isSymbolic()``: symbolic alleles start with '<'."""
    return allele.startswith("<")


def is_span_del(allele: str) -> bool:
    return allele == SPAN_DEL


# ---------------------------------------------------------------------------
# AlignmentUtils.normalizeAlleles (:818-845) and its three predicates
# ---------------------------------------------------------------------------

def _last_base_on_right_is_same(sequences, bounds) -> bool:
    last = sequences[0][bounds[0][1] - 1]
    return all(sequences[n][bounds[n][1] - 1] == last for n in range(len(sequences)))


def _first_base_on_left_is_same(sequences, bounds) -> bool:
    first = sequences[0][bounds[0][0]]
    return all(sequences[n][bounds[n][0]] == first for n in range(len(sequences)))


def _next_base_on_left_is_same(sequences, bounds) -> bool:
    left = sequences[0][bounds[0][0] - 1]
    return all(sequences[n][bounds[n][0] - 1] == left for n in range(len(sequences)))


def normalize_alleles(sequences, bounds, max_shift, trim):
    """Returns ``(startShift, endShift)``; ``bounds`` is mutated in place."""
    assert sequences
    assert len(sequences) == len(bounds)
    assert all(max_shift <= b[0] for b in bounds)

    start_shift = 0
    end_shift = 0
    min_size = min(e - s for s, e in bounds)

    # consume any redundant shared bases at the end of the alleles
    while trim and min_size > 0 and _last_base_on_right_is_same(sequences, bounds):
        for bound in bounds:
            bound[1] -= 1
        min_size -= 1
        end_shift += 1

    while trim and min_size > 0 and _first_base_on_left_is_same(sequences, bounds):
        for bound in bounds:
            bound[0] += 1
        min_size -= 1
        start_shift -= 1

    while (start_shift < max_shift and _next_base_on_left_is_same(sequences, bounds)
           and _last_base_on_right_is_same(sequences, bounds)):
        for bound in bounds:
            bound[0] -= 1
            bound[1] -= 1
        start_shift += 1
        end_shift += 1

    return start_shift, end_shift


# ---------------------------------------------------------------------------
# GATKVariantContextUtils.trimAlleles(:1454-1481) + trimAlleles(:1489-1521)
# ---------------------------------------------------------------------------

def reverse_trim_alleles(ref: str, alts, trim_forward: bool = False,
                         trim_reverse: bool = True):
    """Port of ``GATKVariantContextUtils.reverseTrimAlleles(record)``.

    Returns ``None`` when GATK returns the input record unchanged (either the
    ``:1458`` guard, or ``fwdTrimEnd == -1 && revTrim == 0`` at ``:1492``);
    otherwise returns the new allele list ``[REF, ALT...]``.  The record start is
    never moved (:1515) and the new stop is ``start + len(new REF) - 1`` (:1516).
    """
    alleles = [ref] + list(alts)

    # :1458 -- the guard
    if len(alleles) <= 1:
        return None
    if any(len(a) == 1 and not is_span_del(a) for a in alleles):
        return None

    candidates = [a for a in alleles if not is_symbolic(a) and not is_span_del(a)]
    if not candidates:  # Utils.nonEmpty(sequences) at AlignmentUtils.java:819
        raise ValueError("no non-symbolic, non-'*' allele: GATK would throw")
    sequences = [a.encode("ascii") for a in candidates]
    bounds = [[0, len(s)] for s in sequences]

    shifts = normalize_alleles(sequences, bounds, 0, True)
    end_trim = shifts[1]
    start_trim = -shifts[0]

    empty_allele = any(e - s == 0 for s, e in bounds)
    restore_one_base_at_end = empty_allele and start_trim == 0
    restore_one_base_at_start = empty_allele and start_trim > 0

    end_bases_to_clip = end_trim - 1 if restore_one_base_at_end else end_trim
    start_bases_to_clip = start_trim - 1 if restore_one_base_at_start else start_trim

    # trimAlleles(inputVC, (trimForward ? startBasesToClip : 0) - 1,
    #             trimReverse ? endBasesToClip : 0)
    fwd_trim_end = (start_bases_to_clip if trim_forward else 0) - 1
    rev_trim = end_bases_to_clip if trim_reverse else 0
    if fwd_trim_end == -1 and rev_trim == 0:
        return None

    if rev_trim < 0 or fwd_trim_end + 1 < 0:
        raise ValueError("negative clip: GATK's copyOfRange would throw")

    trimmed = []
    for allele in alleles:
        if is_symbolic(allele) or is_span_del(allele):
            trimmed.append(allele)
        else:
            trimmed.append(allele[fwd_trim_end + 1:len(allele) - rev_trim])
    return trimmed


def fires(ref: str, alts) -> bool:
    """True when GATK's reverse trim rewrites the record."""
    try:
        return reverse_trim_alleles(ref, alts) is not None
    except ValueError:
        return False


def clip_length(ref: str, alts):
    """The number of trailing bases GATK clips (0 when the record is untouched)."""
    result = reverse_trim_alleles(ref, alts)
    if result is None:
        return 0
    return len(ref) - len(result[0])
```

### Per-corpus counts

Aggregate over the 142 real gVCFs (`.diag/round49/per_corpus_table.tsv`,
`.diag/round49/reach_scan.out`):

```
corpora scanned: 142
records: 165527
  single-allele (getNAlleles() <= 1)          0
  one-base-allele-guard (:1458)               165524
  no-shared-trailing-base (endShift == 0)     3
  FIRES                                       0
  counterfactual fires without the :1458 guard: 0
  FIRES among records carrying a concrete ALT:  0
```

The 165,524 guard-blocked records split exactly as
**162,712 blocked by a one-base REF** (98.3 % of all records; 19,929 of these also carry a
one-base ALT, i.e. the ordinary SNP shape) and **2,812 blocked by a one-base ALT only**
(per-corpus split in the table below; e.g. `GnarlyGenotyper/NA12892.chrX.diploid.rb.g.vcf`
31,171 REF-blocked vs 1,328 ALT-blocked).

The three records that do reach `normalizeAlleles()` — every one of them, in every corpus — are:

```
tools/GenomicsDBImport/mnp.input.g.vcf   20:69772  TTA > CTC,<NON_REF>   -> unchanged (A != C)
tools/walkers/CombineGVCFs/mnp.g.vcf     20:69772  TTA > CTC,<NON_REF>   -> unchanged
tools/walkers/GenotypeGVCFs/mnp.input.g.vcf 20:69772 TTA > CTC,<NON_REF> -> unchanged
```

Independent confirmation from the real method: a Java probe over the same files
(`.diag/round49/RevTrimProbe.java`, calling `GATKVariantContextUtils.reverseTrimAlleles` from the
pinned jar) reports `total=165527 changed=0 changedVariant=0 changedNonVariant=0
changedWithSpanDel=0 errors=0` (`.diag/round49/probe-all.tsv`). No record made the real method throw.

### 5 concrete `before -> after` pairs

The corpora contain **no** firing record, so no example can come from them. The five below are the
gate's own measured fixtures plus measured stress-corpus rows; the "after" column is what GATK
actually emitted/returned, not a prediction:

| # | corpus | before (REF > ALTs) | after (REF > ALTs) | clip | source of "after" |
| --- | --- | --- | --- | --- | --- |
| 1 | gate fixture `suffix-substitution-default` | `AAAA > AACA,<NON_REF>` | `AAA > AAC` | 1 | measured GATK row, `.diag/round49/gate-revtrim-run1.log` |
| 2 | gate fixture `emptiness-guard-clips-one-less` | `ACGTACGT > ACGT,<NON_REF>` | `ACGTA > A` | 3 | measured GATK row (run capped at 4 → restore-one-base) |
| 3 | gate fixture `longer-alleles-same-one-base-run` | `AAAAA > AAACA,<NON_REF>` | `AAAA > AAAC` | 1 | measured GATK row |
| 4 | gate fixture `single-concrete-candidate-clipped-to-one-base` (dense) | `AA > *,<NON_REF>` | `A > *` | 1 | measured GATK row (`*` copied untouched) |
| 5 | stress corpus (randomized), real method | `GAAGACGTGCCT > *` | `G > *` | 11 | `.diag/round49/probe-stress.tsv` |

Measured clip-length histogram over the 2,990 firing stress records:
`{1: 1510, 2: 773, 3: 347, 4: 153, 5: 80, 6: 1, 7: 65, 11: 61}` — the clip really is the run
length, not a fixed one base.

## 3. Variant records vs reference blocks (task item 3)

**0** of the firing records are variant records, because **0** records fire at all. Both oracles agree:

* python: `FIRES among records carrying a concrete ALT: 0`; the population of records that carry at
  least one concrete (non-`*`, non-symbolic) ALT is **25,305**.
* real method: `changedVariant=0` and `changedNonVariant=0` on all 142 corpora.

One measurement worth recording because it breaks the obvious shortcut: in the pinned jar,
**`VariantContext.isVariant()` is `true` even for a `<NON_REF>`-only reference block**
(`getType()` = `SYMBOLIC`, measured on 1,034 such records of the 1291-record HC gVCF). It is *not*
usable as the reference-block discriminator at `GenotypeGVCFsEngine.java:153`; the engine's real
discriminator is `isProperlyPolymorphic(regenotypedVC) || includeNonVariants` at `:161`. This report
therefore counts "variant record" by **allele content** (any ALT that is not symbolic and not `*`),
and reports it explicitly rather than leaning on `isVariant()`.

## 4. Cross-checks against GATK (task item 4)

### 4a. The registered gate, run in full

`python3 fastgatk-native/scripts/verify_genotype_gvcf_reverse_trim_gatk_oracle.py`
→ `{"status": "pass", "violations": [], "gatk_version": "4.6.2.0", "mode": "strict"}`,
all 7 cases with `gatk_exit=0 native_exit=0` (`.diag/round49/gate-revtrim-run1.log`, `EXIT=0`).
That is 7 pinned GATK fixtures measured on this machine, not synthetic values.

### 4b. Does the predicate predict those measured rows?

`.diag/round49/gate_crosscheck.py` applies the predicate to each case's input allele list and
compares with the GATK rows the gate run just produced (concrete/`*` alleles compared; GenotypeGVCFs'
output-allele subset drops `<NON_REF>`, which is a different step — the trim copies symbolic alleles
untouched at `:1497-1501`):

| case | measured GATK (REF > ALTs) | predicate | verdict |
| --- | --- | --- | --- |
| suffix-substitution-default | `AAA > AAC` | `AAA > AAC` | MATCH |
| no-shared-trailing-base-unchanged | `AAAA > AACC` | `AAAA > AACC` (no clip) | MATCH |
| one-base-alt-guard | `AAA > A` | `AAA > A` (guard) | MATCH |
| emptiness-guard-clips-one-less | `ACGTA > A` | `ACGTA > A` | MATCH |
| longer-alleles-same-one-base-run | `AAAA > AAAC` | `AAAA > AAAC` | MATCH |
| trim-runs-after-output-allele-subset | `AAA > AAC` | raw `AAAA > AACA` (DIFFER) / after modelling the one-base-ALT prune `AAA > AAC` (MATCH) | MATCH with subset |
| single-concrete-candidate-clipped-to-one-base (dense) | `AAA > A`, `A > *` | same | MATCH |

The one case that needs the subset modelled is the gate's own property 3: the second sample's one-base
ALT is pruned by the standard-confidence subset *before* the trim, which is exactly what the recorded
GATK row shows. So the ordering claim is confirmed by measurement, and for the six other cases the raw
predicate alone reproduces GATK byte for byte. **No case disagreed.**

### 4c. Predicate vs the real GATK method on 5,000 randomized records

Because the corpora are all negatives, the "fires" branch was falsified separately:
`.diag/round49/differential_test.py` generated 5,000 valid VCF records (REF/ALT lengths 1–12, shared
trailing runs, `*`, `<NON_REF>`, multi-allelic, plus all seven gate fixtures), ran the real
`reverseTrimAlleles` from the pinned jar over them, and compared verdict **and the exact trimmed
REF/ALT strings** record by record:

```
records generated=5000 measured=5000
predicate fires=2990  GATK changed rows=2990
clip-length histogram (measured): [(1,1510),(2,773),(3,347),(4,153),(5,80),(6,1),(7,65),(11,61)]
disagreements=0
```

So the predicate is validated on 2,990 positives and 165,527 + 2,010 negatives.

### 4d. Idempotence (used again in §6)

Of the 2,990 GATK-trimmed stress records, **0** satisfy the fires predicate again, i.e. GATK's emitted
form never carries a shared trailing run. This is what makes "GATK's output shows no trace of the
trim, so the trim cannot have changed this fixture input" a measurement rather than a guess.

## 5. Could the pre-existing registered corpus have detected it? (task item 5)

Method: text scan of every gate script and C++ test for embedded VCF data lines (the project's
fixtures are overwhelmingly inline string literals), extracting REF/ALT and applying the predicate;
hits mapped to CTest names from `fastgatk-native/CMakeLists.txt`
(`.diag/round49/registered_corpus_scan.py`). 354 files scanned (338 gate scripts including `.sh`,
6 C++ tests under `fastgatk-native/tests`, plus the `tests/fixtures` entries). A second, broader net
(any `REF<TAB>ALT`-shaped allele pair, no VCF-structure requirement) was run over the 43 registered
scripts whose test names contain `genotype`/`gvcf`; it found **no additional** firing record.

**Records the trim rewrites (REF ≥ 2, all concrete ALTs ≥ 2, shared trailing base) — 7 hits in 4 files:**

| file:line | record | predicted after | registered as |
| --- | --- | --- | --- |
| `fastgatk-native/scripts/verify_genotype_gvcf_reverse_trim_gatk_oracle.py:122` | `chr1:2 AAAA > AACA,<NON_REF>` | `AAA > AAC` | `fastgatk-genotype-gvcf-reverse-trim-gatk-oracle` — **the new gate of this round** |
| `…/verify_genotype_gvcf_reverse_trim_gatk_oracle.py:141` | `chr1:2 ACGTACGT > ACGT,<NON_REF>` | `ACGTA > A` | same (new) |
| `…/verify_genotype_gvcf_reverse_trim_gatk_oracle.py:146` | `chr1:2 AAAAA > AAACA,<NON_REF>` | `AAAA > AAAC` | same (new) |
| `…/verify_genotype_gvcf_reverse_trim_gatk_oracle.py:153` | `chr1:2 AAAA > AACA,<NON_REF>` (2 samples) | `AAA > AAC` | same (new) |
| `fastgatk-native/scripts/verify_genotype_gvcf_dense_materialize_gatk_oracle.py:152` | `chr1:2 AAAA > AACA,<NON_REF>` | `AAA > AAC` | **(not registered)** — added in the 2026-09-11 dense round; its case `suffix-substitution-default-reverse-trim` (`:338-350`) exists specifically to pin this divergence and its `why` text says "Native emits 'AAAA/AACA' untrimmed" |
| `fastgatk-native/scripts/verify_left_align.py:52` | `chr1:7 ATG > ACG` | `AT > AC` | `fastgatk-left-align-trim-contract` — **LeftAlignAndTrimVariants**, input, not GenotypeGVCFs |
| `fastgatk-native/scripts/verify_reblock_gvcf.py:422` | `chr1:50 AC > GC,TC,<NON_REF>` | `A > G,T,<NON_REF>` | `fastgatk-reblock-gvcf-contract` — **ReblockGVCF**, input, not GenotypeGVCFs |

`grep -rn reverseTrimAlleles gatk-source/src/main/java` returns exactly three call sites:
`GenotypeGVCFsEngine.java:167`, `HaplotypeCallerGenotypingEngine.java:549` and
`VCFComparator.java:864`. Neither LeftAlignAndTrimVariants nor ReblockGVCF calls it, so the two
pre-existing hits cannot observe this divergence.

**Near misses (guard-free but no shared trailing base, so no clip) — 26 hits**, including GenotypeGVCFs
gates that come closest without firing:

* `verify_genotype_gvcf_spandel_gatk_oracle.py:818, 828` `chr1:2 AAAA > AAC,<NON_REF>`;
  `:848` `chr1:2 AA > *,AC,<NON_REF>`; `:992` `AAA > AAC`; `:1005` `AA > *,AC` — all `shared_run=0`.
* `verify_hc_chr20_max_mnp_gvcf_gatk_oracle.py:79` `20:10020680 CA > AT`; `:97` `CA > AT,CG` — `shared_run=0`.
* `verify_genotype_gvcf_spanning_deletion_gatk_oracle.py:48,57` `17:69067 T > C,*,<NON_REF>` (one-base REF ⇒ guarded).

**Conclusion:** before this fix, **no registered test exercised a record that the reverse trim
rewrites**, and none could have: the registered GenotypeGVCFs inputs are either inline
single-base-REF synthetic loci or HaplotypeCaller-generated gVCFs from the chr17/chr20 fixtures, and
both generated inputs contain **0** guard-passing records out of 126 and 49 records respectively
(measured in §6). The divergence was invisible to the corpus of that time. The first artifact in the
tree that names it is the *unregistered* dense diagnostic gate from the immediately preceding round.

## 6. Practical significance: does it change any row of the real-data fixtures? (task item 6)

**No.** Three independent lines of evidence.

**(a) The pre-fix run of the real-data parity gate was green.** The repository's own log
`.diag/emitted-ownership-ctest-genotype-gvcf.log` has mtime **2026-09-12T20:45**, which precedes every
artifact of the fix — `genotype_gvcf_tool.cpp` 23:40, the gate script 23:38, `CMakeLists.txt` 23:41,
the binary `fastgatk-native/build/fastgatk-genotype-gvcf` 23:40 (and the fix is an uncommitted
working-tree change on top of `HEAD` = `f0a8277`, 21:52, whose own message records that the reverse
divergence was found only then). That run executes 18 `genotype-gvcf` tests, contains **no**
`reverse-trim` test, and
reports `fastgatk-genotype-gvcf-gatk-oracle ... Passed` with
`100% tests passed, 0 tests failed out of 18` (the same is true of
`.diag/cs-ctest-genotype-gvcf.log`, 12:50, and `.diag/filter-pass-ctest-genotype-gvcf.log`). That gate
(`verify_gatk_genotype_gvcf.py`) keys its comparison on `(chrom, pos, REF, ALT)`
(`record_map()`) plus full FORMAT text, so a trimmed-vs-untrimmed REF/ALT would have shown up as a
missing/extra record and failed the test. It passed pre-fix ⇒ the divergence changed no row there.

**(b) Every emitted row of both real-data fixtures has a one-base REF.** Measured on GATK 4.6.2.0 output
I generated in this round (`.diag/round49/parity/`, `.diag/round49/dense/`):

| fixture | input gVCF | emitted rows | rows with REF ≥ 2 bp | rows carrying `*` |
| --- | --- | --- | --- | --- |
| `fastgatk-genotype-gvcf-gatk-oracle` (chr17 `17:69000-70000`, NA12878 chr17s 69k–70k BAM) default | 126 records, 3 variant records, **0 guard-passing** | 3 (`17:69067 T>G`, `17:69368 G>C`, `17:69631 C>T`) | 0 | 0 |
| same, `--include-non-variant-sites` | — | 1001 | 0 | 0 |
| `fastgatk-hc-dense-gvcf-genotype-gatk-oracle` (chr20 `20:10019901-10020710`, `fixtures/chr20/mnp.bam`) dense | 49 records, 11 variant records, **0 guard-passing** | 810 | 0 | 0 |

**(c) Why that is sufficient.** The trim runs on the *merged* record at
`GenotypeGVCFsEngine.java:167`. `ReferenceConfidenceVariantContextMerger.java:137` +
`determineReferenceAlleleGivenReferenceBase` (`:417-426`) take the merged REF from the contributing
records' REF at a record-start locus, and use a **single** reference base (`Allele.create(refBase, true)`)
when every contributing record is spanning (the synthesized dense rows); the later output-allele subset
can only *drop* alleles. Both fixture inputs are 100 % one-base-REF records, so **every** pre-trim
merged allele list in both fixtures contains a length-1 non-`*` allele ⇒ `:1458` short-circuits ⇒ the
trim is a no-op at every locus in both traversal modes. Reasoning, not instrumentation: I cannot call
the pre-trim merge directly, and (a)+(b) are the measured backstops.

**(d) Multi-sample merging**, the one shape where a merged allele list can differ from every single
input record, was checked separately: across all 142 corpora only **7** extra records share a
`(contig, pos, REF)` with an earlier record (2 corpora: `CombineGVCFs/sample2.MT.g.vcf` 3,
`ValidateVariants/badGVCF.outOfOrder.g.vcf` 4), and the union-of-alleles view fires **0** times.
Derived (union, not the merger's dedup), and it does not change the answer.

## 7. Artifacts

* Report: `fastgatk-native/evidence/2026-09-13-round49/reverse-trim-reach-on-real-corpora.md` (this file).
* Scratch (git-ignored), all under `.diag/round49/`:
  `reverse_trim.py` (the port), `RevTrimProbe.java` + `classes/` (real-method probe),
  `corpora.txt`, `reach_scan.py` + `reach_scan.out`, `per_corpus_table.tsv`, `per_corpus_compact.txt`,
  `differential_test.py`, `probe-all.tsv`, `probe-stress.tsv`, `probe-sanity.tsv`,
  `gate_crosscheck.py`, `registered_corpus_scan.py`, `gate-revtrim-run1.log`,
  `parity/` (chr17 gVCF + GATK joint VCFs), `dense/`, `ceutrio/`.

Reproduction (each is a single command; nothing was rebuilt):

```
python3 fastgatk-native/scripts/verify_genotype_gvcf_reverse_trim_gatk_oracle.py     # 4a: 7/7 pass
python3 .diag/round49/reach_scan.py                                                 # 2/3: per-corpus counts
python3 .diag/round49/differential_test.py                                          # 4c: 5000 records, 0 disagreements
python3 .diag/round49/gate_crosscheck.py                                            # 4b + 4d
python3 .diag/round49/registered_corpus_scan.py                                     # 5
```

## 8. Which claims are measured, which are reasoned

**Measured (instrumented, value recorded):**
per-corpus record counts and guard/fire classification (python, 165,527 records); the same files pushed
through the pinned jar's own `reverseTrimAlleles` (0 changed, 0 errors); the 7 gate cases run against
pinned GATK (`status=pass`, `violations=[]`); the differential test (5,000 records, 2,990 positives,
0 disagreements, exact trimmed allele strings); idempotence (2,990 trimmed rows, 0 re-fires);
the two generated real fixture gVCFs (126 and 49 records, 0 guard-passing) and their GATK
GenotypeGVCFs outputs (3 / 1001 / 810 rows, 0 rows with REF ≥ 2 bp, 0 `*` rows); the 1291-record
file's record count; LFS-pointer detection; the pre-fix green CTest logs; the textual fixture scan
(354 files, 7 hits, 26 near misses).

**Reasoned, not instrumented:**
(i) that the merged pre-trim allele list at a fixture locus equals the input record's allele list —
read from `ReferenceConfidenceVariantContextMerger.java:137,417-426` and backstopped by §6(a);
(ii) that a pre-fix native binary would have emitted the untrimmed merge (read from the working-tree
diff: the trim calls at `genotype_gvcf_tool.cpp:3122,3159,3389` are the entire change), since no
pre-fix binary exists to run;
(iii) that HC's own output is already reverse-trimmed, which explains the zero reach — read from
`HaplotypeCallerGenotypingEngine.java:546-549`, not measured.

## 9. What I could not falsify / limits

* **The pre-fix native binary was never executed** (rebuilding is forbidden and the binary is current),
  so "pre-fix output = GATK output on these fixtures" is reasoned from the working-tree diff plus the
  measured absence of any trimmable locus, and corroborated by the pre-fix green gate log.
* **20 git-lfs gVCFs were not readable** in this checkout, including
  `large/NA12878.prod.chr20snippet.g.vcf.gz` and `large/testProductionGVCF.expected.g.vcf`. If any
  larger production corpus exists, it is not in the tree; the pinned-corpus result (0 of 165,527)
  cannot be extended to it by measurement.
* **The corpus scan is textual.** Fixtures assembled programmatically (e.g. allele lists built from
  computed values) could hide a firing record from the regex; the broad second-pass net bounded that
  risk for the registered GenotypeGVCFs scripts, and the two real-data gates were verified by probing
  the actual generated gVCFs rather than by reading their source.
* **Dense-mode synthesized rows** are judged from the merger source plus the measured 1-base REF of
  GATK's own 1001/810 dense rows; the pre-trim content of a synthesized row was not observed directly.
  The previous round's evidence (`2026-09-11-wave0/round-dense-materialize.md` §1.4) measured the same
  one-base-REF property on synthetic dense fixtures.
* **The empty-allele corner of the port**: my literal port raises if an allele list contains an empty
  concrete allele (VCF cannot express one); the harness hit this only through a greedy-regex artifact
  in the coverage audit, never through a real record. GATK's `copyOfRange` would throw on such an
  input too, so this does not affect any count here.
* **`fires_without_guard()` counterfactual**: it reuses `normalizeAlleles` on the candidate list; it is
  a hypothetical (the guard is real GATK behaviour) and is reported only to show that even removing the
  guard would not reach these corpora.

## 10. Per-corpus table (all 142 real gVCFs)

`guard-passing` = records that survive the `:1458` guard and reach `normalizeAlleles()`;
`FIRES` = records GATK actually rewrites. Only three rows have `guard-passing > 0`, none has `FIRES > 0`.

| # | corpus (under gatk-source/src/test/resources/) | records | concrete-ALT records | guard-passing | FIRES |
| --- | --- | --- | --- | --- | --- |
| 1 | `engine/GenomicsDBIntegration/intervalsRestrictedExpected.g.vcf` | 3 | 1 | 0 | 0 |
| 2 | `engine/GenomicsDBIntegration/tiny.g.vcf` | 11 | 3 | 0 | 0 |
| 3 | `tools/ArtificalPhasedData.1.g.vcf` | 12 | 12 | 0 | 0 |
| 4 | `tools/GenomicsDBImport/expected.testGVCFMode.gatk4.g.vcf` | 1290 | 257 | 0 | 0 |
| 5 | `tools/GenomicsDBImport/expected.testGVCFMode.gatk4.posteriors.g.vcf` | 815 | 260 | 0 | 0 |
| 6 | `tools/GenomicsDBImport/mnp.input.g.vcf` | 1 | 1 | 1 | 0 |
| 7 | `tools/GenomicsDBImport/testHeaderContigLineSorting1.g.vcf` | 1 | 0 | 0 | 0 |
| 8 | `tools/GenomicsDBImport/testHeaderContigLineSorting2.g.vcf` | 1 | 0 | 0 | 0 |
| 9 | `tools/IndexFeatureFile/test_variants_for_index.g.vcf` | 1 | 0 | 0 | 0 |
| 10 | `tools/haplotypecaller/expected.CEUTrio.HiSeq.WGS.b37.NA12878.CONTAMINATED.WITH.HCC1143.NORMALS.15PCT.20.10100000-10150000.gatk3.8-1-1-gdde23f56a6.g.vcf` | 579 | 121 | 0 | 0 |
| 11 | `tools/haplotypecaller/expected.CEUTrio.HiSeq.WGS.b37.NA12878.CONTAMINATED.WITH.HCC1143.NORMALS.15PCT.20.10100000-10150000.postIndelRefConfUpdate.g.vcf` | 477 | 121 | 0 | 0 |
| 12 | `tools/haplotypecaller/expected.CEUTrio.HiSeq.WGS.b37.NA12878.calls.20.10100000-10150000.g.vcf` | 529 | 116 | 0 | 0 |
| 13 | `tools/haplotypecaller/expected.testGVCFMode.3.8-4-g7b0250253f.g.vcf` | 1556 | 260 | 0 | 0 |
| 14 | `tools/haplotypecaller/expected.testGVCFMode.gatk3.5.alleleSpecific.g.vcf` | 1572 | 260 | 0 | 0 |
| 15 | `tools/haplotypecaller/expected.testGVCFMode.gatk3.5.g.vcf` | 1572 | 260 | 0 | 0 |
| 16 | `tools/haplotypecaller/expected.testGVCFMode.gatk3.8-4-g7b0250253f.alleleSpecific.g.vcf` | 1556 | 260 | 0 | 0 |
| 17 | `tools/haplotypecaller/expected.testGVCFMode.gatk4.alleleSpecific.g.vcf` | 1291 | 257 | 0 | 0 |
| 18 | `tools/haplotypecaller/expected.testGVCFMode.gatk4.g.vcf` | 1291 | 257 | 0 | 0 |
| 19 | `tools/haplotypecaller/expected.testGVCFMode.gatk4.posteriors.g.vcf` | 815 | 260 | 0 | 0 |
| 20 | `tools/mitochondria/NA12878.MT.filtered.g.vcf` | 151 | 37 | 0 | 0 |
| 21 | `tools/mitochondria/NA12878.MT.g.vcf` | 6590 | 310 | 0 | 0 |
| 22 | `tools/walkers/CombineGVCFs/NA12878.AS.NON_REF_remap_check.chr20snippet.g.vcf` | 4 | 2 | 0 | 0 |
| 23 | `tools/walkers/CombineGVCFs/NA12878.AS.chr20snippet.g.vcf` | 208 | 26 | 0 | 0 |
| 24 | `tools/walkers/CombineGVCFs/NA12878.MT.filtered.g.vcf` | 6590 | 310 | 0 | 0 |
| 25 | `tools/walkers/CombineGVCFs/NA12891.MT.filtered.g.vcf` | 6215 | 420 | 0 | 0 |
| 26 | `tools/walkers/CombineGVCFs/NA12892.AS.chr20snippet.g.vcf` | 170 | 7 | 0 | 0 |
| 27 | `tools/walkers/CombineGVCFs/NA12892.test.rb.g.vcf` | 50 | 6 | 0 | 0 |
| 28 | `tools/walkers/CombineGVCFs/NA19240.MT.filtered.g.vcf` | 6316 | 516 | 0 | 0 |
| 29 | `tools/walkers/CombineGVCFs/YRIoffspring.chr20snippet.g.vcf` | 10681 | 315 | 0 | 0 |
| 30 | `tools/walkers/CombineGVCFs/gvcfExample1WithTrailingReferenceBlocks.g.vcf` | 13 | 3 | 0 | 0 |
| 31 | `tools/walkers/CombineGVCFs/gvcfExample2WithTrailingReferenceBlocks.g.vcf` | 7 | 0 | 0 | 0 |
| 32 | `tools/walkers/CombineGVCFs/gvcfWithTrailingReferenceBlocksBandedExpected.g.vcf` | 43 | 3 | 0 | 0 |
| 33 | `tools/walkers/CombineGVCFs/gvcfWithTrailingReferenceBlocksExpected.g.vcf` | 20 | 3 | 0 | 0 |
| 34 | `tools/walkers/CombineGVCFs/mnp.g.vcf` | 1 | 1 | 1 | 0 |
| 35 | `tools/walkers/CombineGVCFs/newMQcalc.combined.g.vcf` | 11514 | 371 | 0 | 0 |
| 36 | `tools/walkers/CombineGVCFs/sample1.MT.g.vcf` | 114 | 24 | 0 | 0 |
| 37 | `tools/walkers/CombineGVCFs/sample1.expected.MT.g.vcf` | 114 | 24 | 0 | 0 |
| 38 | `tools/walkers/CombineGVCFs/sample2.MT.g.vcf` | 110 | 23 | 0 | 0 |
| 39 | `tools/walkers/CombineGVCFs/spanningDel.1.g.vcf` | 3 | 1 | 0 | 0 |
| 40 | `tools/walkers/CombineGVCFs/spanningDel.2.g.vcf` | 3 | 1 | 0 | 0 |
| 41 | `tools/walkers/CombineGVCFs/spanningDel.many.g.vcf` | 3 | 1 | 0 | 0 |
| 42 | `tools/walkers/CombineGVCFs/spanningDel.many.haploid.g.vcf` | 3 | 1 | 0 | 0 |
| 43 | `tools/walkers/CombineGVCFs/spanningDel.many.tetraploid.g.vcf` | 3 | 1 | 0 | 0 |
| 44 | `tools/walkers/CombineGVCFs/spanningDeletionBaseExtensionTestExpected.g.vcf` | 14 | 3 | 0 | 0 |
| 45 | `tools/walkers/CombineGVCFs/tetraploidRun.GATK3.g.vcf` | 1839 | 8 | 0 | 0 |
| 46 | `tools/walkers/CombineGVCFs/twoSamples.MT.g.vcf` | 9633 | 761 | 0 | 0 |
| 47 | `tools/walkers/GenomicsDBImport/newMQcalc.combined.g.vcf` | 11512 | 371 | 0 | 0 |
| 48 | `tools/walkers/GenotypeGVCFs/CEUTrio.20.21.missingIndel.g.vcf` | 1 | 1 | 0 | 0 |
| 49 | `tools/walkers/GenotypeGVCFs/NA12878.AS.chr20snippet.g.vcf` | 208 | 26 | 0 | 0 |
| 50 | `tools/walkers/GenotypeGVCFs/anotherError.g.vcf` | 1 | 1 | 0 | 0 |
| 51 | `tools/walkers/GenotypeGVCFs/badHet.g.vcf` | 1 | 1 | 0 | 0 |
| 52 | `tools/walkers/GenotypeGVCFs/chr21.bad.pl.g.vcf` | 17 | 1 | 0 | 0 |
| 53 | `tools/walkers/GenotypeGVCFs/combineReblocked.g.vcf` | 2 | 2 | 0 | 0 |
| 54 | `tools/walkers/GenotypeGVCFs/combined.MT.g.vcf` | 124 | 24 | 0 | 0 |
| 55 | `tools/walkers/GenotypeGVCFs/compareWithoutPLs.g.vcf` | 2 | 2 | 0 | 0 |
| 56 | `tools/walkers/GenotypeGVCFs/homVar1DP.g.vcf` | 1 | 1 | 0 | 0 |
| 57 | `tools/walkers/GenotypeGVCFs/homVarNoDP.g.vcf` | 1 | 1 | 0 | 0 |
| 58 | `tools/walkers/GenotypeGVCFs/leadingDeletion.g.vcf` | 5 | 2 | 0 | 0 |
| 59 | `tools/walkers/GenotypeGVCFs/mixHaploidDiploidHighAlt/haploid.rb.g.vcf` | 2 | 2 | 0 | 0 |
| 60 | `tools/walkers/GenotypeGVCFs/mixHaploidDiploidHighAlt/s01.rb.g.vcf` | 2 | 2 | 0 | 0 |
| 61 | `tools/walkers/GenotypeGVCFs/mixHaploidDiploidHighAlt/s02.rb.g.vcf` | 2 | 2 | 0 | 0 |
| 62 | `tools/walkers/GenotypeGVCFs/mixHaploidDiploidHighAlt/s03.rb.g.vcf` | 2 | 2 | 0 | 0 |
| 63 | `tools/walkers/GenotypeGVCFs/mnp.input.g.vcf` | 1 | 1 | 1 | 0 |
| 64 | `tools/walkers/GenotypeGVCFs/multiSamples.expected.g.vcf` | 5 | 5 | 0 | 0 |
| 65 | `tools/walkers/GenotypeGVCFs/multiSamples.g.vcf` | 80 | 5 | 0 | 0 |
| 66 | `tools/walkers/GenotypeGVCFs/newMQcalc.combined.g.vcf` | 11620 | 373 | 0 | 0 |
| 67 | `tools/walkers/GenotypeGVCFs/noReads.g.vcf` | 1 | 1 | 0 | 0 |
| 68 | `tools/walkers/GenotypeGVCFs/spanningDel.combined.g.vcf` | 6 | 2 | 0 | 0 |
| 69 | `tools/walkers/GenotypeGVCFs/spanningDel.delOnly.g.vcf` | 1 | 1 | 0 | 0 |
| 70 | `tools/walkers/GenotypeGVCFs/spanningDel.depr.delOnly.g.vcf` | 1 | 1 | 0 | 0 |
| 71 | `tools/walkers/GenotypeGVCFs/test.tooManyAltsNoPLs.g.vcf` | 1 | 1 | 0 | 0 |
| 72 | `tools/walkers/GenotypeGVCFs/testAlleleSpecificAnnotations.CombineGVCF.expected.g.vcf` | 24 | 24 | 0 | 0 |
| 73 | `tools/walkers/GenotypeGVCFs/testAlleleSpecificAnnotations.CombineGVCF.output.g.vcf` | 349 | 26 | 0 | 0 |
| 74 | `tools/walkers/GenotypeGVCFs/threeSamples.MT.g.vcf` | 11583 | 1070 | 0 | 0 |
| 75 | `tools/walkers/GenotypeGVCFs/twoReblocked.g.vcf` | 139 | 24 | 0 | 0 |
| 76 | `tools/walkers/GenotypeGVCFs/withOxoGReadCounts.g.vcf` | 448 | 112 | 0 | 0 |
| 77 | `tools/walkers/GnarlyGenotyper/NA12891.chrX.haploid.rb.g.vcf` | 16509 | 6878 | 0 | 0 |
| 78 | `tools/walkers/GnarlyGenotyper/NA12892.chrX.diploid.rb.g.vcf` | 32499 | 9941 | 0 | 0 |
| 79 | `tools/walkers/GnarlyGenotyper/NA20846.rb.g.vcf` | 1 | 1 | 0 | 0 |
| 80 | `tools/walkers/GnarlyGenotyper/NA20890.rb.g.vcf` | 1 | 0 | 0 | 0 |
| 81 | `tools/walkers/GnarlyGenotyper/chrY_haploid_dragen.g.vcf` | 4 | 1 | 0 | 0 |
| 82 | `tools/walkers/GnarlyGenotyper/emptyASAnnotations.g.vcf` | 1 | 1 | 0 | 0 |
| 83 | `tools/walkers/GnarlyGenotyper/fake_sample2.rb.g.vcf` | 2 | 1 | 0 | 0 |
| 84 | `tools/walkers/GnarlyGenotyper/testNoReads.rb.g.vcf` | 1 | 0 | 0 | 0 |
| 85 | `tools/walkers/VCFComparator/actual.NA12878.rb.g.vcf` | 77 | 25 | 0 | 0 |
| 86 | `tools/walkers/VCFComparator/diploid.rb.g.vcf` | 2 | 2 | 0 | 0 |
| 87 | `tools/walkers/VCFComparator/expected.NA12878.rb.g.vcf` | 76 | 25 | 0 | 0 |
| 88 | `tools/walkers/VCFComparator/haploid.rb.g.vcf` | 2 | 2 | 0 | 0 |
| 89 | `tools/walkers/ValidateVariants/NA12891.AS.chr20snippet.BAD_MISSING_NON_REF.g.vcf` | 169 | 21 | 0 | 0 |
| 90 | `tools/walkers/ValidateVariants/NA12891.AS.chr20snippet.g.vcf` | 169 | 22 | 0 | 0 |
| 91 | `tools/walkers/ValidateVariants/NA12891.AS.chr20snippet.missingrefblock.g.vcf` | 168 | 21 | 0 | 0 |
| 92 | `tools/walkers/ValidateVariants/NA12891.AS.chr20snippet_BAD_INCOMPLETE_REGION.g.vcf` | 169 | 21 | 0 | 0 |
| 93 | `tools/walkers/ValidateVariants/NA12891.AS.fullchr20_with_chr20_dict_only.g.vcf` | 169 | 21 | 0 | 0 |
| 94 | `tools/walkers/ValidateVariants/badGVCF.outOfOrder.g.vcf` | 118 | 24 | 0 | 0 |
| 95 | `tools/walkers/ValidateVariants/badGVCF.outOfOrderThreeContigs.g.vcf` | 13 | 3 | 0 | 0 |
| 96 | `tools/walkers/ValidateVariants/badGVCF.outOfOrderTwoContigs.g.vcf` | 8 | 2 | 0 | 0 |
| 97 | `tools/walkers/ValidateVariants/goodGVCF.inOrderThreeContigs.g.vcf` | 11 | 3 | 0 | 0 |
| 98 | `tools/walkers/ValidateVariants/goodGVCF.inOrderTwoContigs.g.vcf` | 7 | 2 | 0 | 0 |
| 99 | `tools/walkers/ValidateVariants/hasRefBlockOverlaps.g.vcf` | 7 | 2 | 0 | 0 |
| 100 | `tools/walkers/ValidateVariants/hasVariantOverlappingRefBlock.g.vcf` | 7 | 2 | 0 | 0 |
| 101 | `tools/walkers/annotator/allelespecific/NA12878.AS.InsertSizeRankSum.chr20snippet.g.vcf` | 181 | 24 | 0 | 0 |
| 102 | `tools/walkers/annotator/allelespecific/NA12878.AS.MateRankSum.chr20snippet.g.vcf` | 181 | 24 | 0 | 0 |
| 103 | `tools/walkers/annotator/allelespecific/NA12878.AS.chr20snippet.g.vcf` | 208 | 26 | 0 | 0 |
| 104 | `tools/walkers/annotator/allelespecific/NA12891.AS.chr20snippet.g.vcf` | 169 | 22 | 0 | 0 |
| 105 | `tools/walkers/annotator/allelespecific/NA12892.AS.InsertSizeRankSum.chr20snippet.g.vcf` | 168 | 6 | 0 | 0 |
| 106 | `tools/walkers/annotator/allelespecific/NA12892.AS.MateRankSum.chr20snippet.g.vcf` | 168 | 6 | 0 | 0 |
| 107 | `tools/walkers/annotator/allelespecific/NA12892.AS.chr20snippet.g.vcf` | 170 | 7 | 0 | 0 |
| 108 | `tools/walkers/variantutils/ReblockGVCF/HG002.snippet.g.vcf` | 1 | 1 | 0 | 0 |
| 109 | `tools/walkers/variantutils/ReblockGVCF/chr20.shard0.g.vcf` | 144 | 0 | 0 | 0 |
| 110 | `tools/walkers/variantutils/ReblockGVCF/chr20.shard1.g.vcf` | 112 | 1 | 0 | 0 |
| 111 | `tools/walkers/variantutils/ReblockGVCF/chr20.shard2.g.vcf` | 209 | 0 | 0 | 0 |
| 112 | `tools/walkers/variantutils/ReblockGVCF/chr20.shard3.g.vcf` | 197 | 0 | 0 | 0 |
| 113 | `tools/walkers/variantutils/ReblockGVCF/dragen.g.vcf` | 5 | 3 | 0 | 0 |
| 114 | `tools/walkers/variantutils/ReblockGVCF/dropGQ0Dels.g.vcf` | 5 | 4 | 0 | 0 |
| 115 | `tools/walkers/variantutils/ReblockGVCF/expected.NA12878.AS.chr20snippet.reblocked.g.vcf` | 63 | 24 | 0 | 0 |
| 116 | `tools/walkers/variantutils/ReblockGVCF/expected.NA12878.AS.chr20snippet.reblocked.hiRes.g.vcf` | 77 | 24 | 0 | 0 |
| 117 | `tools/walkers/variantutils/ReblockGVCF/expected.NA12892.AS.chr20snippet.reblocked.g.vcf` | 18 | 6 | 0 | 0 |
| 118 | `tools/walkers/variantutils/ReblockGVCF/expected.aggressiveQualFiltering.g.vcf` | 8 | 2 | 0 | 0 |
| 119 | `tools/walkers/variantutils/ReblockGVCF/expected.overlappingDeletions.g.vcf` | 6 | 1 | 0 | 0 |
| 120 | `tools/walkers/variantutils/ReblockGVCF/gvcfForReblocking.g.vcf` | 12 | 4 | 0 | 0 |
| 121 | `tools/walkers/variantutils/ReblockGVCF/justHeader.g.vcf` | 0 | 0 | 0 | 0 |
| 122 | `tools/walkers/variantutils/ReblockGVCF/noCallGTs.g.vcf` | 201 | 12 | 0 | 0 |
| 123 | `tools/walkers/variantutils/ReblockGVCF/nonRefAD.g.vcf` | 2 | 2 | 0 | 0 |
| 124 | `tools/walkers/variantutils/ReblockGVCF/overlappingDeletions.hc.g.vcf` | 8 | 7 | 0 | 0 |
| 125 | `tools/walkers/variantutils/ReblockGVCF/prod.chr20snippet.withRawMQ.expected.g.vcf` | 35 | 4 | 0 | 0 |
| 126 | `tools/walkers/variantutils/ReblockGVCF/prod.chr20snippet.withRawMQ.g.vcf` | 1900 | 4 | 0 | 0 |
| 127 | `tools/walkers/variantutils/ReblockGVCF/prodWesInput.g.vcf` | 53 | 2 | 0 | 0 |
| 128 | `tools/walkers/variantutils/ReblockGVCF/prodWesOutput.g.vcf` | 23 | 1 | 0 | 0 |
| 129 | `tools/walkers/variantutils/ReblockGVCF/prodWgsInput.g.vcf` | 695 | 40 | 0 | 0 |
| 130 | `tools/walkers/variantutils/ReblockGVCF/prodWgsOutput.g.vcf` | 195 | 29 | 0 | 0 |
| 131 | `tools/walkers/variantutils/ReblockGVCF/reblock_cleanup_bug_variant.g.vcf` | 14 | 6 | 0 | 0 |
| 132 | `tools/walkers/variantutils/ReblockGVCF/testFirstPositionOnContigNotDropped.g.vcf` | 2 | 2 | 0 | 0 |
| 133 | `tools/walkers/variantutils/ReblockGVCF/testJustOneSample.expected.g.vcf` | 1 | 1 | 0 | 0 |
| 134 | `tools/walkers/variantutils/ReblockGVCF/testLeftContigBoundary.g.vcf` | 4 | 0 | 0 | 0 |
| 135 | `tools/walkers/variantutils/ReblockGVCF/testNonRefADCorrection.expected.g.vcf` | 2 | 1 | 0 | 0 |
| 136 | `tools/walkers/variantutils/ReblockGVCF/testOneSampleAsForGnomAD.expected.g.vcf` | 8 | 2 | 0 | 0 |
| 137 | `tools/walkers/variantutils/ReblockGVCF/treeScoreGvcf.g.vcf` | 12 | 4 | 0 | 0 |
| 138 | `tools/walkers/variantutils/SelectVariants/diploid-multisample-sac.g.vcf` | 56 | 3 | 0 | 0 |
| 139 | `tools/walkers/variantutils/SelectVariants/gvcfExample.g.vcf` | 12 | 4 | 0 | 0 |
| 140 | `tools/walkers/variantutils/SelectVariants/tetraploid-multisample-sac.g.vcf` | 80 | 3 | 0 | 0 |
| 141 | `utils/IndexUtils/test_variants_for_index.g.vcf` | 1 | 0 | 0 | 0 |
| 142 | `utils/variant/writers/small.g.vcf` | 4 | 1 | 0 | 0 |

