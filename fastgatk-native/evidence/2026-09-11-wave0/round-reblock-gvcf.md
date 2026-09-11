# Round report — ReblockGVCF `--drop-low-quals` (CASE A) and triploid reblocking (CASE B)

Scope: two pinned divergences in `fastgatk-native/scripts/verify_reblock_gvcf.py`
(lines 64/82-83 and 283-287). Gate first per case, fix second, one case at a time.

**Outcome up front:** CASE A is **fixed and gated (landed)**. CASE B was
implemented, verified against pinned GATK with two new strict oracles *and* the
pre-existing real-data GATK oracle, then **reverted**: its correct fix makes a
different registered expectation block in the same test file unsatisfiable, and
that block cannot be corrected to GATK truth without unrelated behavioural work
and losing the only coverage of three native features. Details in CASE B below;
the tree contains only the CASE A change.

---

## STEP 1 — measured GATK truth (not taken on faith)

Both cases were re-measured with the pinned jar on the *literal* fixtures and
flags of `verify_reblock_gvcf.py`. GATK 4.6.2.0 refuses a `-V` without `-R`
("Argument reference was missing"), so every GATK run below carries `-R` against
a synthetic 100 bp `chr1` whose bases match the fixture REF alleles (C at 6 and
20, A elsewhere; `reference.fa` + `.fai` + `.dict`; measured: GATK resolves the
dictionary as `<stem>.dict`, not `<stem>.fa.dict`). Inputs were plain
(unindexed) `.vcf`, as instructed — identical record bytes for both tools.

```
third_party/jdk17/bin/java -Xmx1g \
  -jar third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
  ReblockGVCF -R reference.fa -V input.g.vcf \
  -GQB 20 -GQB 100 --drop-low-quals --rgq-threshold 10 --floor-blocks \
  -O gatk.g.vcf
```

### CASE A — the audit fixture (`verify_reblock_gvcf.py:42-49`)

GATK 4.6.2.0, exit 0, **two** data rows:

```
chr1	1	.	A	<NON_REF>	.	PASS	END=10	GT:DP:GQ	0/0:9:0
chr1	12	.	A	G,<NON_REF>	42.64	.	DP=20;RAW_GT_COUNT=0,1,0;RAW_MQandDP=72000,20	GT:AD:DP:GQ:PL	0/1:12,8,0:20:20:50,0,80,99,99,20
```

native (before the fix), same arguments, exit 0, **three** rows:

```
chr1	1	.	A	<NON_REF>	.	.	END=10	GT:DP:AD:GQ	0/0:9:10,0:0
chr1	12	.	A	G,<NON_REF>	.	.	DP=20;RAW_MQandDP=72000,20;RAW_GT_COUNT=0,1,0	GT:DP:AD:PL:GQ	0/1:20:12,8,0:50,0,80,99,99,20:20
chr1	20	.	C	<NON_REF>	.	.	DP=5;END=20	GT:DP:GQ	0/0:5:0
```

**The audit is confirmed on the outcome:** GATK drops the `--drop-low-quals`
site at POS 20 entirely; it emits no third record. (The audit's *reason* is
consistent with the measurement, but see "what remains unproven" — from the
outside I cannot tell which of GATK's two drop paths fires.)

### CASE B — the triploid fixture (`verify_reblock_gvcf.py:269-273`)

GATK 4.6.2.0, no extra flags (`-R` + `-V` only), exit 0, **one** row:

```
chr1	30	.	A	<NON_REF>	.	.	END=30	GT:DP:GQ:MIN_DP:PL	0/0/0:17:10:17:0,10,20,30
```

native (before/after the revert), same arguments, exit 0, **one** row:

```
chr1	30	.	A	C,<NON_REF>	.	.	DP=17;RAW_MQandDP=61200,17;RAW_GT_COUNT=0,1,0	GT:DP:AD:PL:GQ	0/1/1:17:12,5,0:0,10,20,30,100,110,120,160,170,190:10
```

**The audit is confirmed exactly:** GATK reblocks to `GT 0/0/0`, `ALT <NON_REF>`,
`END=30`, and **4** PLs (`0,10,20,30` = `numLikelihoods(2, 3)`), GQ 10,
`MIN_DP 17`, `DP 17`. Note the 4 PLs are a real REF-versus-best-ALT subset
(`0,10,20,30`), not four zeros.

---

## GATK rule, with `gatk-source/` file:line

`ReblockGVCF.regenotypeVC` (`gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/variantutils/ReblockGVCF.java`):

* `:377-407` — passthrough for existing hom-ref blocks (`isHomRefBlock`, `:476`).
* `:415-424` — with `--drop-low-quals`, a concrete variant with non-zero site
  depth is **re-genotyped** through `genotypingEngine.calculateGenotypes`; a
  `null` result returns immediately, i.e. **the record is dropped**. The calling
  confidence is only armed in drop mode: `:327`
  (`standardConfidenceForCalling = dropLowQuals ? genotypeArgs.standardConfidenceForCalling : 0.0`).
* `:427-440` — only *afterwards* is `shouldBeReblocked(result)` consulted.
  `shouldBeReblocked` is `:514-535`:
  * `:520-523` `pls == null` → reblock;
  * `:524` `MathUtils.minElementIndex(pls)` — the **minimum-likelihood** entry,
    *not* the called GT;
  * `:525` `GenotypesCache.get(genotype.getPloidy(), minLikelihoodIndex)` then
    `:527` `alleleCounts.asAlleleList(vc.getAlleles())`;
  * `:528` `pls[0] < rgqThreshold` → reblock (**this is `--rgq-threshold`**);
  * `:529` `!GATKVariantContextUtils.genotypeHasConcreteAlt(finalAlleles)` → reblock
    (`GATKVariantContextUtils.java:751-757`: non-reference, non-symbolic, non-`*`);
  * `:530` the minimum-likelihood genotype contains `<NON_REF>` → reblock;
  * `:531` genotype has neither PL nor GQ → reblock; `:532` TREE_SCORE below
    threshold → reblock.
* `:537-563` `lowQualVariantToGQ0HomRef`; under drop mode it returns `null`
  (drop) unless the variant is a monomorphic hom-ref call with concrete ALTs and
  a called genotype — `:542-545`.
* `:576-631` `changeCallToHomRefVersusNonRef`: `:585-586` chooses the synthetic
  all-zero likelihood vector when genotype 0's `PL[0] != 0`, otherwise subsets
  the likelihoods to REF plus the single most likely ALT
  (`AlleleSubsettingUtils.calculateMostLikelyAlleles(vc, ploidy, 1, true)`,
  `gatk-source/.../genotyper/AlleleSubsettingUtils.java:342-361`, ranked at
  `:374-398` by summed GL difference); `:620-627` takes `DP`/`MIN_DP` from INFO
  `DP` when present, else from `sum(AD)`; `:627` sets `END`; `:629` sets GT to
  the reference repeated at the sample ploidy.

**Native's corresponding divergence (CASE A):** native evaluated the
`--rgq-threshold` / PL[0] conversion *before* its own drop-mode
re-genotyping emulation — the reverse of GATK's `:415` before `:427`. The
POS-20 site therefore got converted into a GQ0 block instead of being dropped.

**Native's corresponding divergence (CASE B):** native classified
variant-versus-block from the record's ALT list
(`reference_block = has_non_ref(record) && !has_concrete_alt(record)`), never
from the minimum-likelihood PL genotype, so clause `:529` was absent.

---

## CASE A — FIXED (landed)

**Gate before the fix** — `fastgatk-native/scripts/verify_reblock_gvcf_droplowqual_gatk_oracle.py`
(new strict gate, modelled on `verify_variant_filtration_flag_only_gatk_oracle.py`):

```
GATK   rows (1):  chr1 1 . A <NON_REF> . . END=10 GT:DP:GQ  0/0:10:20
NATIVE rows (2):  [same row]
                  chr1 20 . C <NON_REF> . . DP=5;END=20 GT:DP:GQ 0/0:5:0
VIOLATION: record count differs: GATK=1 native=2
EXIT=1
```

(gated fixture: the surviving reference block written in the `GT:DP:GQ` form both
tools encode identically, plus the POS-20 site copied verbatim from
`verify_reblock_gvcf.py:48`; flags identical to the test's main case. A second
gated control case `--drop-low-quals` without `--rgq-threshold` passes before and
after, so the gate is not vacuous. The literal four-record audit fixture is
carried as a reported-only diagnostic case.)

**Fix applied** — `fastgatk-native/src/reblock_gvcf_tool.cpp`, one reordering
inside the record loop (no other tool touched):

* the drop-mode re-genotyping block (`--drop-low-quals` + concrete ALT) now runs
  *before* the `convert_low_quality` (PL[0] < `--rgq-threshold`) conversion,
  mirroring `ReblockGVCF.java:415` before `:427`;
* `convert_low_quality` and the conversion branch now test the refreshed
  `now_block` rather than the stale `reference_block`, so a record the
  re-genotyping step already turned into a block is not converted twice;
* the new `regenotyped_to_block` flag keeps the pre-existing
  "drop records with no GQ under `--drop-low-quals`" branch reachable exactly as
  before for records the re-genotyping step did *not* convert (GATK keeps a
  confidently-hom-ref site in that situation, so this is also closer to GATK).

**Gate after the fix** — exit 0, rows byte-identical:

```
[drop-low-quals-with-rgq-threshold]  GATK rows / NATIVE rows:
    chr1 1 . A <NON_REF> . . END=10 GT:DP:GQ 0/0:10:20        (both, 1 row each)
[drop-low-quals-alone]               both: chr1 1 . A <NON_REF> . . END=10 GT:DP:GQ:MIN_DP 0/0:10:20:10
```

**Registered test corrected** (`fastgatk-native/scripts/verify_reblock_gvcf.py`):

| line (before) | change |
| --- | --- |
| 65 | `assert len(records) == 3` → `assert len(records) == 2` + comment citing `ReblockGVCF.java:415-424`, `:427`, `:542-545` |
| 81-83 | `assert records[2][4] == "<NON_REF>"` / `assert "END=20" in records[2][7]` → `assert [record[1] for record in records] == ["1", "12"]` with the GATK-source comment and the measured row set |

Nothing else in that file depends on the dropped record: the later runs reuse
the same source but (a) the plain-output query selects `chr1:12-12`, (b) the
interval/INTERSECTION/tree-score/QUALapprox/multisite blocks omit
`--drop-low-quals` and are unaffected, and (c) `summary["merged_blocks"] == 1`
still holds.

---

## CASE B — implemented, verified, REVERTED (not fixed)

**Gate before the fix** —

```
GATK   rows (1): chr1 30 . A <NON_REF> . . END=30 GT:DP:GQ:MIN_DP:PL 0/0/0:17:10:17:0,10,20,30
NATIVE rows (1): chr1 30 . A C,<NON_REF> . . DP=17;RAW_MQandDP=61200,17;RAW_GT_COUNT=0,1,0 GT:DP:AD:PL:GQ 0/1/1:17:12,5,0:0,10,20,30,100,110,120,160,170,190:10
VIOLATION: row 1 is not byte-identical ...
EXIT=1
```

**Fix implemented and measured** (in `reblock_gvcf_tool.cpp`, then removed):

* a `pl_min_genotype_has_concrete_alt()` helper decoding
  `minElementIndex(PL)` into the genotype at the sample ploidy (the htsjdk
  colexicographic `GenotypeAlleleCounts` index, inverted via
  `rank = Σ_j C(a_j + j, j + 1)`) and testing `genotypeHasConcreteAlt`
  (`ReblockGVCF.java:524-530`);
* a clause in the record loop invoking `changeCallToHomRefVersusNonRef`'s two
  branches per `:585-586` — `PL[0] != 0` → `convert_to_ref_block` (synthetic
  zeros), `PL[0] == 0` → `convert_high_confidence_homref_to_ref_block`
  (REF-versus-best-ALT subsetting);
* a `record_ploidy` derived from the GT width.

Measured with that change in place:

* case-B gate `verify_reblock_gvcf_triploid_gatk_oracle.py`: **exit 0**, both
  cases byte-identical, including the surviving-block control row —
  `chr1 30 . A <NON_REF> . . END=30 GT:DP:GQ:MIN_DP:PL 0/0/0:17:10:17:0,10,20,30`;
* case-A gate: still **exit 0**;
* `fastgatk-reblock-gvcf-gatk-oracle` (real HaplotypeCaller gVCF, 17 sub-cases):
  **Passed** — the faithful clause does not move any real-data record;
* `fastgatk-reblock-gvcf-overlap-gatk-oracle`, `-shards-gatk-oracle`,
  `-java-benchmark-command-contract`: **Passed**.

**Why it was reverted.** With the clause in place,
`verify_reblock_gvcf.py` fails at its *third* fixture, the "reverse trimming /
deletion gap / NON_REF AD" block (`:441-480`), because that fixture has
`PL = 0,10,20,...,90` with `PL[0] = 0` — the minimum-likelihood genotype is
`0/0`, so GATK reblocks it too and never takes the `cleanUpHighQualityVariant`
path the block asserts. Measured with the fix, GATK on that fixture gives

```
chr1	50	.	A	<NON_REF>	.	.	END=51	GT:DP:GQ:MIN_DP:PL	0/0:15:10:15:0,10,20
```

while native gives the same shape but `0/0:19:10:19:0,10,20` (native keeps the
FORMAT `DP=19`; GATK prefers INFO `DP=15`, `ReblockGVCF.java:620-627`). Making
that block pass would require **all** of:

1. a third behavioural change beyond the two pinned cases (the INFO-DP
   preference rule) in the shared ref-block conversion helper;
2. deleting the block's `reverse_allele_trimmed_variants == 1`,
   `deletion_gap_ref_blocks == 1`, `non_ref_ad_zeroed == 1` and
   `non_ref_ad_kernel_*` telemetry expectations (they all become 0 for that
   input) — `verify_reblock_gvcf.py` is the **only** test in the tree covering
   reverse allele trimming, deletion-gap reference blocks and NON_REF AD
   zeroing (checked by grep over `fastgatk-native/scripts` and `tests`);
3. or re-designing that fixture's PL vector and re-deriving ~25 assertion
   values — and even then the block is not GATK-true: measured on a corrected
   `PL` vector, GATK and native still differ on the variant row in FORMAT order
   (`GT:AD:DP:GQ:PL` vs `GT:DP:AD:PL:GQ`), INFO order, `DP` (19 vs 15) and on
   the gap block's QUAL/FORMAT (`GT:DP:GQ:MIN_DP:PL` vs `GT:PL:GQ:DP:MIN_DP`,
   native `QUAL=0`).

Per the round's rule ("If a case's correct fix is not minimal or not safe, STOP
that case"), CASE B is therefore **not landed**. The code change was removed by
targeted reverse edits (helpers, the `record_ploidy` computation, the loop
clause) and both trees rebuilt; `verify_reblock_gvcf.py` keeps its original
triploid assertions, now annotated as a KNOWN DIVERGENCE whose strict gate
fails until the fix lands.

**Recommendation for a follow-up round.** Land the clause as written above
(it is verified correct), and in the same change:

* correct the trim fixture's `PL` so its minimum-likelihood genotype contains a
  concrete ALT (for example `10,0,20,30,40,50,60,70,80,90`), keeping the
  reverse-trim/gap-block/NON_REF-AD coverage alive, and re-derive its
  assertions from a GATK measurement;
* implement `ReblockGVCF.java:620-627` (INFO `DP` wins over the genotype `DP`
  for `DP`/`MIN_DP` in `changeCallToHomRefVersusNonRef`) — the same fixture
  exposes it;
* treat the FORMAT/INFO ordering and gap-block QUAL differences as their own
  items; they are pre-existing and independent of this rule.

---

## STEP 5 — gate results

| step | command | result |
| --- | --- | --- |
| (a) CASE A oracle | `python3 fastgatk-native/scripts/verify_reblock_gvcf_droplowqual_gatk_oracle.py` | **exit 0** (was exit 1 before the fix) |
| (a) CASE B oracle | `python3 fastgatk-native/scripts/verify_reblock_gvcf_triploid_gatk_oracle.py` | **exit 1 — still failing (reported, not registered as passing)** |
| (b) contract | `ctest --test-dir fastgatk-native/build -R 'reblock' -V` | **5/5 Passed** (omp) |
| (b) contract | `ctest --test-dir fastgatk-native/build-serial -R 'reblock' -V` | **5/5 Passed** (serial) |
| (c) registered gates | `ctest -R 'window-invariance-gatk-oracle\|alleles-overlap-gate-oracle\|span-del-qual-gatk-oracle\|gvcf-symbolic-prior-gatk-oracle\|af-zero-format-gatk-oracle\|arbitrary-ploidy-span-del-prior-oracle\|polyploid-gvcf-span-del-prior-oracle\|spanning-prior-genotype-gq-oracle\|multialt-owner-annotation-oracle\|mutect2-recheck\|gvcf-indel-end-gatk-oracle\|culprit-gatk-oracle\|select-variants-refonly-gatk-oracle\|asfilterstatus-gatk-oracle\|variant-filtration-flag-only-gatk-oracle' -V` | **17/17 Passed, exit 0** |
| (d) regression | `fastgatk-native/scripts/run_regression.sh --label reblock-gvcf` | **exit 0 — omp 298/298, serial 298/298** |

Step (c) output (`ctest` exit 0, all `Passed`):

```
 1/17 #83  fastgatk-hc-window-invariance-gatk-oracle .............. Passed  67.00
 2/17 #84  fastgatk-hc-ploidy-window-invariance-gatk-oracle ...... Passed  31.65
 3/17 #85  fastgatk-hc-alleles-overlap-gate-oracle ............... Passed  31.65
 4/17 #86  fastgatk-hc-span-del-qual-gatk-oracle ................. Passed 104.95
 5/17 #87  fastgatk-hc-gvcf-symbolic-prior-gatk-oracle ........... Passed  98.35
 6/17 #88  fastgatk-hc-arbitrary-ploidy-span-del-prior-oracle .... Passed  67.82
 7/17 #89  fastgatk-hc-polyploid-gvcf-span-del-prior-oracle ...... Passed  64.83
 8/17 #90  fastgatk-hc-spanning-prior-genotype-gq-oracle ......... Passed  81.11
 9/17 #91  fastgatk-hc-af-zero-format-gatk-oracle ................ Passed 117.32
10/17 #92  fastgatk-hc-multialt-owner-annotation-oracle .......... Passed  68.76
11/17 #93  fastgatk-mutect2-recheck-normal-replay-oracle ......... Passed  54.96
12/17 #94  fastgatk-mutect2-recheck-joint-assembly-oracle ........ Passed  13.21
13/17 #95  fastgatk-hc-gvcf-indel-end-gatk-oracle ................ Passed  36.73
14/17 #96  fastgatk-variant-recalibrator-culprit-gatk-oracle ..... Passed  46.40
15/17 #97  fastgatk-select-variants-refonly-gatk-oracle .......... Passed 138.48
16/17 #98  fastgatk-variant-filtration-asfilterstatus-gatk-oracle  Passed  39.51
17/17 #99  fastgatk-variant-filtration-flag-only-gatk-oracle ..... Passed  42.84
```

Step (d) output (`fastgatk-native/scripts/run_regression.sh --label reblock-gvcf`,
run once; exit 0; evidence directory `.diag/regression/20260911-213607`,
`FASTGATK_REQUIRE_GATK_ORACLE=1` set by the runner, no staleness warning emitted):

```
| 后端   | 构建目录                        | 结果 | 通过/总数 | 耗时         |
| omp    | OpenMP (fastgatk-native/build)  | 通过 | 298/298  | 1297.97 sec  |
| serial | Serial (fastgatk-native/build-serial) | 通过 | 298/298 | 1324.99 sec |
```

Final re-confirmation after the run: `verify_reblock_gvcf.py` exit 0
(`{"status": "pass", "output_records": 2, "merged_blocks": 1}`), CASE A oracle
exit 0, CASE B oracle exit 1 with 2 violations (unchanged, as required).

`ctest` is not on `PATH` in this environment; the bundled
`third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest` was used.

### Tree state

* `fastgatk-native/src/reblock_gvcf_tool.cpp` — **CONTAINS** the CASE A
  reordering (see `git diff`); the CASE B helpers/clause are **absent**.
* `fastgatk-native/scripts/verify_reblock_gvcf.py` — **CONTAINS** the CASE A
  assertion correction; the CASE B block is the original assertion set plus a
  KNOWN-DIVERGENCE comment.
* `fastgatk-native/scripts/verify_reblock_gvcf_droplowqual_gatk_oracle.py` —
  new file, present, passing.
* `fastgatk-native/scripts/verify_reblock_gvcf_triploid_gatk_oracle.py` —
  new file, present, **failing** (deliberately, per step 5a).
* Neither new script is registered in `CMakeLists.txt` (constraint 1 reserves
  registration to the orchestrator).

---

## What remains unproven (and what is speculation)

* **Speculation:** which of GATK's two drop paths removes the POS-20 site.
  `ReblockGVCF.java:415-424` (re-genotyping returns `null`) and `:542-545`
  (`lowQualVariantToGQ0HomRef` returns `null` after re-genotyping left no
  concrete ALT) both produce the measured output, and I found no fixture that
  separates them externally. The audit's `:542-545` citation is therefore
  plausible but not established; what *is* established is the ordering,
  `:415` before `:427`.
* **Unproven:** native's drop-mode re-genotyping is a heuristic
  (PL[0] sum < 30, PL-derived GQ < 30), not GATK's genotyping engine. It
  reproduces GATK on every fixture exercised here (the real HaplotypeCaller
  oracle included), but the threshold rule itself is an approximation and was
  not re-derived in this round.
* **Unproven:** CASE B's fix was verified only on the fixtures in the two new
  gates plus the real-data reblock oracle; it was not run through
  `run_regression.sh`, because it was reverted. A follow-up round that lands it
  must run the full suite.
* **Known, untouched divergences exposed while measuring** (all pre-existing and
  independent of the two pinned cases): native leaves `QUAL` where GATK writes
  the re-genotyped value on drop-mode variants; INFO and FORMAT key ordering
  differ; ref-block `FILTER` is `.` in native where GATK writes `PASS` when the
  input record did; the `--rgq-threshold`-only conversion emits
  `INFO=DP=5;END=20` / `GT:DP:PL:GQ:MIN_DP` where GATK emits `END=20` /
  `GT:DP:GQ:MIN_DP:PL`; block `DP`/`MIN_DP` prefer the FORMAT field where GATK
  prefers INFO `DP`.
