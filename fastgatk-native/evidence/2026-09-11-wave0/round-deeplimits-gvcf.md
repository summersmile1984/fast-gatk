# Round: `verify_hc_alleles_deep_limits.py` gVCF field-level divergences

Date: 2026-09-11. Repo: `/home/turing-agents/Documents/fast-gatk` (worktree only; no commit, no
branch, no CMake edit, no root `*.md` touched).

Scope: the two non-exit-code cases that `fastgatk-native/scripts/verify_hc_alleles_deep_limits.py`
reported as `diverged` with **equal row counts**, i.e. pure field-level gVCF differences.

| # | case | outcome of this round |
|---|------|------------------------|
| A | `gvcf-homopolymer-insertion` (forced homopolymer insertion beside `<NON_REF>`) | **NOT CONTAINED** — root cause found, minimal fix written and proven by a new gate, then **reverted** because it reddens the registered native contract test `fastgatk-indel-contract` (§1.5). |
| B | `gvcf-max-alt-alleles-1` (gVCF reference-confidence blocks under a forced-alt budget) | **NOT FIXED** — reproduced, localized with four controls, reported with a recommendation; its gate stays red (§2). |
| — | `max-genotype-count-2` | out of scope by instruction (GATK exit 3 / native exit 0); untouched. |

Final tree state: **no production-source change survives**; two new strict gates and the round's
diagnostic scratch are added, and `verify_hc_alleles_deep_limits.py` still exits 1 with
`diverged_cases = ['max-genotype-count-2', 'gvcf-homopolymer-insertion', 'gvcf-max-alt-alleles-1']`,
i.e. exactly the round's starting state.

Tooling: pinned `third_party/jdk17/bin/java -Xmx1g -jar
third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar HaplotypeCaller`;
`fastgatk-native/build/fastgatk-hc-call` (Kokkos/OpenMP) and
`fastgatk-native/build-serial/fastgatk-hc-call` (Serial). Every oracle builds its own fixture under
`tempfile.TemporaryDirectory`.

---

## 1. Target A — `gvcf-homopolymer-insertion`: native-only `INFO/END` on a concrete-indel record

### 1.1 Exact field-level diff (reproduced pre-fix)

`verify_hc_alleles_deep_limits.py --case gvcf-homopolymer-insertion` reports
`first_diff_row=6`, `first_diff_field=7` (INFO), 37 rows on both sides. Only row index 6 differs:

```
GATK    chr1  697  .  TC  T,<NON_REF>  98.60  .  BaseQRankSum=0.000;DP=36;ExcessHet=0.0000;MLEAC=1,0;MLEAF=0.500,0.00;MQRankSum=0.000;RAW_MQandDP=129600,36;ReadPosRankSum=3.417   GT:AD:DP:GQ:PL:SB  0/1:30,6,0:36:99:106,0,961,196,979,1176:30,0,6,0
NATIVE  chr1  697  .  TC  T,<NON_REF>  98.60  .  END=698;BaseQRankSum=0.000;DP=36;ExcessHet=0.0000;MLEAC=1,0;MLEAF=0.500,0.00;MQRankSum=0.000;RAW_MQandDP=129600,36;ReadPosRankSum=3.417  GT:AD:DP:GQ:PL:SB  0/1:30,6,0:36:99:106,0,961,196,979,1176:30,0,6,0
```

* **POS `chr1:697`; column INFO; key `END`; GATK: absent; native: `END=698`.**
* QUAL `98.60`, every other INFO key and its order, all six FORMAT keys and every value, the row
  order and all 36 other rows (including each `<NON_REF>` block and its own `END=`) are identical.

The same single-key difference reproduces **without `--alleles`** (the assembled `chr1:697 TC>T`
deletion) and on a forced insertion in the tandem-repeat window (`END=900`), so it is not an
`--alleles` divergence — it is the previously recorded `O1` in `.diag/track-b-alleles-findings.md`.

### 1.2 GATK's rule (authoritative, quoted)

* A gVCF run hands ordinary `VariantContext`s to the block combiner; a *variant* record is emitted
  verbatim:
  `gatk-source/src/main/java/org/broadinstitute/hellbender/utils/variant/writers/GVCFBlockCombiner.java:203-209`
  — `} else { emitCurrentBlock(); nextAvailableStart = vc.getEnd(); contigOfNextAvailableStart = vc.getContig(); toOutput.add(vc); }`.
* `END` is only ever *set* on hom-ref blocks:
  `.../utils/variant/writers/GVCFBlock.java:53` — `vcb.attribute(VCFConstants.END_KEY, getEnd());`.
* The combiner only *reads* `END` when feeding a hom-ref site into a block
  (`GVCFBlockCombiner.java:131`: `currentBlock.add(vc.getStart(), vc.getAttributeAsInt(VCFConstants.END_KEY, vc.getStart()), g);`),
  i.e. `VariantContext.getEnd()` (which does span the deletion) is used privately for contiguity and
  is **not serialized**.
* No other `END`-setting site exists in the HaplotypeCaller path (`ReblockingGVCFBlockCombiner`,
  `CombineGVCFs`, `ReblockGVCF`, `VariantDataManager` and the SV/CNV composers are the only other
  writers of the key).

### 1.3 Native site and the minimal fix (written, gated, then reverted — see §1.5)

`fastgatk-native/src/hc_call.cpp`'s two concrete-gVCF-candidate writers each contained an
unconditional END emission (`has_indel` = REF length != ALT length):

```cpp
if (!bp_resolution && has_indel) { out << "END=" << reference_end; info_started = true; }
```

The minimal fix (12 insertions / 12 deletions, saved verbatim at
`.diag/scratch-deeplimits/end_fix_reverted.diff`) removed both blocks plus the then-unused
`has_indel` helper and `reference_end` locals, and added a comment citing the rule above. Nothing
else was touched; the hom-ref block writer's own `END=` (which GATK does emit, and which
`verify_hc_chr20_min_pruning_gvcf_gatk_oracle.py:81` pins) stays.

### 1.4 Gate (new) and its measured results

New strict oracle `fastgatk-native/scripts/verify_hc_gvcf_indel_end_gatk_oracle.py`
(`main()->int`, `--native` / `FASTGATK_HC_BINARY`, `--expect-divergence`, `TemporaryDirectory`
scratch, final JSON status line; fixture shared with `verify_hc_alleles_deep_boundary`).
Four cases are run, **three gated with no canonicalization of any kind**, one probe (§4):

| gate run | exit | literal result |
|----------|------|----------------|
| before any fix | **1** | `gvcf-homopolymer-insertion: row 6 POS=697: INFO END None -> '698'`; `gvcf-assembled-deletion: row 6 POS=697: INFO END None -> '698'`; `gvcf-tandem-insertion: row 27 POS=900: INFO END None -> '900'` |
| with the §1.3 fix applied (OpenMP **and** Serial) | **0** | all three cases `data_rows_byte_identical=True` (37/37, 37/37, 29/29) |
| final tree (fix reverted) | **1** | the same three violations as "before any fix" |

With the fix applied, `verify_hc_alleles_deep_limits.py` also moved that case to `match`
(37 vs 37 rows) — i.e. the fix demonstrably closes target A.

### 1.5 Why the fix is **not** in the tree (criterion failure → revert, as instructed)

The full dual-backend regression **with** the fix applied was **292/293 on both backends**
(`.diag/regression/20260911-115942/{omp,serial}.log`), failing exactly one test:

```
242/293 Test #51: fastgatk-indel-contract ...................................***Failed
  File ".../scripts/verify_indel.py", line 173, in main
    assert int(info["END"]) == 11
KeyError: 'END'
```

`fastgatk-native/scripts/verify_indel.py:173` is a **native-internal contract** on a synthetic
`17:1-32` `10M1D10M` fixture driven with native-only flags (`--min-candidates` etc.; no GATK run in
that block). It asserts that the concrete deletion candidate *does* carry `END=11`. That assertion
dates from the baseline commit (`180c2bf`, where the file is added whole) and pins native's own — now
known to be non-GATK — behaviour; the neighbouring reference-block overlap assertion
(`verify_indel.py:174-180`) does **not** depend on the candidate's `END`.

*GATK-verified* behaviour is the opposite of that assertion — the new gate above measures it against
the pinned jar on two fixtures. Still, this round may not edit or weaken an existing oracle
(constraint 3), and step 4(d) requires 293/293 on both backends. The change was therefore reverted
with a targeted restore of the single file that carried it
(`git checkout -- fastgatk-native/src/hc_call.cpp`; no other edit lived in that file), both backends
were rebuilt, and the regression was re-run: **293/293 on both backends** (§3).

**Recommendation for target A.** Land the saved patch together with an update of the stale
`verify_indel.py:173` assertion (assert that the candidate carries *no* `END`, matching GATK, while
`verify_indel.py:174-180` continues to assert the reference-block span exclusion). That is a two-file
change whose second half is outside this round's authority, which is why it is reported rather than
applied. The gate `verify_hc_gvcf_indel_end_gatk_oracle.py` is ready to register as soon as that
happens.

---

## 2. Target B — `gvcf-max-alt-alleles-1`: reference-block `PL`/`GQ` under forced alleles (NOT FIXED)

### 2.1 Exact field-level diff (unchanged by this round)

`verify_hc_alleles_deep_limits.py --case gvcf-max-alt-alleles-1` (42 rows on both sides) reports only
the first differing field (`first_diff_row=6`, `first_diff_field=9`, the sample column). The complete
field-level diff:

```
row  6  POS=621  GATK   chr1 621 . A <NON_REF> . . END=639  GT:DP:GQ:MIN_DP:PL  0/0:36:99:36:0,108,1612
                 NATIVE chr1 621 . A <NON_REF> . . END=639  GT:DP:GQ:MIN_DP:PL  0/0:36:90:36:0,90,1350
row  8  POS=641  GATK   0/0:36:99:36:0,108,1612            NATIVE 0/0:36:90:36:0,90,1350
row 10  POS=661  GATK   0/0:36:99:36:0,108,1612            NATIVE 0/0:36:90:36:0,90,1350
row 11  POS=697  GATK   ... BaseQRankSum=0.000;DP=36;...   NATIVE ... END=698;BaseQRankSum=0.000;...
```

* POS 621/641/661, **sample column**: `FORMAT GQ` **99 → 90**, `FORMAT PL` **`0,108,1612` → `0,90,1350`**.
  `GT`, `AD`, `DP`, `MIN_DP`, the block boundaries (`END=639/659/696`) and every INFO field of those
  rows are identical.
* POS 697, INFO `END` — the target-A defect, which §1.3 fixes (not in the tree).
* The three differing rows are exactly the reference blocks *between* the forced SNPs
  (`chr1:620/640/660`); the blocks before the first forced SNP and after the deletion agree.

Both vectors are members of GATK's own model family, so this is *model selection*, not formatting:
`0,90,1350` is the indel-cache row for `nInformativeReads = 30` (`10*log10(2)*30 ≈ 90.3`,
`45*30 = 1350`); `0,108,1612` is the 36-read Ref-vs-Any SNP row (`3.010156*36 ≈ 108.4`,
`44.7712*36 ≈ 1611.8`) — **not** the 36-read indel row, which would be `0,108,1620`. Both callers
emit `0,90,1350` on the same fixture **without** `--alleles`, and the divergence also appears without
`--max-alternate-alleles` (case `forced-snps-gvcf` of the new gate), so the alt budget is not what
triggers it.

### 2.2 GATK's rule (authoritative, quoted)

`gatk-source/.../walkers/haplotypecaller/ReferenceConfidenceModel.java:301-317`
(`doIndelRefConfCalc`) keeps the **least confident** of the SNP Ref-vs-Any model and an independent
indel model keyed on the number of indel-informative reads:

```java
final GenotypeLikelihoods snpGLs = GenotypeLikelihoods.fromLog10Likelihoods(
        homRefCalc.getGenotypeLikelihoodsCappedByHomRefLikelihood());
final int nIndelInformativeReads = calcNReadsWithNoPlausibleIndelsReads(
        pileup, refOffset, ref, indelInformativeDepthIndelSize);
final GenotypeLikelihoods indelGLs = getIndelPLs(ploidy, nIndelInformativeReads);
final GenotypeLikelihoods leastConfidenceGLs = getGLwithWorstGQ(indelGLs, snpGLs);
homRefCalc.finalPhredScaledGenotypeLikelihoods = leastConfidenceGLs.getAsPLs();
```

`getGLwithWorstGQ` (`ReferenceConfidenceModel.java:331-342`) returns the indel model only when its
hom-ref log-GQ is the *larger* (less confident) one; at depth 36 the two models differ by ~0.001
phred, so which vector is emitted is decided entirely by `nIndelInformativeReads` (30 → indel row,
36 → SNP row). The count comes from `readHasNoPlausibleIdealsOfSize` /
`traverseEndOfReadForIndelMismatches` (`ReferenceConfidenceModel.java:594-768`, driven from
`calcNReadsWithNoPlausibleIndelsReads`, lines 785-800), evaluated over the **RCM pileup's reads**,
which GATK builds from `readLikelihoods` *after*
`AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype` (`HaplotypeCallerEngine.java:964` →
`calculateRefConfidence(..., readLikelihoods, ...)`, line 1019; `getPileupsOverReference(...,
readLikelihoods, samples)`, `ReferenceConfidenceModel.java:228`).

### 2.3 Localization (measured this round)

Native's counterparts are `aligned_view_is_indel_informative` / `read_is_indel_informative`
(`fastgatk-native/src/calling_pipeline.cpp:2403-2496`) reached from
`calculate_reference_confidence` (`calling_pipeline.cpp:3144-3281`), fed per reference-confidence
segment by `append_segment` (`calling_pipeline.cpp:13543-13561`).

Diagnosis used temporary instrumentation (env var `FASTGATK_DEBUG_RCM_WINDOW`, **since reverted**;
`git diff fastgatk-native/src/calling_pipeline.cpp` is empty). Measurements:

* native's per-locus `nInformativeReads` in the `--alleles` run drops **36 → 30 exactly at
  `chr1:620`**, the first locus of the event-trimmed segment whose reference window is `[599,775]`;
  the six reads that flip are exactly the fixture's six 1 bp-frameshift deletion reads
  (records 30-35, at `read_offset = pos - 560`, i.e. their raw-alignment offsets).
* in the same fixture **without** `--alleles`, native's drop is at `chr1:644`, inside a segment with
  window `[399,752]` — and GATK agrees there (both emit `0,90,1350`, blocks `568-643`/`644-696`).
* native's own summary JSON differs between the two runs exactly in the realignment accounting:
  `rcm_realigned_observations` = **6119** with `--alleles` versus **5435** without (both with
  `rcm_haplotype_realignment_used = true`, `pairhmm_used = true`,
  `pairhmm_skip_reason = executed`), consistent with the forced alleles changing the best-haplotype
  assignment of those six reads.
* GATK's own window for the alleles run — **derived from the source, not observed** (GATK exposes no
  hook for it) — is `paddedVariantSpan` = `[599,775]` (`AssemblyRegionTrimmer.trim`,
  `AssemblyRegionTrimmer.java:186-204`; `Result.getVariantRegion()` →
  `originalRegion.trim(variantSpan, paddedVariantSpan)`), indexed by
  `curPos - activeRegion.getPaddedSpan().getStart()` (`ReferenceConfidenceModel.java:233/239`,
  `globalRefOffset`). **If that derivation is right, the reference windows agree**, and the residual
  difference must be in the reads or in the traversal — the assumption behind H1/H2 below.

Four controls isolate it (all four are cases of the gate in §2.4):

| control | measured result |
|---------|-----------------|
| `--indel-size-to-eliminate-in-ref-model 0` (plausible-indel search vacuous in both implementations) | **byte-identical** — both sides move to `0,108,1612`; only the (separate) `END` differed |
| same fixture with the six deletion reads removed | **byte-identical** (21 rows) |
| same six reads re-encoded as a real `139M1D161M` CIGAR deletion | **byte-identical** (only `END` differed) |
| same fixture **without** `--alleles` | **byte-identical** blocks (both `0,90,1350`) |

So the divergence is confined to the plausible-indel informativeness classification of six
pathologically-encoded (frameshift) reads, reached through the reference-confidence segmentation that
the forced alleles select.

**Why no fix was applied.** The two candidate mechanisms this round did **not** separate:

* **H1 (leading; structural evidence).** Native evaluates the informativeness test on the *original*
  read payload while the observations carry *realigned* offsets: `calling_pipeline.cpp:17146-17197`
  remaps `reference_confidence_observations` onto `pairhmm.rcm_realigned_observations`, taking each
  base from `active_reads.bases[begin + observation->read_offset]` (the raw read), and
  `calculate_reference_confidence` then calls `read_is_indel_informative(reads, references, record,
  observations.read_offset[observation], ...)` (`calling_pipeline.cpp:3170`) on that same raw
  `ReadBatch`. GATK runs the test on the *realigned* reads (`HaplotypeCallerEngine.java:964` →
  `ReferenceConfidenceModel.java:228/795`), whose CIGAR for a frameshift read can carry the deletion
  — `AlignmentUtils.getBasesAndBaseQualitiesAlignedOneToOne` then inserts a gap character, and
  `isMismatchAndNotAnAlignmentGap` ignores gaps — which is exactly what turns those reads
  informative. Native has a purpose-built path for this (`projected_read_is_indel_informative`,
  `calling_pipeline.cpp:2750`), but the per-segment `realigned_projections` map is **empty** in this
  run (`[FASTGATK_RCM_CALL] ... projections=0` for all four segments, with and without `--alleles`),
  so the raw-read path is always taken.
* **H2 (not excluded; speculation).** A transcription detail of the traversal itself that only
  becomes observable with this window/read combination. A line-by-line comparison of
  `aligned_view_is_indel_informative` against `traverseEndOfReadForIndelMismatches` showed no obvious
  mismatch (loop bounds, the `referenceWasShorter` branch selection and last-base reset, the flip
  range and the early-exit conditions all correspond), but the traversal was **not** re-derived
  numerically for the failing (read, position, window) triple, so H2 cannot be excluded.

Forcing either vector, special-casing frameshift reads, or retuning the segment window so this fixture
passes would break the passing controls in §2.4 and in `verify_hc_span_del_qual_oracle.py` /
`verify_hc_gvcf_symbolic_prior_gatk_oracle.py` / the polyploid gVCF gates — i.e. it would not be
minimal or safe. **Recommendation:** pursue H1 — make the informativeness test consume the same
realigned read view that the observation remap uses (or populate `realigned_projections` from
`pairhmm.rcm_realigned_observations` when the owner path does not) — and keep this oracle as the
acceptance test. That is a change in the RCM↔realignment plumbing, hence reported rather than
attempted.

### 2.4 Gate (new) and its exit statuses

New strict oracle `fastgatk-native/scripts/verify_hc_gvcf_forced_allele_refblock_gatk_oracle.py`
(same conventions; no canonicalization; prints both sides' reference-block `PL` lists and the per-row
field diff). Four gated cases:

| case | before | final tree |
|------|--------|------------|
| `forced-snps-gvcf` (3 forced SNPs, `-ERC GVCF`) | FAIL: rows 6/8/10 `FORMAT GQ '99' -> '90'`, `FORMAT PL '0,108,1612' -> '0,90,1350'`; row 11 `INFO END None -> '698'` | **FAIL** — same five rows |
| `forced-snps-gvcf-max-alt-1` (the deep-limits case) | FAIL: identical diff | **FAIL** — identical diff |
| `forced-snps-gvcf-indel-size-0` (control) | FAIL only on row 11 `END` | **FAIL** only on row 11 `END` (passes once §1.3 lands) |
| `forced-snps-gvcf-nodeletion` (control) | **PASS** (21/21) | **PASS** (21/21) |

Exit status: **1 before, 1 at the final tree** — the honest state for an unfixed target.
**Do not register this gate as a passing test yet**; `--expect-divergence` gives the diagnostic form
that always exits 0.

---

## 3. Step 4 — gate results

All numbers below are for the **final tree state** (reverted production source, both backends rebuilt
from it).

| step | command | result |
|------|---------|--------|
| 4a | `cmake --build fastgatk-native/build --target fastgatk-hc-call -j 16` | **exit 0** (57 s with the fix applied; rebuilds after the revert also exit 0) |
| 4b | `ctest --test-dir fastgatk-native/build -R '<the 10-name regex>' -V` run standalone | **exit 0 — 12/12 passed**, 725 s (`.diag/scratch-deeplimits/ctest_step4b.log`) |
| 4b (re-confirmed) | the same 12 tests inside the final full regression | **12/12 Passed on both backends**: `fastgatk-hc-window-invariance`, `-ploidy-window-invariance`, `-alleles-overlap-gate`, `-span-del-qual`, `-gvcf-symbolic-prior`, `-arbitrary-ploidy-span-del-prior`, `-polyploid-gvcf-span-del-prior`, `-spanning-prior-genotype-gq`, `-af-zero-format`, `-multialt-owner-annotation`, `mutect2-recheck-normal-replay`, `mutect2-recheck-joint-assembly` |
| 4c | `python3 fastgatk-native/scripts/verify_hc_alleles_gatk_oracle.py` | **exit 0**, `"status": "pass"`, `"multi_alt_and_anchored_indel_gatk_exact": true` |
| 4d | Serial `cmake --build fastgatk-native/build-serial --target fastgatk-hc-call -j 16`, then `run_regression.sh --label 'deep-limits gVCF 分诊后（A 项回滚，B 项未修）' -j 12` | **exit 0 — omp 293/293, serial 293/293** (901.95 s / 905.55 s); `.diag/regression/20260911-121618/{omp,serial}.log` |

Transient measurement recorded while the §1.3 fix was applied (the reason it was reverted):
`run_regression.sh` = **292/293 on both backends**, sole failure `fastgatk-indel-contract`
(`verify_indel.py:173`, `KeyError: 'END'`), `.diag/regression/20260911-115942/{omp,serial}.log`.

New gates at the final tree state: `verify_hc_gvcf_indel_end_gatk_oracle.py` **exit 1** (3 gated
violations, 1 probe note); `verify_hc_gvcf_forced_allele_refblock_gatk_oracle.py` **exit 1**
(2 divergent cases, 2 passing controls); `verify_hc_alleles_deep_limits.py` **exit 1** with the
round's starting `diverged_cases`.

---

## 4. New, unrelated observation found while gating (reported, not fixed, not hidden)

The `bp-resolution-control` case (same forced-homopolymer-insertion fixture at
`-ERC BP_RESOLUTION`) exposes a **pre-existing** divergence unrelated to the `END` rule and to
target B — on the *reference* row for the base deleted by the six frameshift reads:

```
POS=698  GATK   ... GT:AD:DP:GQ:PL  0/0:30,6:36:0:0,0,1135
         NATIVE ... GT:AD:DP:GQ:PL  0/0:30,0:30:90:0,90,1343
```

GATK counts the six deletion reads as covering the deleted base (`AD=30,6`, `DP=36`, hence `GQ=0`);
native reports only the 30 reference reads (`AD=30,0`, `DP=30`, `GQ=90`). It is recorded as an
out-of-scope **probe note** in `verify_hc_gvcf_indel_end_gatk_oracle.py` so that gate stays a strict
oracle for the rule it names. It is not covered by `verify_hc_alleles_deep_limits.py` (which never
runs `BP_RESOLUTION`) and does not appear in any earlier report in `.diag/`.

---

## 5. What is contained in the tree, and what remains unproven

### 5.1 Contained (surviving changes; no commit, no branch)

| path | change |
|------|--------|
| `fastgatk-native/scripts/verify_hc_gvcf_indel_end_gatk_oracle.py` | new strict gate for target A. Currently **exit 1** (target A unfixed in the tree); turns green with the saved patch of §1.3. |
| `fastgatk-native/scripts/verify_hc_gvcf_forced_allele_refblock_gatk_oracle.py` | new strict gate for target B. Currently **exit 1**; do not register until B is fixed. |
| `.diag/round-deeplimits-gvcf.md` | this report (`.diag/` is git-ignored). |
| `.diag/scratch-deeplimits/*` | round diagnostics: `end_fix_reverted.diff` (the exact reverted patch), `probe.py`/`probe2.py`/`probe3.py` (controls), captured stdout/stderr of every run quoted above. Scratch only, no production effect. |

**Production source: unchanged.** `git status --porcelain` lists only the two new scripts; the sole
modified file, `fastgatk-native/src/hc_call.cpp`, was restored to `bb06ea4`, and the temporary
`calling_pipeline.cpp` instrumentation was reverted before the first rebuild (that file is clean).
Not touched: `fastgatk-native/CMakeLists.txt`, any root `*.md`, Mutect2, the EventMap / PairHMM / AF /
prior machinery, `verify_hc_alleles_deep_limits.py`, and every pre-existing oracle or test script.

### 5.2 Proven by this round

* Target A's exact field-level difference, GATK's rule for it (with citations), the minimal native fix,
  and the fact that the fix makes the reported case byte-identical against the pinned jar (gate 1 → 0;
  the deep-limits script's case matched).
* Target A **cannot** be landed in this round: it reddens the registered native contract test
  `fastgatk-indel-contract`, whose `verify_indel.py:173` pin is itself the non-GATK behaviour.
* Target B is a distinct defect (reference-block `PL`/`GQ` model selection), confined to the
  plausible-indel informativeness classification of the fixture's six frameshift reads under forced
  alleles — established by four controls, not by argument.
* The registered gates, `verify_hc_alleles_gatk_oracle.py`, and the full suite are green at the final
  tree state (293/293 on both backends).

### 5.3 Unproven / speculation

* **Which of H1/H2 causes target B** — not separated. H1 is the leading hypothesis on structural
  grounds (the observation-remap / `projections=0` asymmetry); H2 is **speculation** that remains open
  because the traversal was not re-derived numerically.
* **GATK's internal reference window** for the alleles run is *derived* (`[599,775]`), not observed.
  If the derivation is wrong, the localization moves to `append_segment` / `AssemblyRegionTrimmer`
  parity rather than to the read payload.
* **Whether the six frameshift reads are realigned differently by the two callers** is inferred from
  `rcm_realigned_observations` (6119 vs 5435) and does not prove which read payload GATK's
  informativeness test finally sees.
* Target B is **not fixed**, so `verify_hc_alleles_deep_limits.py` still exits 1 with
  `diverged_cases = ['max-genotype-count-2', 'gvcf-homopolymer-insertion', 'gvcf-max-alt-alleles-1']`
  — the first is the separately-tracked exit-code case (untouched), the second is fixed-but-reverted
  (§1.5), the third is this round's unfixed finding.
* The §4 `BP_RESOLUTION` reference-row divergence is a single-fixture observation; it was not narrowed
  further and is gated nowhere.
