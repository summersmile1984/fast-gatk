# Round G2 — census of the "owner first-success" pattern (read-only)

**Question.** D2 (a published gVCF allele annotated from an AssemblyRegion owner that never
requested PairHMM for that allele, so `likelihood_candidate_read_context_ordinals` holds the
`UINT32_MAX` sentinel) has now been found **twice** — once in the diploid gVCF annotation path and
once in the arbitrary-ploidy path. Enumerate every site in `fastgatk-native/src/hc_call.cpp` (and
`calling_pipeline.cpp`) that iterates `result.assembly_region_likelihood_results` or otherwise
selects among AssemblyRegion owners, and classify each: **(A)** same shape *and* can produce a
sentinel-owner-dependent wrong result, **(B)** same shape but cannot (with the reason), **(C)**
different shape / not applicable.

**Constraints honoured.** No file under `fastgatk-native/src`, `fastgatk-native/include`,
`fastgatk-kernels/*` or any `CMakeLists.txt` was edited; nothing was rebuilt; no `git` state was
changed; no root `*.md` was touched. Every experiment ran inside a single `mktemp -d` per bash
invocation. The only file created by this round is this report.

---

## 0. TL;DR

* **15 owner-selection sites** exist: **12** in `hc_call.cpp` (10 lambdas/loops + 2 already-fixed
  `ordered_owners` builders) and **3** in `calling_pipeline.cpp` (one of them the owner-pool
  plumb-through). A further ~20 analogous loops live in `mutect2_tool.cpp` — **out of scope, not
  audited here** (see §5).
* **Exactly one unfixed (A) site remains: `hc_call.cpp:3521`** — the ordinary-VCF
  side of `annotations_from_owner`. It has the identical shape as the two fixed sites (loops over
  owners, returns on the first owner whose `calculate_output_variant_annotations` succeeds, with no
  `owner_has_pairhmm_context` preference). It is *reachable in principle* but **not reachable with
  the `fixtures/chr20` fixture** (measured, see E6), so this classification is **reasoned, not
  measured**.
* The two known instances — `hc_call.cpp:4894` (diploid gVCF annotation) and `hc_call.cpp:4564`
  (arbitrary-ploidy gVCF annotation) — are **(A) but already fixed**; they keep a residual
  first-success recursion *within* each preference group (see R-A2/R-A3).
* **Every other site is (B) or (C), and the (B) classification is not an opinion**: it follows from
  a structural invariant of the native `Result` (§1), because those sites' success predicate
  requires a *finite likelihood row*, and a sentinel candidate's rows are `-inf`/`NaN` by
  construction. Only `calculate_output_variant_annotations` has a success predicate that a sentinel
  owner satisfies (index-bounds only) while silently degrading the evidence.

---

## 1. The classification rule (why some sites can be wrong and others cannot)

Three facts about `fastgatk::calling::Result`, all verified by reading the code and partially by
measurement:

**INV-1 (sentinel ⇔ no PairHMM request).** `candidate_read_context_ordinals` is initialised to
`UINT32_MAX` (`calling_pipeline.cpp:10547-10548`, `missing_context`) and is only overwritten when a
`candidate_reads` entry (a PairHMM request row) exists for that candidate
(`calling_pipeline.cpp:10556-10568`). Therefore, for a given owner `Result` **and** for the
flattened `result`:

> `likelihood_candidate_read_context_ordinals[i] == UINT32_MAX` ⟺ that `Result` holds **no PairHMM
> request** for candidate `i`.

**INV-2 (no request ⇒ non-finite rows).** Row storage is allocated as `-inf` for
`allele_read_likelihoods` / `reference_read_likelihoods` (`calling_pipeline.cpp:7211-7216`) and as
`NaN` for `non_ref_read_likelihoods` (`7218-7221`), and is written only from a request row
(`11943-11950`, and the `<NON_REF>` analogue). So INV-1 implies **every row of that candidate is
non-finite** in that `Result`.

**CONSEQUENCE-1 — row-consuming owner probes fail loudly.** A probe whose success test needs a
finite row returns `false` / `nullopt` for a sentinel owner, so a `for (owner : owners) { if
(probe(owner)) return; }` loop simply moves on to the next owner:

| probe | success test | sentinel owner |
| --- | --- | --- |
| `derive_multiallelic_depths` (`hc_call.cpp:502-571`) | `depth > 0` (finite best **and** second-best margin per read, `is.finite` gate) | `nullopt` |
| `joint_candidate_pl` (`577-652`) | `informative_reads > 0` (finite REF row **and** finite ALT row) | `false` |
| `joint_candidate_pl_with_nonref` (`717-788`) | `non_ref_reads > 0` (finite `<NON_REF>` row) | `false` |
| `derive_pairhmm_nonref_pl` (`calling_pipeline.cpp:6002-6030`) | `informative > 0` with REF, ALT **and** `<NON_REF>` all finite | `false` |

**CONSEQUENCE-2 — the annotation probe succeeds silently.** `calculate_output_variant_annotations`
(`calling_pipeline.cpp:13762-13820`) only checks that the published alleles are *present by allele
identity* in the owner and that their row **indices are in bounds** (`has_rows`, `13811-13815`). It
never inspects finiteness and never inspects the context ordinal. It then calls
`calculate_variant_annotations`, which — when it sees the sentinel ordinal — leaves
`have_context_mapping_evidence == false` (`4110-4152`) and falls back to the raw-overlap predicate as
the MQ gate (`4293-4297`). The result is a **non-nullopt annotation built from the wrong read
population**, i.e. exactly the D2 signature (`RAW_MQandDP` / `SB` change with `-L`).

> **The classification criterion used below:** a site is **(A)** iff it can return/accept an
> owner-selected value that is produced by the sentinel mechanism, i.e. iff it consumes
> `likelihood_candidate_read_context_ordinals` (directly, or through the annotation boundary). A
> site that consumes only likelihood rows is **(B)** for the sentinel mechanism.

Measured backing for CONSEQUENCE-1/-2 is in §3 (E2, E3).

**Owner order is `-L`-dependent.** The owner vector is built in partition-merge outcome order
(`calling_pipeline.cpp:15398`, `15545-15548`) and each outcome is one AssemblyRegion of the
traversal, so the number and order of owners depends on the analysis window (the D2 round measured
3 owners in window A vs 2 in window B for the same locus). That is what turns "first-success" into
`-L` dependence.

---

## 2. Site census

`hc_call.cpp` (all line numbers are the current tree; `result` = the flattened/aggregated `Result`):

| # | line | what it decides | class | reason / evidence |
| --- | --- | --- | --- | --- |
| 1 | `3134-3137` (lambda `3128`, also used at `3475` for `output_calls`) | Ordinary-VCF FORMAT AD/DP for the published ALT set (`derive_multiallelic_depths`) | **B** | Row-consuming: needs a finite, informative-margin read (CONSEQUENCE-1). Also `result` is probed **first** (`3131`), then owners. A sentinel owner returns `nullopt` → loop continues. |
| 2 | `3227-3230` (`joint_from_owner`, arbitrary ploidy, incl. `*`) | Ordinary-VCF joint PL for `sample_ploidy != 2` | **B** | Row-consuming (`joint_candidate_pl` needs `informative_reads > 0`, `has_alt` uses `isfinite`). `result` first (`3225`). |
| 3 | `3240-3243` (`joint_concrete_pl`) | Ordinary-VCF diploid concrete joint PL | **B** | Same predicate as #2. |
| 4 | **`3521-3528`** (`annotations_from_owner`) | **Ordinary-VCF annotations for the published ALT subset** (`calculate_output_variant_annotations`) | **A — unfixed** | Ordinal-consuming; success = index bounds only (CONSEQUENCE-2). No `owner_has_pairhmm_context` preference, and owners are tried **before** the `result` fallback (`3530`). Guarded by `annotation_reads != nullptr && (max_allele_subset_changed \|\| joint_af_allele_subset_changed)` (`3503-3504`) → needs a group with ≥2 concrete ALTs whose published set was reduced (`3301-3328` = `--max-alternate-alleles`; `3420-3464` = joint-AF plausibility subset). Measured unreachable in this fixture (E6). |
| 5 | `4021-4027` (`load_symbolic_pl`) | gVCF `--max-alternate-alleles` symbolic-PL load for ALT subsetting | **B** | Row-consuming (`joint_candidate_pl_with_nonref`, `<NON_REF>` rows are `NaN`-initialised). If nothing loads it throws `BACKEND_UNAVAILABLE` rather than degrading. |
| 6 | `4116-4117` (`phase_contexts`) | gVCF physical phasing: collects haplotype membership from `result` **and** every owner | **C** | Different shape: **no early return** — each owner is consumed independently and `haplotype_map` is rebuilt per context (only `physical_phases`, declared outside the loop at `4106`, persists). It consults `somatic_candidate_haplotype_indices` / `calls` / `gvcf_hom_ref_calls`, **never** context ordinals or likelihood rows → not sentinel-dependent. (The cross-context accumulation is the separately classified D1/window-dependence family, already ruled GATK-consistent by R3/Wave 1.) |
| 7 | `4347-4352` (`joint_non_ref_from_owner`, arbitrary ploidy) | Arbitrary-ploidy gVCF joint symbolic PL (+ `*` when covered) | **B** | Row-consuming (`joint_candidate_pl_with_nonref`). |
| 8 | `4535-4539` | Arbitrary-ploidy gVCF FORMAT AD/DP | **B** | Row-consuming (`derive_multiallelic_depths`, see #1). |
| 9 | `4564-4578` (fixed at commit `a942fb1`) | Arbitrary-ploidy gVCF annotations | **(A) fixed, residual** | Same shape as #4 and it *was* the live second instance (`--sample-ploidy 3` pre-fix: `97200,27 / 0,0,0,0` in windows B–E). Now `ordered_owners` puts `owner_has_pairhmm_context` owners first. Residual: still first-success **within** each group, and the non-context group can still return a sentinel-owner annotation if every context owner returns `nullopt` (R-A3). |
| 10 | `4680-4685` (`joint_non_ref_from_owner`, diploid) | Diploid gVCF joint symbolic PL | **B** | Row-consuming (as #7). |
| 11 | `4721-4725` | Diploid gVCF FORMAT AD/DP | **B** | Row-consuming (as #1). |
| 12 | `4894-4909` (fixed at commit `c577da8`) | Diploid gVCF annotations | **(A) fixed, residual** | The first D2 instance; windows `10020381/10020391/10020401/10020411` moved from `97200,27 / 0,0,0,0` to GATK's `28800,8 / 0,0,3,3`. Same residual as #9 (R-A2). |

`calling_pipeline.cpp`:

| # | line | what it decides | class | reason / evidence |
| --- | --- | --- | --- | --- |
| 13 | `6125-6129` | `gvcf_candidate_is_retained` owner fallback for the "monomorphic → fold into a reference block" decision | **B** | Row-consuming: `gvcf_candidate_is_non_monomorphic` → `derive_pairhmm_nonref_pl` requires finite REF+ALT+`<NON_REF>` (`6002-6030`) → `nullopt` on a sentinel owner → loop continues. If nothing resolves, it returns **`true` (conservative retain)**, not a degraded annotation — a different failure mode from D2. |
| 14 | `13123` / `13174-13185` | `build_profile_local_reference_blocks`: picks the owner whose read population feeds a variant-bearing RCM segment | **C** | Different shape: the owner is matched by **region-coordinate identity** (`reference_confidence_regions` tid/start/end/active_start/active_end, `13176-13184`), not by "the first owner whose call succeeds". Sentinel-ness cannot change which owner matches, so no owner-order dependence. *Watch item (unmeasured):* the per-candidate lookups **inside** the matched owner still take the first duplicate by `(tid,pos,ref,alt)`, and `gvcf_candidate_is_retained(*owner_result, …)` then degrades to the conservative `true`; that is a within-owner variant of the duplicate-candidate hazard, not an owner-order one. |
| 15 | `15398`, `15545-15548`, `15679` | Owner-pool construction (clear / reserve+`push_back` in outcome order / plumb into RCM) | **C** | No selection at all. It is what *establishes* the `-L`-dependent owner order consumed by #1–#12. |

### Shape summary

* Same shape (loop over owners, return on first success): **#1–#5, #7–#12, #13** = 12 sites.
* Of those, **ordinal-consuming** (can be fooled by the sentinel): **#4, #9, #12** = 3 sites —
  **two fixed, one open**.
* Row-consuming (cannot be fooled, by INV-1/INV-2 + CONSEQUENCE-1): **9 sites**.
* Different shape: **#6, #14, #15**.

---

## 3. Empirical work (measured, with commands)

All runs: `fastgatk-native/build/fastgatk-hc-call` (OpenMP, current tree, working tree clean at
`a942fb1`), fixture `fixtures/chr20/ref20mnp.fasta` + `fixtures/chr20/mnp.bam`,
`--max-mnp-distance 1 --threads 1 --add-output-vcf-command-line false`, one `mktemp -d` per
command, `timeout` on every invocation.

**E1 — post-fix baseline (diploid gVCF).** Windows A `20:10019901-10020710` (48 rows) and B
`20:10020381-10020710` (13 rows) both emit

```
20  10020680  .  CA  AT,*,<NON_REF>  154.25  .  DP=8;…;RAW_MQandDP=28800,8  … :0,0,3,3
```

i.e. GATK's `28800,8` / `0,0,3,3` in both windows. This reproduces the D2 fix report's acceptance
table and confirms the built binary contains the fix.

**E2 — measured proof that a sentinel owner's annotation *succeeds* with degraded evidence**
(this is CONSEQUENCE-2). `FASTGATK_DEBUG_ANNOTATION_POSITION=10020679`
(= 1-based 10020680, the divergent row), diploid gVCF:

```
=== window 20:10019901-10020710
[SUMMARY] candidate=5 … ref_f=0 ref_r=0 alt_f=3 alt_r=1 …
[SUMMARY] candidate=5 … ref_f=0 ref_r=0 alt_f=3 alt_r=3 …
=== window 20:10020381-10020710
[SUMMARY] candidate=5 … ref_f=0 ref_r=0 alt_f=0 alt_r=0 …   <-- sentinel owner: reached, all-zero evidence
[SUMMARY] candidate=1 … ref_f=0 ref_r=0 alt_f=3 alt_r=1 …   <-- identical-allele twin owner
[SUMMARY] candidate=1 … ref_f=0 ref_r=0 alt_f=3 alt_r=3 …
```

In window B the annotation boundary is entered **three** times for the same published allele: the
first pass is for candidate index 5 — the index that the D2 round's own instrumentation measured as
the sentinel duplicate in exactly this window — and its evidence is entirely empty (chain counts all
zero). The next two passes come from candidate index 1, the identical-allele twin whose ordinal is
`0`. That all-zero pass is a *silent success* — precisely the property that makes class (A)
dangerous — and it is the same window in which the pre-fix published row read `SB=0,0,0,0` /
`RAW_MQandDP=97200,27`. (In window A the only passes are for index 5, and they carry real evidence —
consistent with the D2 record that window A's owner holds `idx=5 … ord=0`.)

**E3 — the row-consuming probes really do reject the sentinel owner.** Recorded pre-fix divergence
for that row (`gvcf-stream-by-region-divergence-20260910.md`, 证据 2) lists **only** `RAW_MQandDP`
and `SB` as differing; `QUAL`, `PL`, `GT`, `GQ`, `PGT/PID/PS`, `AD` and `DP` were unchanged and
equal to GATK. If any of the row-consuming first-success probes (#1, #2, #3, #7, #8, #10, #11, #13)
had accepted the sentinel owner, AD/DP/PL for that row would have collapsed (all-zero / absent
evidence) as well. They did not. This is the measured counterpart of CONSEQUENCE-1 on the very
locus of the defect. *(Reasoned extension: the same argument transfers to the gVCF `*`/symbolic PL
probes #5/#7/#10, which additionally would have thrown `BACKEND_UNAVAILABLE` rather than degrade.)*

**E4 — ploidy sweep on the previously-failing windows (the discriminator that found instance #2).**
gVCF, windows `10020381/10020391/10020401/10020411` plus control `10019901`, at `--sample-ploidy 1`
and `3` (runs the arbitrary-ploidy annotation branch #9):

| ploidy | windows B–E, POS 10020680 | control window A | verdict |
| --- | --- | --- | --- |
| 1 | `RAW_MQandDP=28800,8`, `SB` col `0,0,3,1`, `AD=1,4,0`, `DP=5` | identical | window-invariant |
| 3 | `RAW_MQandDP=28800,8`, `SB` col `0,0,3,3`, `AD=0,4,2,0`, `DP=6` | identical | window-invariant |

and the E2-style trace at ploidy 1/3 shows the sentinel pass (all-zero) still occurring in windows
B–E while the published row keeps the twin-owner evidence. Pre-fix, that same ploidy-3 sweep read
`97200,27 / 0,0,0,0` in all four windows (recorded in the ploidy oracle's header, commit `a942fb1`).

**E5 — ordinary VCF (`vcf_text`, sites #1–#4) is window-invariant at this locus, and matches pinned
GATK.**
* POS 10020680 *is* emitted by ordinary VCF mode (no `--emit-ref-confidence`), and its row is
  **byte-identical in all five window starts** (`10019901/10020381/10020391/10020401/10020411`):
  `…DP=8;…MQ=60.00;QD=30.85;SOR=1.609  GT:AD:DP:GQ:PL  0/1:1,4:5:15:165,0,15`.
* The five windows do differ at `10020429/10020431/10020434/10020438` (`DP=31/32` vs `12/13`,
  `QD`, `FS`, `SOR`, `BaseQRankSum`, `PL`). Pinned GATK 4.6.2.0 produces **the same window-dependent
  values** (verified live for windows `10019901` and `10020381`: e.g. GATK `10020431` is
  `DP=32;QD=16.46;SOR=0.608` for start 10019901 and `DP=13;QD=17.20;SOR=0.132` for start 10020381,
  field-for-field equal to native in each window). So those differences are GATK's own
  window dependence (the already-classified D1/D3 family), **not** an owner-selection defect.
* E6 — reachability probe for the open (A) site #4: for **every** position emitted in ordinary-VCF
  mode in windows A and B, exactly **one** `[FASTGATK_ANNOTATION_SUMMARY]` pass occurs
  (`annotation_passes=1` for all 10 + 6 positions), i.e. the render-time `annotations_from_owner`
  loop at `3521` is never entered because the `3503-3504` guard is false. Consistently, **zero**
  multi-ALT records are emitted in ordinary-VCF mode for 3 ploidies (1/2/3) × 6 window starts, so
  neither `max_allele_subset_changed` nor `joint_af_allele_subset_changed` can become true here.

*Precision note on E2:* the trace names candidate **indices**, not the `Result` that supplied them,
so it cannot by itself distinguish "owner 0's sentinel duplicate" from "the flattened `result`'s
sentinel duplicate" — both are ordinal-sentinel candidates for the same allele in window B (the D2
round's instrumentation located the sentinel at `owner 0 idx=5` and the twin at `owner 1 idx=1`).
What the trace *does* establish, and all that the (A) classification needs, is: **an
ordinal-sentinel candidate for a published allele reaches the annotation boundary and returns a
non-empty, all-zero-evidence annotation instead of failing.**

**E7 — the repository's own D2 gates re-run live in this round** (each re-establishes the GATK
expectation live, prints JSON, and exits non-zero on any violation; both use only a
`TemporaryDirectory` and write nothing into the repo):

* `python3 fastgatk-native/scripts/verify_hc_window_invariance_gatk_oracle.py` (diploid,
  12 windows = 10 gating + 2 structural): `"status": "pass"`, `"strict_mode": true`,
  `"violations": []`, `ORACLE_EXIT=0`. All 10 gating windows report
  `data_rows_byte_identical=True` and `pos_10020680_identical=True` with
  `RAW_MQandDP=28800,8` / `SB=0,0,3,3` on both sides; both structural windows also report
  `pos_10020680_identical=True` with the *GATK-own* `97200,27` / `1,2,6,16`; `gatk_baseline_drift`
  is `false` everywhere.
* `python3 fastgatk-native/scripts/verify_hc_ploidy_window_invariance_gatk_oracle.py`
  (`--sample-ploidy 3`, 7 windows = 5 gating + 2 structural): `"status": "pass"`,
  `"violations": []`, `ORACLE_EXIT=0`. All five gating windows are byte-identical to live GATK
  (window A: 94 rows both sides; windows B–E: 13 rows both sides) and POS 10020680 is identical in
  **all seven** windows.

Interpretation for this census: the two fixed (A) sites are the *only* owner-selection sites whose
degradation was observable on this fixture, and they are now window-invariant. The gates cannot say
anything about `hc_call.cpp:3521` (ordinary VCF mode, single `-L` per gated window set and no
multi-ALT group) — which is exactly why that site stays classified as reasoned-only.

---

## 4. (A) sites, ranked by likely impact, with the cheapest next test

**R-A1 — `hc_call.cpp:3521` (open; highest remaining D2-shaped risk).**
*Impact:* ordinary-VCF (non-gVCF) runs in which a published ALT set was reduced — i.e. real
`--max-alternate-alleles` use, or a multi-ALT locus where joint AF dropped an allele. The affected
fields are the same family as D2: the MQ gate (`have_context_mapping_evidence` fallback) →
`MQ`/`RMSMappingQuality`-derived fields and `SOR`/`FS`, plus rank sums whenever the sentinel owner's
read population differs. A scatter by `-L` (the intended production usage) is exactly the trigger.
*Why not confirmed:* measured unreachable with `fixtures/chr20` (E6) — the guard needs a multi-ALT
group this fixture never produces.
*Cheapest next test (no code edit):* find/produce a fixture whose **ordinary-VCF** output contains a
multi-ALT group that survives a subset reduction, e.g. add `--max-alternate-alleles 1` (accepted,
`hc_call.cpp:1215-1223`) on a tri-allelic locus, or a locus where joint AF drops one of ≥2 concrete
ALTs, at a **duplicated-candidate region boundary** (the D2 precondition — i.e. the same locus must
be produced by ≥2 AssemblyRegion owners, one of them sentinel). Then run the same command with two
different `-L` starts and diff the shared rows; and/or run with
`FASTGATK_DEBUG_ANNOTATION_POSITION=<0-based pos>` — a **second** `[FASTGATK_ANNOTATION_SUMMARY]`
pass for that position means the `3521` loop ran, and an all-zero-evidence pass means it selected a
sentinel owner. Compare the same windows against pinned GATK to decide which side is right.
*Cheapest fix shape (for a future round, not applied here):* make `annotations_from_owner` use the
same `ordered_owners` preference as `4564`/`4894`, or simply call `owner_has_pairhmm_context` as a
pre-filter. Note the guard itself is orthogonal: it decides *whether* the loop runs, not *which*
owner it picks.

**R-A2 — `hc_call.cpp:4894` (fixed; residual).**
*Residual impact:* low. `ordered_owners` guarantees no sentinel-carrying owner is consulted before a
context-carrying one, so a sentinel-owner annotation can only be produced when **every** context
owner returns `nullopt` (e.g. an owner that holds the twin but not *all* published candidates of a
multi-ALT group — `calculate_output_variant_annotations` returns `nullopt` for a missing sibling,
`13786-13798`). In that configuration the old `-L` dependence can reappear.
*Cheapest next test:* re-run the E4 ploidy/window sweep on a **multi-ALT gVCF** locus where
`group.candidates` has ≥2 alleles and the group's owners are heterogeneous, asserting
`owner_has_pairhmm_context` coverage; plus a temporary-instrumentation-free check: assert that for
every published gVCF group at least one owner passes `owner_has_pairhmm_context` (a Host-side
invariant that can be checked from the debug trace).

**R-A3 — `hc_call.cpp:4564` (fixed; residual).** Same as R-A2, on the arbitrary-ploidy branch; the
`--sample-ploidy 3` oracle (E7) is the ready-made gate. Its residual test is identical: a
`sample_ploidy != 2` run whose group is multi-ALT with heterogeneous owners.

---

## 5. Reasoned vs measured (explicit)

**Measured (observed in this round, commands in §3).**
* Post-fix window invariance at ploidy 2, 1 and 3 in the previously-failing windows (E1, E4).
* A sentinel-carrying owner reaches the annotation boundary and returns an all-zero-evidence
  annotation (E2) — the "silent success" that defines class (A).
* Pre-fix divergence on that row was annotation-only (`RAW_MQandDP`, `SB`), i.e. row-consuming owner
  probes did not select the sentinel owner (E3).
* Ordinary-VCF mode: POS 10020680 is window-invariant; the `10020429-10020438` differences are
  reproduced by pinned GATK too; no multi-ALT records at any tested ploidy; one annotation pass per
  emitted position (E5, E6).
* Both repository D2 oracles re-run live in this round: **strict pass, `violations: []`, exit 0**
  (E7) — diploid (10 gating windows byte-identical) and `--sample-ploidy 3` (5 gating windows
  byte-identical, POS 10020680 identical in all 7 windows of each gate).

**Reasoned from code (read-only; not measured end-to-end).**
* INV-1/INV-2 and hence the (B) classification of every row-consuming site. This is a structural
  argument (row initialisation + write-only-from-request), corroborated but not exhaustively tested
  for every backend/option combination.
* The (A) classification of `hc_call.cpp:3521` — same shape as two proven instances, ordinal-
  consuming, no owner preference. **Its reachability is unproven and measured to be false for
  `fixtures/chr20`.**
* The residuals R-A2/R-A3 (the fallback group can still yield a sentinel-owner annotation).
* Site #14's within-owner duplicate-lookup hazard.

**Not audited (out of scope).** `mutect2_tool.cpp` contains ~20 owner loops. It does **not** call
`calculate_output_variant_annotations` (no ordinal-consuming annotation boundary of the D2 shape),
and it already contains a *different* sentinel mitigation — it scans for a non-sentinel sibling
within one `Result` (`mutect2_tool.cpp:3772-3790`) instead of switching owner. It was not classified
here; if the census should be extended, that file is the recommended G3 target (it also has
`result.assembly_region_likelihood_results[part_index]` index-paired loops, e.g. `5316-5323`, which
would need their own shape analysis).

---

## 6. Files

* This report: `.diag/round-g2-owner-pattern-census.md`.
* No other file was created, modified or deleted by this round. (Scratch lived in one
  `mktemp -d` per bash command; both oracle scripts write only to their own
  `TemporaryDirectory`.)

**Incidental observation (not caused by this round).** `git status --short` is otherwise empty, but
`fastgatk-native/scripts/verify_hc_ploidy_window_invariance_gatk_oracle.py` — the gate for the
arbitrary-ploidy D2 instance — is **untracked** (`git log -- <path>` returns nothing, while the
diploid twin was added in `9ac589d`). The file was already present and untracked when this round
started; this round only read and ran it. Whoever owns the D2 fix may want to commit it, since the
second instance currently has no committed regression guard.
