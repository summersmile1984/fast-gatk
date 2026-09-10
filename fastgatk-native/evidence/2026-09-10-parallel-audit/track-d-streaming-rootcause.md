# Track D — Root-cause localization for the `--stream-by-region` gVCF divergence

Scope: **read-only root-cause analysis**. No production code was modified; the only files
created are this report and `fastgatk-native/scripts/diagnose_stream_by_region_rootcause.py`.
No `*.md` at the repo root was touched. All scratch output went to `TemporaryDirectory`.

Repro command set (unless stated otherwise):

```bash
B=fastgatk-native/build/fastgatk-hc-call
COMMON=(-R fixtures/chr20/ref20mnp.fasta -I fixtures/chr20/mnp.bam \
        -L 20:10019901-10020710 --emit-ref-confidence GVCF --max-mnp-distance 1 \
        --threads 1 --add-output-vcf-command-line false)
```

Baseline **re-verified in this session**: `$B "${COMMON[@]}"` produces 48 records and
`diff` against pinned GATK 4.6.2.0 on the same interval is **empty (0 diff lines)**. The
non-streamed path is the correct one, so every deviation below is attributed to the
streamed side.

## Headline result

**None of D1/D2/D3 is caused by the streaming plumbing.** D2 is reproduced *exactly* by a
plain non-streamed run with a different `-L` start — no `--stream-by-region` involved at
all. D3's per-locus reference-confidence values inside a tile are degenerate
(`informative=0`), i.e. the tile's own numerics are wrong, not the stitcher. D1 lives in the
phasing block of `gvcf()`, whose input haplotype sets are assembled per tile.

The single variable that governs D2 (and explains all 8 tile sizes quantitatively) is the
`(interval_start, interval_end, assembly_region_padding)` triple that
`run_region_streaming` hands to `fastgatk::calling::run` at `hc_call.cpp:6658-6666`:

```cpp
tile_options.locus_intervals = decoded.locus_intervals;   // = the halo intervals
tile_options.intervals       = decoded.locus_intervals;
tile_options.interval_start  = tile_options.intervals.front().start;
tile_options.interval_end    = tile_options.intervals.back().end;
auto result = fastgatk::calling::run(decoded.reads, references, tile_options);
```

---

## D2 — annotation evidence scope on the spanning-deletion record (`20:10020680`)

Wrong value: `RAW_MQandDP=97200,27`, `SB=0,0,0,0`; correct (GATK and non-streamed native):
`28800,8`, `0,0,3,3`. Both have MQ=60, so the difference is purely the **evidence count**
(27 vs 8 reads entering `mapping_square_sum`/`mapping_count`).

### Discriminating experiments and observed results

**E2.1 — non-streamed `-L` start sweep (no streaming at all), pad fixed at 100.**
Holds everything constant except the read-window start; `20:10020680` is inside every window.

| `-L 20:<start>-10020710` | `10020680` RAW_MQandDP | SB |
| ---: | --- | --- |
| 10020101 | `28800,8` | `0,0,3,3` |
| 10020201 | `28800,8` | `0,0,3,3` |
| 10020301 | `28800,8` | `0,0,3,3` |
| 10020351 | `28800,8` | `0,0,3,3` |
| **10020381 … 10020411** | **`97200,27`** | **`0,0,0,0`** |
| 10020421 / 10020431 | `97200,27` | `1,2,6,16` |
| 10020451 / 10020501 / 10020601 | `28800,8` | `0,0,3,3` |

This is a **pure interval-window effect with `--stream-by-region` absent**. The switch is
non-monotone in window size (a 10020381 start is "wrong") and there is a third, distinct
wrong state (`SB=1,2,6,16`) for starts 10020421–10020431, so this is not a
"more reads → more evidence" effect.

**E2.2 — tile sweep, and 8/8 agreement with the halo-window rule.**
The tile containing `10020680` has core start `c` and halo start `c − context`, with
`context = max(padding=100, 100, indel_padding)` (`hc_call.cpp:6457-6462`):

| tile | containing core | halo start | observed | matches `-L` sweep prediction |
| ---: | ---: | ---: | --- | --- |
| 100 | 10020601 | 10020501 | `28800,8` correct | yes |
| 200 | 10020501 | 10020401 | `97200,27` wrong | yes |
| 300 | 10020501 | 10020401 | `97200,27` wrong | yes |
| 405 | 10020306 | 10020206 | `28800,8` correct | yes |
| 500 | 10020401 | 10020301 | `28800,8` correct | yes |
| 600 | 10020501 | 10020401 | `97200,27` wrong | yes |
| 700 | 10020601 | 10020501 | `28800,8` correct | yes |
| 810 | 10019901 | 10019801 | `28800,8` correct | yes |

The streamed value at the locus equals the non-streamed value for the **same read-window
start**, for all 8 tile sizes.

**E2.3 — padding sweep at fixed tile (500):** pad 100/150 → correct, pad 200/300 → wrong.
Not a counterexample but a *confirmation of the same mechanism*: raising padding changes both
the halo window (via `context`) and `tile_options.assembly_region_padding` inside
`calling::run`, so it moves the same argument triple. (At pad=200 the tile-500 halo start
10020201 coincides with a *correct* non-streamed window, but the non-streamed comparator
uses pad=100; the extra variable is the per-tile `assembly_region_padding`, which changes
AssemblyRegion partitioning inside the tile.)

**E2.4 — the deviation is confined to the `*` (spanning-deletion) record.**
At the D2 locus only the record with `*` in ALT changes; `10020679 AC>TA,<NON_REF>` keeps
`RAW_MQandDP=25200,7 SB=3,1,0,2` across all window starts. With
`--max-mnp-distance 0` the `*` record disappears entirely and the window-start sensitivity
vanishes: both `-L 10019901` and `-L 10020401` give `10020680 C>A,<NON_REF> … 25200,7`.

**E2.5 — `FASTGATK_DEBUG_ANNOTATION_POSITION` trace.**
With the gate on the candidate's 0-based position (10020679) the same 6 reads appear in
**both** modes (`20GAVAAXX100126:8:44:10223:188523`, `…:6:3:6181:93202`, `…:1:62:21171:66429`,
`20FUKAAXX100202:4:43:2630:81606`, `…:3:45:1441:112187`, `…:1:8:20035:107624`), so the
printed evidence rows are *not* what changes — the realignment-context population is. The
`[FASTGATK_ANNOTATION_SUMMARY]` lines for position 10020679 are:
non-streamed → `(ref_f=0,ref_r=0,alt_f=3,alt_r=1)` then `(0,0,3,3)` (the emitted `*` record
takes the second); tile=300 → `(0,0,0,0)`, `(0,0,3,1)`, `(0,0,0,0)` (the emitted `*` record
takes the all-zero one).

### Responsible code location and mechanism

1. The `*` record's annotations are **recomputed at render time**, in `gvcf()`, by the branch
   guarded on `(include_spanning_deletion || group.max_alt_subset) && annotation_reads != nullptr`:
   * `fastgatk-native/src/hc_call.cpp:4513-4531` (candidate-site class) and
   * `fastgatk-native/src/hc_call.cpp:4830-4848` (the second candidate-site class),
   both calling `calculate_output_variant_annotations(*annotation_reads, *owner, group.candidates, …)`
   over `result` and each `result.assembly_region_likelihood_results` owner.
   `annotation_reads` is non-null in both modes (`hc_call.cpp:6679` for streamed,
   `hc_call.cpp:6986` for non-streamed), so the recompute itself is not stream-specific.
2. The evidence population is built in
   `fastgatk-native/src/calling_pipeline.cpp:4109-4152` (`context_mapping_evidence`): for each
   likelihood index it walks the candidate's `likelihood_candidate_read_context_ordinals` →
   `likelihood_candidate_read_realignments` entries, keeps those with `context->qualified`
   whose *realigned* interval overlaps the event, and marks those physical records. `mapping_count`
   = number of marked records → **27 instead of 8**.
3. `calculate_output_variant_annotations` (`calling_pipeline.cpp:13766-13831`) resolves the
   candidate to `all_indices`/`published_indices` inside a given `Result` and calls
   `calculate_variant_annotations` with `spanning_deletion_is_output=true`; the BestAllele
   classification (`calling_pipeline.cpp:4235-4270`) decides the strand counts. In the wrong
   state every read fails `likelihood_informative && likelihood_best_allele_is_output`, so all
   four strand counters stay 0 (→ `SB=0,0,0,0`) while `context_mapping_evidence` still counts 27.

**Mechanism (as supported):** the candidate's AssemblyRegion read-realignment context is a
function of the interval span handed to `calling::run` (`hc_call.cpp:6658-6665`), because the
comment there explicitly notes that constraining `intervals` "changes its parent span and
therefore its AssemblyRegion partitioning". A shifted `-L` start (equivalently, a shifted tile
halo) re-partitions the regions, which changes which reads land in the candidate's context
ordinal and therefore the D2 evidence count and BestAllele outcome. The streaming path is a
*trigger*, not the defect.

### What remains unproven, and the confirming experiment

Unproven: **which** of the two inputs (interval span vs `assembly_region_padding`) is
dominant when they disagree (the pad=200 tile-500 case), and the exact reason the context
population jumps 8 → 27.

Next single experiment (no rebuild required): run the annotation audit on the two adjacent
non-streamed windows that straddle the E2.1 switch, i.e.
`FASTGATK_DEBUG_ANNOTATION_POSITION=10020679` with `-L 20:10020351-10020710` (correct) versus
`-L 20:10020381-10020710` (wrong), and compare the `[FASTGATK_ANNOTATION_SUMMARY]` line set
(`candidate=`, `grouped_candidates=`, `output_candidates=`). If a second candidate or a larger
`output_candidates` appears at the switch, the mechanism is candidate grouping/RegionEventMap
membership at the D2 locus; if the counts are unchanged, the mechanism is the context-ordinal →
realignment mapping in `calling_pipeline.cpp:4109-4152`. The invocation sequence is already
implemented in the diagnostic script (`--window-starts` sweep + the phase/RCM traces).

---

## D1 — spurious phasing at `20:10020228` / `20:10020229`

Wrong: `GT=0|1`/`1|0` plus `PGT=10020228_G_A`, `PID`, `PS`. Correct (GATK and non-streamed
native): unphased `0/1`, no `PGT/PID/PS`.

### Discriminating experiments and observed results

**E1.1 — tile sweep, D1 column (from the diagnostic script):**

| tile | `10020228` GT | D1 present |
| ---: | --- | --- |
| 100 | `0|1` | yes |
| 200 | `0|1` | yes |
| 300 | `0/1` | no |
| 405 | `0|1` | yes |
| 500 | `0|1` | yes |
| 600 | `0/1` | no |
| 700 | `0/1` | no |
| 810 | `0/1` | no |

D1 is non-monotone in tile size and appears at tiles that are *larger* than others where it is
absent (405/500 have it, 600/700 do not) — so it is not a size/batch effect.

**E1.2 — `FASTGATK_DEBUG_PHASE=1` trace (`hc_call.cpp:4069, 4140, 4236`).** This is the
decisive evidence. Per-group haplotype sets fed to the phasing block:

| mode | `10020228` haps (per context, in order) | `10020229` haps | `PHASE_OUTPUT` for 10020228 |
| --- | --- | --- | --- |
| whole | `''`, then `0,1,` | `''`, then `2,3,` | **absent** |
| tile=300 | `0,1,` (single context) | `2,3,` | **absent** |
| tile=500 | `''`, `0,`, `''` (three contexts) | `''`, `1,`, `''` | **present (`0|1 pid=10020228_G_A`)** |

`[FASTGATK_PHASE_INPUT]` also differs: whole = 1 call (10 candidates/groups, 16 events, 8
outputs); tile=500 = 2 calls (4 and 6 candidates/groups, 18 events, 10 outputs).

**E1.3 — decoupling from D2.** pad sweep at fixed tile=500: D1 stays `0|1` for pad 100/150/200/300
while D2 flips at pad 200. Tiles 405/500/100 have D1 present with D2 correct; tiles 300/600 have
D1 absent with D2 wrong.

### Responsible code location and mechanism

`fastgatk-native/src/hc_call.cpp:4085-4195`, inside `gvcf()`. The block is the only place that
emits `GT` with `|` plus `PGT/PID/PS`, so the divergence provably enters here. Its inputs are:

* `haplotype_map[group_index]` — built by iterating **`phase_contexts` = `{&result}` ∪
  `{result.assembly_region_likelihood_results[*]}` and doing
  `haplotype_map[group_index].insert(members.begin(), members.end())`
  (`hc_call.cpp` ~4135-4139). Each owner contributes its own **locally-numbered**
  `somatic_candidate_haplotype_indices`, so haplotype ordinals from different owners are
  unions in a single namespace.
* `total_available_haplotypes` = the union over **all** groups (`hc_call.cpp` ~4143-4146).
* The decision predicates `always_together` / `always_apart` (`hc_call.cpp` ~4150-4195) compare
  set sizes against `total_available_haplotypes`, and `always_apart` requires
  `|A| + |B| == total_available_haplotypes` with disjoint sets.

Because the number of phase contexts and the per-owner haplotype ordinals are a function of how
the tile's `calling::run` partitioned the interval into AssemblyRegions, the *same* per-locus
genotype evidence yields a different phasing decision per chunk. In whole/tile=300 the phasing
set `{0,1}` vs `{2,3}` does not satisfy the size test; in tile=500 the sets become `{0}` vs `{1}`
while the total shrinks, and the block declares them phased.

### What remains unproven, and the confirming experiment

Unproven: the exact predicate that flips (I did not reproduce the grouping-loop arithmetic from
the observed sets; the observed sets alone do not let me replay `phase_group`/`phase_01`
without printing `total_available_haplotypes` and the intermediate `phase_group`). Also unproven:
that merging owners' haplotype ordinals is the *only* cause — the differing `[PHASE_INPUT]`
call counts show candidates/groups also differ.

Next single experiment: print `total_available_haplotypes` and the per-group decision inside the
loop (a temporary `FASTGATK_DEBUG_PHASE` line — requires a build, which this track is not allowed
to do), for whole vs tile=500 vs tile=300. Alternatively, a rebuild-free partial check: confirm
that D1 tracks the *number of AssemblyRegion owners* rather than the tile boundary position by
comparing tiles 405 and 600 (both have their D1 locus in a mid-interval core but differ in
whether D1 appears).

---

## D3 — reference-block granularity (`20:10020230-10020428`)

Correct (GATK and non-streamed native, byte-identical): 29 small blocks over that span, e.g.
`10020230 0/0:12:36:12:0,36,494`, `10020231-10020233 0/0:13:39:13:0,39,526`, …,
`10020427-10020428 0/0:31:90:30:0,90,1250`. Streamed tile=500: a single
`10020230 <NON_REF> END=10020428 GT:DP:GQ:MIN_DP:PL = 0/0:26:0:21:0,0,0`.

### Discriminating experiments and observed results

**E3.1 — block count in the span vs a tile core boundary falling inside it.**
Cores are `[interval_start + k·tile, +tile]`, independent of padding.

| tile | core boundaries inside `(10020230, 10020428)` | blocks in span |
| ---: | --- | ---: |
| 100 | 10020301, 10020401 | 1 |
| 200 | 10020301 | 22 |
| 300 | none (10020201, 10020501) | **29 (correct)** |
| 405 | 10020306 | 22 |
| 500 | 10020401 | 1 |
| 600 | none (10020501) | **29** |
| 700 | none (10020601) | **29** |
| 810 | none | **29** |

Perfect correlation: D3 appears exactly when a tile core boundary strictly inside the span. With
tile=200 the degenerate block is exactly the tail of the earlier core
(`10020230 <NON_REF> END=10020300 0/0:24:0:21:0,0,0`) and every block after the boundary is
healthy — so the degeneracy is the **terminal part of a tile's core**, not the span as a whole.

**E3.2 — `FASTGATK_DEBUG_RCM=1` per-locus reference confidence (the decisive measurement).**

| position | whole (non-streamed = GATK) | tile=500 |
| --- | --- | --- |
| 10020230 | `depth=12 informative=12 gq=36 pl=0,36,494` | `depth=8 informative=0 gq=0 pl=0,0,0` |
| 10020231 | `depth=13 informative=13 gq=39 pl=0,39,562` | `depth=8 informative=0 gq=0 pl=0,0,0` |
| 10020401 | `depth=26 informative=22 gq=66 pl=0,66,990` | `depth=26 informative=0 gq=0 pl=0,0,0` |
| 10020428 | `depth=31 informative=30 gq=90 pl=0,90,1350` | `depth=30 informative=0 gq=0 pl=0,0,0` |
| 10020680 (D2 locus) | `6 6 0 …` / `25 0 0 …` | identical (SAME) |

1608 of 7459 joined loci have a different `(depth, informative, gq, pl)` tuple between whole and
tile=500. The tile genuinely computes `informative=0` — and a lower depth (8 vs 12) — for these
loci, producing `gq=0 / pl=0,0,0`, which then collapses into one block. `10020680`'s RCM values are
identical in both modes, confirming D2 and D3 are in different subsystems (candidate-site
annotation vs reference confidence).

**E3.3 — `GQ=0 / PL=0,0,0` is not by itself a defect.** Pinned GATK 4.6.2.0 emits three such long
blocks on this very interval — `10019970-10020227 0/0:24:0:18:0,0,0`,
`10020439-10020678 0/0:26:0:20:0,0,0`, `10020682-10020710 0/0:26:0:25:0,0,0` — and native
non-streamed reproduces the full 39-line `<NON_REF>` block list identically. The defect is the
*placement* of the degenerate values (they swallow 29 healthy blocks), not their existence.

### Responsible code location and mechanism

* The value is produced by the per-tile reference-confidence build
  `fastgatk-native/src/calling_pipeline.cpp:13111` `build_profile_local_reference_blocks(...)`
  (invoked on `activity.profile_regions` + `active_reads` from `calling_pipeline.cpp:15206` when
  `calling_regions.empty()`, and from `calling_pipeline.cpp:15672` on the partition-merge path).
  Since tile=500's degeneracy is already present per locus with `informative=0`, the defect is in
  the observation set this builder receives for the **tail of the tile's core**, i.e. in the
  tile's `profile_regions`/`active_reads`, which are derived from the interval span set at
  `hc_call.cpp:6658-6665` (the same triple as D2).
* The structural amplification — 29 blocks → 1 — is `restrict_stream_gvcf_to_core()`
  (`hc_call.cpp:6668`) truncating the tile's blocks at the core end, followed by
  `RegionGvcfStitcher::merge()` (`hc_call.cpp:5864` per the prior evidence doc; `pl` element-wise
  min), which joins the degenerate tail block with the next tile's healthy block. Confirmed to be
  an amplifier only: the degenerate `gq=0/pl=0,0,0` values exist before merging.

### What remains unproven, and the confirming experiment

Unproven: (a) which of the two builders actually emitted the degenerate pass (both call sites are
live; the per-pass `FASTGATK_DEBUG_RCM` output does not tag the builder on the `[FASTGATK_RCM]`
lines), and (b) whether the D3 trigger is the interval span or `assembly_region_padding` (same
ambiguity as D2). Also unproven: that `-L 20:10020301-10020710` non-streamed reproduces the
collapse *independently* — that observation is **confounded**, because that `-L` start also
removes the earlier records from the output interval; every confound-free start (≤ 10020230)
gives 29 blocks.

Next single experiment: run `FASTGATK_DEBUG_RCM=1` for whole vs tile=500 vs tile=300 and compare
(a) the number of `[FASTGATK_RCM_BUILD_INPUT] profiles=N` pass headers per mode, which identifies
which builder produced the degenerate pass, and (b) the profile-region extents in that header for
the pass containing position 10020300. This distinguishes "the tile's ActivityProfile/profile
regions stop short of the core end" (profile-region defect) from "the builder drops observations
inside a correctly-sized profile" (observation-set defect) without touching any source file.

---

## Hypotheses REFUTED (with the evidence that refuted them)

1. **"D2 is caused by the streaming plumbing (per-tile `calling::run`, core restriction, stitcher,
   or batch size)."** REFUTED. A non-streamed run with `-L 20:10020381-10020710` (and 10020391,
   10020400, 10020401, …, 10020431) produces the identical wrong tuple `RAW_MQandDP=97200,27`,
   `SB=0,0,0,0` with no `--stream-by-region` at all (E2.1). Additionally tile=810's batch is larger
   than tile=300's yet gives the correct value (E2.2), so "batch too large" is wrong.
2. **"Raising `--assembly-region-padding` repairs or explains D2."** REFUTED as an explanation.
   pad=100/150 → correct, pad=200/300 → wrong at fixed tile=500 (E2.3). Padding is a second way to
   move the same `(interval, padding)` argument pair, not a fix.
3. **"D1 and D2 are the same root cause."** REFUTED. Tile=500 (and 405, 100) has D1 present with D2
   correct; tile=300 (and 600) has D1 absent with D2 wrong (E1.3, E2.2). They also live in
   different subsystems: D1 only changes `GT`/`PGT`/`PID`/`PS` after genotyping, while D2 only
   changes INFO/FORMAT annotations of the `*` record whose GT/PL/QUAL are unchanged; and D2's
   per-locus RCM values are identical across modes while D1's phasing differs (E3.2).
4. **"The degenerate `PL=0,0,0` block is created by `RegionGvcfStitcher::merge()`."** REFUTED at
   the per-locus level. `FASTGATK_DEBUG_RCM` shows the tile itself computes
   `informative=0 gq=0 pl=0,0,0` at 10020230/10020231/10020401/10020428 (E3.2). `merge()` only
   joins already-degenerate values — it amplifies the *structure* (29 blocks → 1), it does not
   create the values.
5. **"A `GQ=0 / PL=0,0,0` reference block is inherently wrong output."** REFUTED. GATK 4.6.2.0
   itself emits `10019970-10020227`, `10020439-10020678`, `10020682-10020710` as
   `GQ=0 / PL=0,0,0` on this interval, and native non-streamed matches all 39 `<NON_REF>` lines
   exactly (E3.3). So the D3 defect is the *placement* of degenerate values, not their existence.
6. **"The `--stream-by-region` value is reproduced by the streamed code path specifically."**
   REFUTED for D2, unproven for D3. D2's wrong value is bit-identical to a non-streamed
   `-L`-shifted run (E2.1). For D3, the non-streamed `-L 10020301` collapse is **confounded**
   (that start also drops the earlier records from the output interval); every confound-free start
   gives 29 blocks, so D3's non-streamed reproduction remains unproven — what IS proven is the
   core-boundary correlation (E3.1) and the per-locus `informative=0` degeneracy (E3.2).

## Honest incompleteness

* All three defects are narrowed to specific functions and mechanisms, but **none is proven to the
  line of code that computes the wrong number**: D2's 8 → 27 evidence jump and D1's predicate flip
  both need one instrumentation line, which this track was not permitted to build.
* The D2 "which input dominates" question (interval span vs `assembly_region_padding`) is open.
* Behaviour outside `20:10019901-10020710` on `fixtures/chr20/mnp.bam`, and combinations with
  `--stream-by-contig` / `--max-alternate-alleles` / `--read-filter` options, were not exercised.
* D3's degenerate-block mechanism was verified on tile=200 and tile=500 only.

## Reproducible script

`fastgatk-native/scripts/diagnose_stream_by_region_rootcause.py` — diagnostic (not a pass/fail
oracle). Runs the tile sweep, the non-streamed `-L` window sweep, the `--max-mnp-distance 0`
control, and captures the `FASTGATK_DEBUG_PHASE` traces; prints a JSON summary ending with
`"status": "diagnostic-complete"`. All scratch under `TemporaryDirectory`; honours
`FASTGATK_HC_BINARY` / `--native`. Runtime ≈ 30 s on this fixture set. Not registered in CTest.
