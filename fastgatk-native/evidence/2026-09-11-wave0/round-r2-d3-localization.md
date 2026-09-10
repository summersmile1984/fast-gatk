# R2 — Localization of D3 (tile-local reference-confidence degeneracy), `--stream-by-region`

Scope: **read-only**. No production source was edited, nothing was rebuilt, no `*.md` at the repo root was
touched. The only file created by me is this report plus scratch evidence under `.diag/r2-*`.
All experiments were capped with `timeout 600` and run inside single `bash` invocations.

Repro base:

```bash
B=fastgatk-native/build/fastgatk-hc-call
COMMON=(-R fixtures/chr20/ref20mnp.fasta -I fixtures/chr20/mnp.bam \
        -L 20:10019901-10020710 --emit-ref-confidence GVCF --max-mnp-distance 1 \
        --threads 1 --add-output-vcf-command-line false)
```

---

## 0. Result in one paragraph

The tile-local reference confidence is degenerate because **`--stream-by-region` replaces the user's `-L`
interval with the tile's *halo* interval when it calls `fastgatk::calling::run`**
(`fastgatk-native/src/hc_call.cpp:6658-6665`, halo computed by `input_halo_intervals` at
`hc_call.cpp:6459-6480`). That different interval changes the *ActivityProfile → AssemblyRegion
partition* inside `run`, which changes the **active span of the EventMap group** that owns the locus.
`build_profile_local_reference_blocks` (`calling_pipeline.cpp:13111`) resolves the RCM evidence owner
per *segment* by requiring the segment to be **contained in a group's active span**
(`calling_pipeline.cpp:13527-13532`); every position between the two shortened groups has no owner, so it
takes the **no-variation flank branch** (`calling_pipeline.cpp:13534-13540`), which never receives the
PairHMM realignment projection (`calling_pipeline.cpp:13173` guard + `13334-13440`). Without the
projection the host indel-informativeness test counts **0 informative reads**
(`calling_pipeline.cpp:3150-3180`), so the all-zero indel model wins `getGLwithWorstGQ`
(`calling_pipeline.cpp:3184-3195`) → `GQ=0 / PL=0,0,0` → the block combiner merges 29 healthy blocks into 1.

**The decisive observation**: a **non-streamed** run with `-L 20:10019901-10020500` (exactly the halo
window of the tile that owns 20:10020230 at `--stream-by-region 500`) reproduces the streamed tile's rows
**byte-for-byte**, including the `10020230 <NON_REF> … 0/0:26:0:21:0,0,0` block, and produces the identical
`activity.profile_regions` list. D3 is therefore a pure function of the interval window handed to
`calling::run` — not of the stitcher, not of a different builder, not of the read batch size.

---

## 1. Commands run and observed traces

### 1.1 Tool: `FASTGATK_DEBUG_RCM=1` now tags which pass produced which locus

`FASTGATK_DEBUG_RCM=1` emits four kinds of lines (`grep -n FASTGATK_RCM fastgatk-native/src/calling_pipeline.cpp`):

| line | site | meaning |
| --- | --- | --- |
| `[FASTGATK_RCM_BUILD_INPUT] profiles=N event_trims=M` | `calling_pipeline.cpp:13135` | start of one `build_profile_local_reference_blocks` call |
| `[FASTGATK_RCM_PROFILE] tid=0 span=A-B event_trims=K` | `calling_pipeline.cpp:13498/13505` | one `activity.profile_regions` entry (`A`/`B` are 0-based inclusive) |
| `[FASTGATK_RCM_EVENT_TRIM] padded=… event=…` | `calling_pipeline.cpp:7847` (inside `run_pairhmm`) | one EventMap **group** after `AssemblyRegionTrimmer`: `padded` = trimmed region, `event` = the group's **active span** |
| `[FASTGATK_RCM] tid=0 pos=P depth=D informative=I gq=G pl=…` | `calling_pipeline.cpp:3283` (inside `calculate_reference_confidence`) | one locus of one computed pass |

Note: `[FASTGATK_RCM]` lines appear **before** the first `BUILD_INPUT` header. Those come from the
per-AssemblyRegion recursive child runs (`calling_pipeline.cpp:15364` calls `run` again with
`region_options.assembly_region_independent_pass = true`, `calling_pipeline.cpp:15313`, so the children take
the `build_reference_blocks` branch at `16507`/`17633`, which prints no header). Only the pass after the
header is the one that ends up in `result.reference_blocks`, because the merge path replaces them
(`calling_pipeline.cpp:15671-15681`, `build_profile_local_reference_blocks` starts with `blocks.clear()`).

### 1.2 Pass / profile / group trace, whole vs tile=500 vs tile=300

```bash
FASTGATK_DEBUG_RCM=1 timeout 600 $B "${COMMON[@]}"                 -O whole.g.vcf   # 48 records
FASTGATK_DEBUG_RCM=1 timeout 600 $B "${COMMON[@]}" --stream-by-region 500 -O t500.g.vcf  # 20 records
FASTGATK_DEBUG_RCM=1 timeout 600 $B "${COMMON[@]}" --stream-by-region 300 -O t300.g.vcf  # 48 records
```

Profile lists (0-based, as printed) and group active spans:

| mode | `activity.profile_regions` (the pass after `BUILD_INPUT`) | EventMap group `event=` spans |
| --- | --- | --- |
| whole (= GATK) | `10019900-10020198`, **`10020199-10020488`**, `10020489-10020709` | `10019966-10019968`, **`10020227-10020437`**, `10020678-10020680` |
| tile=500 tile#1 (0-based core `[10019900,10020400)`) | `10019900-10020198`, **`10020199-10020283`**, **`10020284-10020499`** | `10020227-10020228`, `10020428-10020437`  ← **group split** |
| tile=500 tile#2 (0-based core `[10020400,10020710)`) | `10020300-10020488`, `10020489-10020709` | `10020428-10020437` (only) |
| tile=300 tile#2 (0-based core `[10020200,10020500)`) | `10019900-10020198`, **`10020199-10020488`**, `10020489-10020599` | **`10020227-10020437`** (unsplit) |
| tile=300 tile#3 (0-based core `[10020500,10020710)`) | `10020400-10020657`, `10020658-10020709` | `10020428-10020437` |

### 1.3 Per-locus RCM at 20:10020230 (1-based; = 0-based 10020229)

| mode | all `[FASTGATK_RCM] pos=10020230` lines in order | healthy line present? |
| --- | --- | --- |
| whole | `depth=23 informative=0`, `depth=23 informative=0`, **`depth=12 informative=12 gq=36 pl=0,36,494`** | yes (last pass) |
| tile=300 | `23/0` ×5, **`12/12 36 0,36,494`** | yes |
| tile=405 | `23/0` ×4, **`12/12 36 0,36,494`** | yes |
| tile=500 | `23/0` ×4, `depth=8 informative=0` | **no — every pass degenerate** |
| tile=100 | `15/0`, `23/0` ×9, `depth=8 informative=0` | **no** |

Emitted non-reference rows in `[10020200,10020440]`:

| tile | degenerate block(s) | healthy rows resume at | in-span non-ref rows |
| ---: | --- | ---: | ---: |
| 100 | `10020230-10020428 0/0:26:0:21:0,0,0` | 10020430 | 6 |
| 200 | `10020230-10020300 0/0:24:0:21:0,0,0` | 10020301 | 27 |
| **300** | **none** | — | 34 (identical to whole) |
| 405 | `10020230-10020305 0/0:24:0:21:0,0,0` | 10020306 | 27 |
| 500 | `10020230-10020428 0/0:26:0:21:0,0,0` | 10020430 | 6 |
| 600/700/810 | none | — | 34 (identical to whole) |

**The degenerate span always starts at 10020230 and ends either at the owning tile's core boundary
(200 → 10020300, 405 → 10020305) or where the next non-degenerate group begins (500 → 10020428).**
The start is *not* a tile boundary: 10020230 is the first 1-based position the tile's shortened group does
not cover (see §2). With tile=100, several consecutive 100 bp cores are all degenerate and the stitcher
merges them into one block, which is why the block extends past its own core end (10020428, same
DP/MIN_DP as the tile=500 block).

### 1.4 The decisive experiment — non-streamed reproduction from the halo window

Tile geometry: `context = max(assembly_region_padding=100, max_probability_propagation_distance+50=100,
indel_padding_for_genotyping=75) = 100` (`hc_call.cpp:6459-6463`); halo = `core ± context` intersected with
`-L` (`hc_call.cpp:6465-6480`). Cores are `[10019900 + k·tile, +tile)` (0-based), so for
`--stream-by-region 500` the tile that owns 10020230 has core `[10019900,10020400)` and halo
**`[10019900,10020500)` = `-L 20:10019901-10020500`**.

```bash
FASTGATK_DEBUG_RCM=1 timeout 600 $B -R fixtures/chr20/ref20mnp.fasta -I fixtures/chr20/mnp.bam \
  -L 20:10019901-10020500 --emit-ref-confidence GVCF --max-mnp-distance 1 --threads 1 \
  --add-output-vcf-command-line false -O halo.g.vcf          # NO --stream-by-region at all
```

Observed:

* `activity.profile_regions` = `10019900-10020198`, `10020199-10020283`, `10020284-10020499` —
  **identical to the tile#1 profile list** of `--stream-by-region 500` (§1.2).
* Emitted rows with `POS ∈ [10020201,10020400]` are **byte-identical, as complete VCF lines**, to the
  streamed tile's rows (3 rows each; Python `a == b` over the raw lines is `True`):

```text
H/S 20  10020228  .  G  A,<NON_REF>  140.64  .  BaseQRankSum=1.834;DP=6;ExcessHet=0.0000;MLEAC=1,0;MLEAF=0.500,0.00;MQRankSum=0.000;RAW_MQandDP=21600,6;ReadPosRankSum=0.842  GT:AD:DP:GQ:PGT:PID:PL:PS:SB  0|1:2,4,0:6:57:0|1:10020228_G_A:148,0,57,154,69,223:10020228:2,0,3,1
H/S 20  10020229  .  T  G,<NON_REF>   16.63  .  BaseQRankSum=0.431;DP=7;ExcessHet=0.0000;MLEAC=1,0;MLEAF=0.500,0.00;MQRankSum=0.000;RAW_MQandDP=25200,7;ReadPosRankSum=-0.210  GT:AD:DP:GQ:PGT:PID:PL:PS:SB  1|0:5,1,0:6:24:1|0:10020228_G_A:24,0,202,39,205,243:10020228:4,1,1,0
H/S 20  10020230  .  A  <NON_REF>       .    .  END=10020428  GT:DP:GQ:MIN_DP:PL  0/0:26:0:21:0,0,0
```

(The same run also reproduces the D1 spurious phasing at 10020228/10020229 — noted, not my task.)

### 1.5 Window sweep: the split is a pure function of the traversal window

Non-streamed, `-L 20:10019901-<end>`, varying only the window end (start fixed at the `-L` start):

| window end | `activity.profile_regions` | RCM @10020230 | emitted non-ref rows in span |
| ---: | --- | --- | --- |
| 10020488 … 10020499 | `…, 10020199-<end-1>` | **12/12 36 0,36,494** | healthy |
| **10020500 … 10020529** | `…, 10020199-10020283, 10020284-<end-1>` | degenerate (all lines `informative=0`) | **degenerate block 10020230-10020428** |
| 10020530 | `…, 10020199-10020494, 10020495-10020529` | **12/12 36 0,36,494** | healthy |
| 10020540 … 10020710 | `…, 10020199-10020488, 10020489-<end-1>` | **12/12 36 0,36,494** | healthy (identical to whole) |

Tile halo windows (0-based halo `[core_start-100, core_end+100)` expressed as `-L`, i.e. 1-based
`[start+1, end]`) vs the table above — the `-L` ends I ran are exact where noted:

| tile | owning core (0-based) | tile halo as `-L` | `-L` run | profiles | group split? | D3? |
| ---: | --- | --- | --- | --- | --- | --- |
| 100 | `[10020200,10020300)` | `20:10020101-10020400` | `…-10020401` (+1 bp) | `…, 10020199-10020400` | yes (tail candidates out of window) | yes |
| 200 | `[10020100,10020300)` | `20:10020001-10020400` | `…-10020400` **exact** | `…, 10020199-10020399` | yes (`event=10020227-10020228`) | yes (core ends 10020300) |
| **300** | `[10020200,10020500)` | `20:10020101-10020600` | `…-10020601` (+1 bp) | `…, 10020199-10020488, 10020489-10020600` | **no** (`event=10020227-10020437`) | no |
| 405 | `[10019900,10020306)` | `20:10019901-10020405` | `…-10020406` (+1 bp) | `…, 10020199-10020405` | yes | yes (core ends 10020305) |
| **500** | `[10019900,10020400)` | `20:10019901-10020500` | `…-10020500` **exact** | `…, 10020199-10020283, 10020284-10020499` | yes (see §4) | yes |
| **600** | `[10019900,10020500)` | `20:10019901-10020600` | `…-10020601` (+1 bp) | `…, 10020199-10020488, 10020489-10020600` | **no** | no |
| 700 | `[10019900,10020600)` | `20:10019901-10020700` | `…-10020701` (+1 bp) | `…, 10020199-10020488, 10020489-10020700` | no | no |
| 810 | `[10019900,10020710)` | `20:10019901-10020710` | `…-10020710` **exact** (= full `-L`) | identical to whole | no | no |

This is exactly the E3.1 tile/BLOCKS correlation, re-derived and now explained.

**Strength of the reproduction (honest statement):**

* tile=500: **exact** — the halo run's rows in the core range are byte-identical full VCF lines to the
  streamed rows, including the block `10020230 … END=10020428 0/0:26:0:21:0,0,0` (§1.4).
* tile=200: agreement is **per locus**, not per merged block. The exact halo run (verified) yields
  `depth=23 informative=0 gq=0 pl=0,0,0` at 10020230/10020231 and `27/0` at 10020300 — identical to the
  tile's own degenerate passes — but its emitted block is `10020230-10020400 0/0:25:0:21:0,0,0` while the
  streamed output shows `10020230-10020300 0/0:24:0:21:0,0,0`, because the streamed block is truncated at
  the tile's core end (10020300) and the remainder of the span is covered by the *next* tile, whose window
  is not split. The per-locus degeneracy (not the merged statistics) is what the halo window determines.

---

## 2. Mechanism (code path, line by line)

1. `hc_call.cpp:6459-6464` computes `context = 100`; `hc_call.cpp:6465-6480` (`input_halo_intervals`)
   builds `halo = core ± context` intersected with the user's `-L` intervals.
2. `hc_call.cpp:6658-6665` assigns **that halo window** to the tile's call options:

   ```cpp
   tile_options.locus_intervals = decoded.locus_intervals;   // the halo intervals
   tile_options.intervals       = decoded.locus_intervals;
   tile_options.interval_start  = tile_options.intervals.front().start;   // 0-based halo start
   tile_options.interval_end    = tile_options.intervals.back().end;
   auto result = fastgatk::calling::run(decoded.reads, references, tile_options);   // 6666
   ```

3. Inside `run`, those bounds become the ActivityProfile traversal span
   (`calling_pipeline.cpp:14940-14950` `activity_options`), and the traversal span bounds the profile's
   zero shoulder/region materialization (`fastgatk-kernels/src/activity_profile.cpp:1104-1157`,
   `append_zero_profile_regions`, and the final `flush_profile(tid, true)`).
   The region segmentation itself is `BandPassProfile::pop_ready_regions`
   (`activity_profile.cpp:424`, `pop_next` at `437`), whose two branches matter:
   `activity_profile.cpp:447-448` (buffer-size gate) and `activity_profile.cpp:450-452`
   (`find_first_activity_boundary` capped at `max_region_size=300`, then
   `find_best_cut_site(end, min_region_size=50)` at `494-505` when the cap is reached).
4. Consequence (observed, §1.2/§1.5): with halo end 10020500 the region
   `[10020199,10020488]` becomes `[10020199,10020283] + [10020284,10020499]`, and each of the two
   resulting AssemblyRegions is assembled separately (`calling_pipeline.cpp:15319-15364` recursive
   `run(...)` children with `force_single_calling_region=true`), so the EventMap group's **active span
   becomes `10020227-10020228`** (plus `10020428-10020437` when that is inside the window) instead of
   `10020227-10020437` (`[FASTGATK_RCM_EVENT_TRIM]`, `calling_pipeline.cpp:7842-7847`).
5. RCM construction: `build_profile_local_reference_blocks` (`calling_pipeline.cpp:13111`) walks the tile's
   `profile_regions` and, per profile, builds the segment boundary list from
   `{profile.start, region.active_start, region.active_end+1, profile.end+1}` (`13513-13522`). For tile#1
   profile `[10020199,10020283]` with the group active span `[10020227,10020228]` the segments are
   `[10020199,10020226]`, `[10020227,10020228]`, `[10020229,10020283]`.
6. Segment → owner resolution (`13527-13540`):

   ```cpp
   const auto owner = std::find_if(relevant_regions.begin(), relevant_regions.end(),
       [&](const auto* region) {                     // 13527-13532
           return output_start >= region->active_start && output_end <= region->active_end;
       });
   if (owner != relevant_regions.end())
       append_segment(output_start, output_end, (*owner)->start, (*owner)->end, *owner);
   else
       append_segment(output_start, output_end, padded_start, padded_end);   // 13540: NO event_owner
   ```

   Segment `[10020229,10020283]` fails the containment test for **every** relevant region (the only one in
   that profile is `[10020227,10020228]`), so it takes the flank branch at `13540` with
   `event_owner == nullptr`. The same happens in tile#1 profile `[10020284,10020499]` for the segment
   `[10020284,10020427]` (the profile's only relevant region is `[10020428,10020437]`).
   The union of those two flank segments is exactly 0-based `[10020229,10020427]` = 1-based
   **`10020230-10020428`** — the observed degenerate block, to the base.
7. `append_segment` (`13155`) only looks up the owning PairHMM `Result` when it was given an
   `event_owner` (`13173`: `if (event_owner != nullptr && event_region_results != nullptr)`). With
   `event_owner == nullptr`, `owner_result` stays null, therefore:
   * the read set stays the ordinary `RegionReadSelector::select(padded)` pileup (no
     `EVENT_READS`/`EVENT_REALIGNED` filter), and
   * **no `RealignedReadProjectionByRecord` is built** (`13334-13440`) and none is passed to
     `build_reference_blocks` (`13471-13477`), so `calculate_reference_confidence` takes the
     non-projected `read_is_indel_informative` branch (`calling_pipeline.cpp:3164-3173`).
8. Numerics: with no realignment the host indel-informativeness loop (`3150-3180`) yields
   `informative == 0` (`depth` is still non-zero — the pileup exists). At `3184-3195` the indel model with
   `nInformativeReads == 0` produces an all-zero PL row and hom-ref log-score 0, which is *less confident*
   than the SNP row, so `getGLwithWorstGQ` selects it → `GQ=0`, `PL=0,0,0`. Because every position in the
   gap has the same `gq_band`, `append_profile_reference_block` (`13052-13062`) merges them into a single
   block; `restrict_stream_gvcf_to_core` (`hc_call.cpp:6668`) and `RegionGvcfStitcher::merge`
   (`hc_call.cpp:5864`) then extend/join it across the core boundary.
9. GATK never takes that flank path here: the equivalent Java transaction is one AssemblyRegion whose
   active span is `10020227-10020437` (`-L 20:10019901-10020710`), which is exactly what the non-streamed
   native run reproduces byte-for-byte (48 rows).

**Responsible code:**

* **proximate site** (where the degenerate numbers are produced):
  `fastgatk-native/src/calling_pipeline.cpp:13527-13540` — the segment→owner containment lookup and its
  `else` flank branch at **line 13540** (guarded by `13173`, fed by the boundary list at `13513-13522`);
  the numeric degeneration itself is `calling_pipeline.cpp:3150-3180` → `3184-3195` (`informative == 0` ⇒
  all-zero indel PL wins `getGLwithWorstGQ` ⇒ `GQ=0 / PL=0,0,0`).
* **root cause** (why the owner spans are wrong for a tile):
  `fastgatk-native/src/hc_call.cpp:6658-6665` (tile options carry the *halo* interval), sourced from
  `hc_call.cpp:6459-6480` `input_halo_intervals`, acting through the region segmentation
  `fastgatk-kernels/src/activity_profile.cpp:424/437-452/494-505`.

---

## 3. Hypotheses REFUTED (with the evidence that refuted them)

1. **"D3 is caused by `build_profile_local_reference_blocks` being an approximate/different builder than
   the locus-scoped `build_reference_blocks`."** REFUTED. Both modes use the same
   `build_profile_local_reference_blocks`: in the non-streamed run `calling_regions` is non-empty and
   `assembly_region_independent_pass` is false, so the partition-merge site
   `calling_pipeline.cpp:15671-15681` is the RCM builder (confirmed by the single
   `BUILD_INPUT profiles=3` header after the children's header-less passes). Feeding that *same* builder
   the tile's window reproduces the tile's output exactly (§1.4), so the builder is not the discriminator.
2. **"The tile's read batch is wrong (too small / truncated halo reads), so its RCM is degenerate."**
   REFUTED. The tile's read batch and the `-L 20:10019901-10020500` run's read batch are the same
   interval-scoped batch (both are reads overlapping the halo window), yet the full-interval run over the
   same reads gives healthy values; and tile=810 (largest batch) is healthy while tile=500 (smaller batch)
   is degenerate. The controlling variable is the **interval window**, not the reads.
3. **"`GQ=0 / PL=0,0,0` over this span is a legitimate no-variation result, so the tile is not wrong."**
   REFUTED by the GATK baseline. GATK 4.6.2.0 and non-streamed native both emit the 29 healthy blocks over
   `10020230-10020428` (48 rows total, byte-identical); the degenerate block appears only when a tile core
   boundary (or a shortened group) falls in the span.
4. **"The degenerate values are created by `RegionGvcfStitcher::merge` / `restrict_stream_gvcf_to_core`
   acting on otherwise-healthy blocks."** REFUTED (extends track-D's per-locus measurement). A
   **non-streamed** run — no stitcher, no core restriction — emits the same degenerate block
   `10020230-10020428 0/0:26:0:21:0,0,0` (§1.4, §1.5). The stitcher only extends/joins blocks that are
   already degenerate.
5. **"`--assembly-region-padding` is an independent cause."** REFUTED as an independent cause. Padding
   enters only through `context` at `hc_call.cpp:6459-6463`, i.e. it moves the same halo window; the whole
   effect is reproducible with `-L` alone and no streaming at all (§1.5).
6. **"The tile's `tile_options.interval_start/interval_end` are consumed by the RCM builder as its
   locus/evidence domain."** REFUTED. `build_profile_local_reference_blocks` never reads
   `options.interval_start/interval_end`; per segment it *rewrites* them from the segment bounds
   (`calling_pipeline.cpp:13300-13307`,
   `local_options.intervals = {profile.tid, output_start, output_end+1}`). The tile bounds act **only
   indirectly**, through the ActivityProfile/AssemblyRegion partition — which is why the failure shows up
   as a *missing event owner*, not as an out-of-range locus domain.
7. **"The tile's per-locus activity values differ from the whole run's, so the region split is a numerical
   activity difference."** REFUTED for the loci that are observable: with
   `FASTGATK_DEBUG_REGION_SCHEDULING=1`, the `[FASTGATK_ACTIVITY_SEED]` signal values are identical in the
   two windows over the relevant range (positions 10020470-10020510 all `1.0`; the only differences are in
   the tail beyond the window end, e.g. 10020591+). The window-dependent split is therefore a
   **segmentation/pop-timing** effect of `BandPassProfile`, not a difference in activity values at the
   split position.
   *(Caveat: `ACTIVITY_SEED` prints only loci with `active != 0`, so the band-pass state at
   non-picket positions such as 10020489 is not directly observable without instrumentation — see §4.)*

---

## 4. What remains unproven, and the single next experiment

**Established:** the degenerate span; the owner-miss → flank → no-realignment → `informative=0` →
`GQ=0/PL=0,0,0` chain; the window-dependence of the whole effect (byte-identical reproduction).

**Not established:** *which* branch of `BandPassProfile::pop_next`
(`fastgatk-kernels/src/activity_profile.cpp:437-460`) fires for which window. The `-L` end sweep shows a
non-monotone flip — degenerate for ends `10020500-10020529`, healthy again from `10020530` (with a
different cut at `10020494`), healthy below `10020500` — which is the signature of the two branches
(`activity_profile.cpp:447-448` buffer-size gate vs `450` boundary scan vs `452` `find_best_cut_site`),
not of a change in state values. I could not observe those decisions without instrumentation, and I am not
permitted to rebuild.

**Single next experiment (requires one temporary instrumentation line + a rebuild of the HC binary):**
print, inside `BandPassProfile::pop_next` (`activity_profile.cpp:437-460`), one line per pop:
`raw_start_`, `states_.size()`, `is_active`, `force`, whether `find_best_cut_site` was used, and the
returned `(first_start, last_start)`. Run it on the three non-streamed windows
`-L 20:10019901-10020520` (expected: cut at `10020283`), `-L 20:10019901-10020530` (expected: boundary at
`10020494`) and `-L 20:10019901-10020710` (expected: boundary at `10020488`). This decides whether the
tile's split is (a) the 300-state cap + `find_best_cut_site` firing because the truncated window makes the
buffer long enough to reach the cap, or (b) a genuine inactive state that only exists in the shorter
window. That in turn tells whether the defect must be fixed in the *segmentation* (make region boundaries
traversal-window-independent) or in the *RCM owner matching* (accept a group as owner for positions inside
its padded span rather than its active span).

**Suggested rebuild-free partial check** (weaker, but doable now): run
`FASTGATK_DEBUG_REGION_SCHEDULING=1` with `-L 20:10019901-10020520` and `-L 20:10019901-10020530` and
diff the `[FASTGATK_REGION_BEGIN] … core=…` lines — if the *group* cores are identical while the *profile*
boundaries differ, the split is a pure segmentation artifact of the truncated traversal window.

---

## 5. Raw evidence kept

* `.diag/r2-tiles/t{100,200,300,405,500,600,700,810}.{vcf,err}` — tile sweep with `FASTGATK_DEBUG_RCM`.
* `.diag/r2-halo/L_full.*`, `L_halo.*`, `L_halo2.*` (= `-L 20:10019901-10020500`, byte-identical to
  tile=500), `end*.{vcf,err}` — the window-end sweep; `t200_exact_halo.*`, `L_halo2.*` — exact tile-halo reproductions.
* `.diag/r2-halo/{whole,t500_halo,t600_halo}.seed` — `FASTGATK_DEBUG_REGION_SCHEDULING` activity seeds.
* `.diag/r2-halo/*.igv` — `--assembly-region-out` dumps (aggregate path only; note `hc_call.cpp:6834`
  rejects it under `--stream-by-region`, which is why the halo-window equivalence was used instead).
* `.diag/r2-scratch/`, `.diag/r2-scratch2/` — first-pass traces (whole/t300/t500 RCM + topology).

No production file was modified; `git status` was not used to change anything.
