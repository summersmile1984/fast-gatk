# R3 — Independent verification of the claimed D1 phasing mechanism

Scope: **read-only**. No production source file was edited, nothing was rebuilt, no `*.md` at the
repo root was touched. The only file created is this report. All scratch lived inside a single
`mktemp -d` per command (per the sandbox note that `/tmp` does not survive between bash calls).

Subject under test — the claim recorded in
`fastgatk-native/evidence/2026-09-10-parallel-audit/track-d-streaming-rootcause.md` (D1 section):

> (C1) the `haplotype_map` **unions** `somatic_candidate_haplotype_indices` **from `result` and from
> every `assembly_region_likelihood_results` owner into one locally-numbered namespace**;
> (C2) `total_available_haplotypes` (used by an `always_apart` size test) **changes with the tile**;
> (C3) the phasing block is pinned to `gvcf()` at `fastgatk-native/src/hc_call.cpp:4085-4195`.

## Verdict (short)

| Sub-claim | Verdict |
| --- | --- |
| **C1** — cross-owner *union* of haplotype ordinals into one namespace | **REFUTED** (source structure + trace + replay simulation, three independent lines) |
| **C2** — `total_available_haplotypes` feeding `always_apart` is the quantity that flips, and it varies with the tiling | **CONFIRMED**, and isolated in a controlled A/B where *only* that quantity differs |
| **C3** — the phasing block lives in `gvcf()` near 4085–4195 | **CONFIRMED with a precision defect**: the block is `hc_call.cpp:4068-4246`; the cited range 4085–4195 ends at line 4195 (`continue`) and therefore **excludes the very `always_apart` predicate** it names, which is defined at 4197–4198 |

**Overall: PARTIALLY CONFIRMED.** The causal chain "tile → different haplotype sets/total →
`always_apart` size test flips → spurious `|` + PGT/PID/PS" is real, but the mechanism attributed to
it is **not** a cross-owner union — the map is rebuilt per phase context (no union ever happens), and
`physical_phases` (not `haplotype_map`) is the only cross-context accumulator. `result` is skipped as
a phase context in every whole-interval run (`haplotype_rows=0` vs `candidates=10`); it does pass the
guard in some tiles (e.g. tile=100: `candidates=2 haplotype_rows=2`) but that did not change any
phasing outcome measured here.

Two further findings that change how D1 should be described are in §5 and §6: D1 is **reproducible
with no streaming at all**, and **pinned GATK emits the identical phased row** for the
tile-equivalent window. D1 is therefore not "invented phasing" in the sense of violating GATK's rule;
it is the *interval/window dependence of the phase context*, the same class of defect as D2.

---

## 1. Reproduction and command set

```bash
B=fastgatk-native/build/fastgatk-hc-call
COMMON=(-R fixtures/chr20/ref20mnp.fasta -I fixtures/chr20/mnp.bam \
        -L 20:10019901-10020710 --emit-ref-confidence GVCF --max-mnp-distance 1 \
        --threads 1 --add-output-vcf-command-line false)
FASTGATK_DEBUG_PHASE=1 $B "${COMMON[@]}"                          -O whole.g.vcf   # 48 records
FASTGATK_DEBUG_PHASE=1 $B "${COMMON[@]}" --stream-by-region 500   -O t500.g.vcf
# documented D1 repro (pad=300):
FASTGATK_DEBUG_PHASE=1 $B "${COMMON[@]}" --stream-by-region 500 --assembly-region-padding 300 -O p300.g.vcf
```

D1 reproduced at `--stream-by-region 500` (default pad=100) — **and the row differs in more than
GT**, contrary to the impression given by the evidence doc's D1 example:

```text
non-streamed / GATK  20 10020228 . G A,<NON_REF>  ... 0/1:4,8,0:12:99:302,0,129,314,153,467:3,1,5,3
tile=500 (pad=100)   20 10020228 . G A,<NON_REF>  ... 0|1:2,4,0:6:57:0|1:10020228_G_A:148,0,57,154,69,223:10020228:2,0,3,1
tile=500 (pad=300)   20 10020228 . G A,<NON_REF>  RAW_MQandDP=43200,12 ... 0|1:4,8,0:12:99:0|1:10020228_G_A:302,0,129,314,153,467:10020228:3,1,5,3
```

The documented D1 example is the **pad=300** case (AD/DP/GQ/PL unchanged, only phasing added) — it is
reproduced exactly. At pad=100 the same tile also changes AD/DP/GQ/PL, which is a *different*
(collateral) divergence, not part of D1.

## 2. The exact code (read-only inspection)

`gvcf()` is declared at `hc_call.cpp:3688`; the phasing block is `hc_call.cpp:4068-4246`.
Byte-exact structure (whitespace as in the file — note the misleading 8-space indent at 4085):

```text
4068:    std::vector<PhysicalPhase> physical_phases(candidate_groups.size());   // OUTSIDE the context loop
 4076:    if (sample_ploidy == 2) {
 4077:        std::vector<const fastgatk::calling::Result*> phase_contexts;
 4078:        phase_contexts.push_back(&result);
 4079:        for (const auto& owner : result.assembly_region_likelihood_results)
 4079:            if (owner != nullptr) phase_contexts.push_back(owner.get());
 4080:        for (const auto* phase_context : phase_contexts) {          // <-- opens context loop
 4081:            if (phase_context == nullptr ||
 4082:                phase_context->somatic_candidate_haplotype_indices.size() !=
 4083:                    phase_context->candidates.size())
 4084:                continue;                                          // skips mismatched contexts
 4085:        std::vector<std::set<std::uint32_t>> haplotype_map(candidate_groups.size());   // <-- INSIDE the loop
 4086:        for (std::size_t group_index = 0; ...) { ... insert members ... }   // 4138-4139
 4150:        std::set<std::uint32_t> all_called_haplotypes;              // INSIDE the loop
 4153:        const auto total_available_haplotypes = all_called_haplotypes.size();
 4197:                const bool always_apart = call_haplotypes.size() + comp_haplotypes.size() ==
 4198:                    total_available_haplotypes && ...
 4244:        }    // closes group-decision loop
 4245:        }    // closes the phase_context loop  <-- haplotype_map / total / phase_group die here
 4246:    }        // closes if (sample_ploidy == 2)
```

Consequences, all of which contradict C1:

1. `haplotype_map` (4085) is **re-created on every iteration of the context loop**, and the whole
   decision (`all_called_haplotypes` 4150, `total_available_haplotypes` 4153, `phase_group`,
   `phase_01`, `unique_counter`, `unphasable`) is inside that iteration. **There is no cross-owner
   union of haplotype ordinals** — each context runs the algorithm on *its own* candidates and
   *its own* locally-numbered haplotypes.
2. `result` participates only if `result.somatic_candidate_haplotype_indices.size() ==
   result.candidates.size()`. In every whole-interval run measured here the diagnostic prints
   `haplotype_rows=0` against `candidates=10`, so **`result` is skipped by the 4081-4084 guard**.
   (It does pass in some tiles, e.g. tile=100 prints `candidates=2 ... haplotype_rows=2`.)
3. The only vector that accumulates across contexts is `physical_phases` (4068), which is *outside*
   the loop. This looks like the object the claim conflated with `haplotype_map`.

## 3. Empirical refutation of C1 — decode 1: the trace shows the map is reset

`FASTGATK_DEBUG_PHASE=1`, tile=500, **single `gvcf()` call** (one `[FASTGATK_PHASE_INPUT]` header,
`candidates=4`), raw stderr excerpts in order:

```text
[FASTGATK_PHASE_EVENT] pos=10019967 ref=C alt=G haps=1,
[FASTGATK_PHASE_EVENT] pos=10019969 ref=T alt=G haps=1,
[FASTGATK_PHASE_EVENT] pos=10020228 ref=G alt=A haps=
[FASTGATK_PHASE_EVENT] pos=10020229 ref=T alt=G haps=
[FASTGATK_PHASE_OUTPUT] pos=10019967 phase=0|1 pid=10019967_C_G
[FASTGATK_PHASE_OUTPUT] pos=10019969 phase=0|1 pid=10019967_C_G
[FASTGATK_PHASE_EVENT] pos=10020228 ref=G alt=A haps=0,       <-- pass 2: NON-EMPTY
[FASTGATK_PHASE_EVENT] pos=10020229 ref=T alt=G haps=1,
[FASTGATK_PHASE_OUTPUT] pos=10020228 phase=0|1 pid=10020228_G_A
[FASTGATK_PHASE_OUTPUT] pos=10020229 phase=1|0 pid=10020228_G_A
[FASTGATK_PHASE_EVENT] pos=10020228 ref=G alt=A haps=       <-- pass 3: EMPTY AGAIN
[FASTGATK_PHASE_EVENT] pos=10020229 ref=T alt=G haps=
```

`std::set::insert` can only grow a set. A group printing `haps=0,` and then, **later in the same
`gvcf()` call**, `haps=` (empty) is only possible if the container was **destroyed and re-created
between the two prints** — i.e. `haplotype_map` is per-context. The same pattern appears at
tile=500/pad=300 (`10020228 haps=0,1,` then `haps=`) and at tile=700/810
(`10020429..10020438 haps=1,3,5,` then `haps=`).

⇒ **Under the claimed union mechanism the third print would have to read `haps=0,`.** It does not.

## 4. Empirical refutation of C1 — decode 2: replay simulation of both models

I replayed the §2 algorithm exactly (transcribed from source) against the observed per-pass haplotype
sets, in two variants:

* **PER-CONTEXT** — decision run independently inside each phase context, using that context's own
  total (`total_available_haplotypes` = union over the groups that context owns). Predictions from
  different contexts accumulate in `physical_phases`.
* **UNION-OWNERS** — the claimed model: one `haplotype_map` accumulating every context's inserts,
  one decision pass, one total.

Pass boundaries were reconstructed from the trace by the position-wrap rule (the group loop walks
candidate groups in ascending position order, so a non-increasing position starts a new context
pass). Compared against the observed `[FASTGATK_PHASE_OUTPUT]` sequence:

| config | observed PHASE_OUTPUT | PER-CONTEXT | UNION-OWNERS |
| --- | --- | --- | --- |
| whole (non-streamed) | 10019967, 10019969, 10020429/31/34/38, 10020679, 10020680 | **MATCH** | **MISMATCH** (drops 10020679/10020680) |
| tile=100 | as above **+** 10020228 `0\|1`, 10020229 `1\|0` | **MATCH** | MATCH (coincidence) |
| tile=200 | as tile=100 | **MATCH** | MATCH (coincidence) |
| tile=300 | as whole | **MATCH** | MATCH (coincidence) |
| tile=405 | as tile=100 | **MATCH** | **MISMATCH** |
| tile=500 | as tile=100 | **MATCH** | **MISMATCH** |
| tile=600 | as whole | **MATCH** | MATCH (coincidence) |
| tile=700 | as whole | **MATCH** | **MISMATCH** |
| tile=810 | as whole | **MATCH** | **MISMATCH** |
| tile=500 + pad=300 | as tile=100 | **MATCH** | **MISMATCH** |

Per-context: 10/10 exact matches, including phase polarity and relative order.
Union-owners: 6/10 mismatches.

For tile=500/pad=300 the replay prints, verbatim:

```text
  pass total=4 sets=[(10020228, [0, 1]), (10020229, [2, 3])] -> [(10020228, '0|1'), (10020229, '1|0')]
  pass total=3 sets=[(10020429, [1, 3, 5]), (10020431, [1, 3, 5]), (10020434, [1, 3, 5]), (10020438, [1, 3, 5])] -> [...]
  pad300 own : [(10019967,'0|1'),(10019969,'0|1'),(10020228,'0|1'),(10020229,'1|0'),...]  PER-CONTEXT: MATCH   UNION: MISMATCH
```

The **whole-interval (non-streamed) mismatch is self-contained and decisive**: the union of all
passes gives `{0,1} ∪ {2,3} ∪ {1,3,5} ∪ {1} ∪ {2} ∪ {0} = {0,1,2,3,5}`, total **5**; with total 5 the
pair `10020679={2}` / `10020680={0}` satisfies neither `always_together` nor `always_apart`, so the
union model predicts **no phasing there**. But native *and* pinned GATK both emit
`10020679 0|1` / `10020680 1|0` (PGT `10020679_AC_TA`) for that very run. Under the union model that
output is unreachable. Note this argument does not depend on my wrap heuristic: with a single
decision pass over the whole run, the pair can never be "apart" against total 5.

## 5. C2 CONFIRMED, with a controlled isolation where **only** `total_available_haplotypes` changes

Two configurations produce **identical haplotype sets for the D1 pair** but different totals:

```text
non-streamed,  -L 20:10019901-10020710        pass total=5   sets={10020228:[0,1], 10020229:[2,3], 10020429..10020438:[1,3,5]x4}
tile=500 pad=300                              pass total=4   sets={10020228:[0,1], 10020229:[2,3]}
```

Replay of `always_apart` on the D1 pair (`A={0,1}`, `B={2,3}`, disjoint):

* total = 5 → `|A|+|B| = 4 != 5` → **`always_apart` false** → no `PHASE_OUTPUT` → `GT 0/1`.
* total = 4 → `|A|+|B| = 4 == 4` and disjoint → **`always_apart` true** → `PHASE_OUTPUT`
  `10020228 0|1` / `10020229 1|0`, pid `10020228_G_A` → emitted as `0|1` + `PGT/PID/PS`.

Observed: the whole-interval run emits `0/1` (no PGT/PID/PS) and the pad=300 run emits
`0|1 … 0|1:10020228_G_A … :10020228: …`, exactly as the replay predicts. Because the two sets are
bit-identical, the only variable that can flip the predicate is `total_available_haplotypes`
(line 4153); the "locally-numbered namespace" of C1 plays no role here.

**Both configurations are matched by the PER-CONTEXT model 3/3 passes**; the union model is not.

### Tile-vs-D1 correlation explained by the same quantity

`D1 present` ⟺ the context owning `10020228`/`10020229` does **not** also own `10020429..10020438`:

| tile | sets in the D1 pair's context | total | D1 |
| ---: | --- | ---: | --- |
| whole | `{0,1}`, `{2,3}` **+ `{1,3,5}`×4** | 5 | no |
| 100 | `{0}`, `{1}` | 2 | yes |
| 200 | `{0}`, `{1}` | 2 | yes |
| 300 | `{0,1}`, `{2,3}` **+ `{1,3,5}`×4** | 5 | no |
| 405 | `{0}`, `{1}` | 2 | yes |
| 500 | `{0}`, `{1}` | 2 | yes |
| 500 + pad=300 | `{0,1}`, `{2,3}` | 4 | yes |
| 600 | `{0,1}`, `{2,3}` **+ `{1,3,5}`×4** | 5 | no |
| 700 | `{0,1}`, `{2,3}` **+ `{1,3,5}`×4** | 5 | no |
| 810 | `{0,1}`, `{2,3}` **+ `{1,3,5}`×4** | 5 | no |

The non-monotonicity over tile size (405/500 yes, 600/700 no) is fully accounted for: it is not a
size effect but the question of whether the tile's AssemblyRegion boundaries put the
`10020429..10020438` candidate group in the *same* owner/context as `10020228`/`10020229`.

## 6. Two findings that reclassify D1 (not in the audited claim)

**(a) D1 is reproducible with no streaming whatsoever.** Plain non-streamed runs:

| `-L` (no `--stream-by-region`) | `10020228` FORMAT | records |
| --- | --- | --- |
| `20:10019901-10020710` | `0/1:4,8,0:12:99:302,0,129,314,153,467:3,1,5,3` | 48 |
| `20:10019901-10020428` | `0\|1:2,4,0:6:57:0\|1:10020228_G_A:148,0,57,154,69,223:10020228:2,0,3,1` | 8 |
| `20:10019901-10020400` | `0\|1:…` (identical to the tile=500/pad=100 row) | 8 |
| `20:10019901-10020350` | `0\|1:…` (identical) | 8 |
| `-L` start sweep 10019801…10020501, end fixed 10020710 | always `0/1` | 4–48 |

So D1 is an **interval-extent dependence** exactly like D2, not a streaming-plumbing defect; the
streamed tile merely selects a window that triggers it. (Same pattern as the D2 correction already
recorded at the top of `gvcf-stream-by-region-divergence-20260910.md`.)

**(b) Pinned GATK produces the same phased row for that window, byte-for-byte.**

```bash
third_party/jdk17/bin/java -Xmx1g -jar third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
  HaplotypeCaller -R fixtures/chr20/ref20mnp.fasta -I fixtures/chr20/mnp.bam -L 20:10019901-10020400 \
  -ERC GVCF --max-mnp-distance 1 -O gatk_trunc.g.vcf --native-pair-hmm-threads 1 \
  --create-output-variant-index false --add-output-vcf-command-line false
```

* native **non-streamed** `-L 20:10019901-10020400` vs GATK same window: `diff` of the data lines is
  **empty (8/8 records byte-identical)**.
* GATK's row at `10020228` is `0|1:2,4,0:6:57:0|1:10020228_G_A:148,0,57,154,69,223:10020228:2,0,3,1`,
  i.e. **identical to native's streamed tile=100/405/500 row**.

Therefore the phasing decision itself is *GATK-faithful for the haplotype set the tile actually sees*;
what deviates is the **scope** of that set relative to the un-tiled interval. Calling D1 "凭空产生
phasing / spurious phasing invented by streaming" is a correct statement only relative to comparing
two different windows; per GATK semantics the `0|1` is a legitimate output of the smaller context,
and GATK's own full-interval answer (`0/1`) is equally legitimate for the full context. The real
parity defect is that `--stream-by-region` **re-scopes the phase context**, so streamed ≠ non-streamed
for the same `-L` — which is precisely the same "chunking invariance" defect as D2/D3.

## 7. Responsible location (my conclusion)

* **Function:** `std::string gvcf(const fastgatk::io::HtsReader& reader, …)` — `fastgatk-native/src/hc_call.cpp:3688`; phasing block `4068-4246`.
* **The quantity that flips:** `total_available_haplotypes` — `hc_call.cpp:4150-4153`
  (`all_called_haplotypes` → `.size()`), consumed by
  * `always_apart` — `hc_call.cpp:4197-4202` (the size test `|A| + |B| == total`, plus disjointness),
    and
  * `always_together` via `call_is_on_all_alt_haplotypes` / `comp_is_on_all_alt_haplotypes`
    — `hc_call.cpp:4161-4170`.
* **What makes it vary:** which candidate groups each phase context owns — i.e. `phase_context->candidates`
  membership and `phase_context->somatic_candidate_haplotype_indices[candidate_index]`
  (`hc_call.cpp:4130-4139`), where `phase_context` ranges over `&result` plus each
  `result.assembly_region_likelihood_results[*]` owner (`hc_call.cpp:4077-4080`). That membership is
  a function of how `fastgatk::calling::run` partitioned the interval into AssemblyRegions:
  the per-tile call at `hc_call.cpp:6666` inside `run_region_streaming` (non-streamed counterpart at
  `hc_call.cpp:6976`).
* **Not responsible:** a cross-owner haplotype-ordinal union. No such union exists.

## 8. What is proven vs. unproven

**Proven (this session, observed or read from the shipped source)**

* C2: `total_available_haplotypes` entering `always_apart` flips the D1 decision, isolated in an A/B
  (`whole`, total 5 → `0/1`) vs (`tile=500 pad=300`, total 4, identical sets → `0|1` + PGT/PID/PS).
* C1 refuted by (i) brace structure (`haplotype_map` at 4085 is inside the loop opened at 4080 and
  closed at 4245), (ii) the trace showing `haps=0,` → `haps=` for the same group inside one `gvcf()`
  call, (iii) a replay in which the union model cannot produce the non-streamed baseline's own
  `10020679`/`10020680` phasing. `result` is skipped (`haplotype_rows=0` vs `candidates=10`) in all
  whole-interval runs examined.
* The PER-CONTEXT replay reproduces the observed `PHASE_OUTPUT` sequence exactly in all 10
  configurations tested (whole, tiles 100/200/300/405/500/600/700/810, and tile=500+pad=300).
* C3 approximately right but self-inconsistent: cited 4085–4195 excludes `always_apart` (4197).
* D1 is reproducible non-streamed (`-L 20:10019901-10020400`), and native non-streamed == GATK
  byte-identical on that window; GATK's row there == the streamed D1 row.
* The documented D1 repro (tile=500 **pad=300**) reproduces exactly; at pad=100 the same tile
  additionally changes AD/DP/GQ/PL, which the evidence doc's D1 example does not show.

**Not proven / residual uncertainty**

* I could not print `total_available_haplotypes` directly (a rebuild was forbidden). The totals in
  §4/§5 are **inferred** by replay from the `[FASTGATK_PHASE_EVENT]` sets. The inference is
  over-determined in the cases that matter (the whole-run `10020679`/`10020680` argument holds under
  any pass splitting), but the individual totals for tiles 100/200/300/600/700/810 rest on the
  position-wrap pass reconstruction and on the assumption that `candidate_groups` is ordered by
  ascending position (observed consistent in every pass printed, not asserted by a printed index).
* **Why** a given tile's AssemblyRegion partitioning puts `10020228`/`10020229` in a different owner
  from `10020429..10020438` is not explained here — that is the same open question as D2/D3 (the
  interval-span/AssemblyRegion-partitioning dependence of `calling::run`), and it needs the
  activity-profile/region-construction trace, not the phasing gate.
* Not tested: other fixtures/contigs, `--stream-by-contig`, `--max-alternate-alleles`,
  non-diploid ploidy, and multi-sample layouts. The `result`-as-phase-context branch
  (`haplotype_rows == candidates`, seen at tile=100) was observed to fire but I did not construct a
  case where it changes a phasing outcome.
* Whether `physical_phases` cross-context accumulation (a context may phase a group and later
  contexts may overwrite it) can itself produce an inconsistency was not exercised — no case with two
  phasing passes over the same position appeared in this fixture.

No patch is proposed (per task constraint).
