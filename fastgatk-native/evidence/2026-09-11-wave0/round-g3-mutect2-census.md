# Round G3 — Mutect2 owner-selection census (read-only)

**Question.** The HaplotypeCaller side of the "published record annotated from an AssemblyRegion owner
that never requested PairHMM for that allele" family (D2) has been censused (15 sites; 13 safe by
CONSEQUENCE-1, one vulnerable, two fixed). The Mutect2 side was never audited. Enumerate **every**
site in `fastgatk-native/src/mutect2_tool.cpp` (plus the helpers it calls) that selects among the
owners in `Result::assembly_region_likelihood_results`; classify each **(A) VULNERABLE** or
**(B) SAFE** under the established criterion; verify the "no `calculate_output_variant_annotations`
in Mutect2" claim; confirm at least one classification empirically if cheap.

**Constraints honoured.** No file under `fastgatk-native/src`, `fastgatk-native/include`,
`fastgatk-kernels/*` or any `CMakeLists.txt` or repo-root `*.md` was edited; nothing was rebuilt; no
`git` state was changed; the ctest suite was not run (targeted runs only). Every experiment ran inside
one `bash` invocation with an explicit timeout and all scratch in `tempfile.TemporaryDirectory`
(`/tmp` does not persist between invocations). The only file created by this round is this report.

---

## 0. TL;DR

* **14 owner-selection sites** in `mutect2_tool.cpp` (§2). Mutect2 uses a *different* within-`Result`
  sibling scan from HaplotypeCaller: it resolves identity through `candidate_index()`
  (`mutect2_tool.cpp:2368-2400`), which is **ownership-priority aware** (it prefers the
  `likelihood_candidate_event_map_owned` row over a raw sibling), and every owner loop prefers the
  owner whose own *active core* contains the locus (priority 2) over a padded halo (priority 1).
* **5 sites are (A)** — `coverage_depth_for_candidate` (3460), `apply_grouped_somatic_depth` (3244),
  `emitted_somatic_fragment_depth` (3323), `somatic_annotation_counts` (4116) and the two reducer
  copies in `calculate_somatic` (4463) / `calculate_posterior` (5316). All of them accept an owner (or
  an identity-keyed sibling entry) whose likelihood rows were never written, and **none of them
  rejects the resulting all-zero/0/-inf table** — the sentinel is silently published as `0`.
* **9 sites are (B)** for a specific reason each; the two most important are
  `emitted_somatic_read_depth` (4388, safe by CONSEQUENCE-1 — its local returns `nullopt` without a
  finite row) and `somatic_annotation_counts`'s no-emitted-candidates branch (4213, **dead code**:
  every call site passes a non-null emitted list). `orientation_counts`'s owner loop (4912) is
  (A)-shaped but unreachable because the flattening reducer returns first.
* **The claim "mutect2_tool.cpp has no `calculate_output_variant_annotations` call" HOLDS** (verified
  by whole-tree grep and by enumerating every `fastgatk::calling::` symbol the tool uses). Mutect2
  uses its own annotation family instead; that family **is** owner-selecting and **can** consume a
  sentinel owner in the narrow case described in §3/R3.
* **Empirical result: a clean negative.** Native Mutect2 is byte-identical to pinned GATK 4.6.2.0 on
  **44/44** data rows of the whole HCC1143 chr20 fixture (a flattened `Result` with **590 owners**),
  identical across 5 different `-L` starts, identical under `--force-active true`, identical on the
  two record clusters, and **every** published record's locus produced populated annotation evidence
  (no missing `-1.0e299` REF row, no all-zero `SB`/`F1R2` published). I could not provoke the defect
  family on any Mutect2 fixture. Details and the exact falsifiers are in §4.

---

## 1. The criterion, and the Mutect2-specific structural facts it has to be applied to

**The criterion (as given).** A site is **(A)** iff it can select an owner using a probe that does not
consume likelihood rows (or consumes them only in a way a sentinel owner still passes); otherwise
**(B)**, with the invariant that makes it safe.

**INV-1 / INV-2 (unchanged, re-verified by reading):** `likelihood_candidate_read_context_ordinals` is
initialised to `UINT32_MAX` and written only from a PairHMM request
(`calling_pipeline.cpp:10547-10568`, inside `run_pairhmm` at `7167`); the rows are initialised to
`-inf` (`7211-7216`) and written only from a request (`11943-11950`). So *sentinel ⟺ this `Result`
holds no PairHMM request for that candidate*, and then every row of it is non-finite.

Four Mutect2-specific facts decide the classifications below. All four are **reasoned from code**
(not measured) unless marked.

**R1 — in somatic mode a candidate is requested by exactly one region.** The PairHMM group builder
assigns a candidate to its group only when `active_owner[candidate] == region_index`, i.e. the
region whose *active core* owns it (`calling_pipeline.cpp:7445-7465`, `7489-7497`), with the explicit
comment that "a candidate belongs to exactly one genotyping region". Every *other* AssemblyRegion
whose candidate list also contains that allele (halo-expanded neighbours) therefore has
`context_ordinal == UINT32_MAX` for it.

**R2 — `likelihood_candidate_event_map_owned` does NOT imply "requested".** It is computed per region
purely from that region's own EventMap
(`calling_pipeline.cpp:7406-7411` → `somatic_event_map_owns_candidate`), independently of the group
assignment in R1. A halo region's EventMap *does* contain the halo allele, so a **sentinel entry can
also be event-map-owned**. Consequence: `candidate_index()`'s `+16` ownership preference
(`mutect2_tool.cpp:2386-2388`) is a strong heuristic but is **not a sentinel filter**.

**R3 — the priority rule is the actual sentinel filter.** Every Mutect2 owner loop prefers the owner
whose own `calling_regions` *contain* the locus (`priority = 2`) over a halo (`priority = 1`). Each
owner `Result` is a child `run()` invocation with `force_single_calling_region`
(`calling_pipeline.cpp:15319-15352`, `force_single_calling_region = true` at `15320`), so
`part.calling_regions == {that region}` with
`active_start/active_end` the *unpadded core* (`15187-15188`). Active cores are disjoint, so for a
SNV at most one owner can be priority 2, and by R1 that owner is the same region that was assigned
the request. This is the mitigation the HaplotypeCaller `annotations_from_owner` path lacks — and it
is why most Mutect2 sites are safe *in practice*.

**R4 — the residual hole.** The priority-2 owner can itself be sentinel: an EventMap-owned candidate
inside its own core that received **zero** qualifying PairHMM requests (no read overlapped, all rows
filtered). Then the selected owner publishes an all-zero / `0` / `-inf` table. A second, narrower
variant of the same hole applies to the *merge-style* reducers (`calculate_somatic`,
`calculate_posterior`, `apply_grouped_somatic_depth`), which copy per identity key **without any
finiteness gate**: if one part holds two entries with the same `(tid,pos,ref,alt)` — the documented
"raw pileup description + EventMap-materialized description of the same allele"
(`mutect2_tool.cpp:2379-2384`) — the entry that is written last wins, and the raw one is the
sentinel one. The documented list order (raw first) makes the modeled row win; that ordering, not a
check, is what currently protects those three sites.

---

## 2. Site census — every owner-selection site in `mutect2_tool.cpp`

Line numbers are the current tree. "probe" = what the site uses to *accept* an owner/entry.
`candidate_index` = `mutect2_tool.cpp:2368-2400` (identity match, tie-broken by event-map ownership,
`graph_derived`, `forced_by_tumor`, recovery flags — it reads **no** likelihood row).
`priority` = the active-core preference of R3 (a calling-region containment test — also no row).

| # | line(s) | what it decides | probe used to accept | class | reason |
| --- | --- | --- | --- | --- | --- |
| 1 | `2798-2846` (loop `2807`) `belongs_to_active_event_map` inside `materialize_reference_confidence_eventmap_calls` | whether a reference-confidence EventMap call is *created* | `any_of` over **all** owners: identity match **and** `likelihood_candidate_event_map_owned[local] != 0` **and** containment in that part's own active core | **B** | No owner value is consumed at all (the decision is a pure ownership/containment predicate; `early true` on the first match, no value copied). The numeric payload is `somatic->tlod[index]` from the *flattened* result and is finite-gated (`2795-2798`), so a sentinel owner cannot inject `-inf`/`0` here. A sentinel halo owner can never satisfy the containment clause (R1/R3). |
| 2 | `3229-3283` (loop `3244`) `apply_grouped_somatic_depth` | FORMAT/FAD source scalars `somatic_fragment_depth/reference_count/alternate_count` for every merged candidate | `candidate_index(result, part.candidates[local])` + priority; copy is `if (priority < ownership) continue;` (equal priority ⇒ **last writer wins**) | **A** | The owned value is the part's fragment-depth table, computed from that part's own grouped rows. A sentinel entry contributes concrete **zeros** (`grouped_somatic_fragment_depth` with no rows ⇒ `active_groups` empty ⇒ all-zero output, `3202-3226`) and the copy is unconditional — no finiteness/emptiness check exists. Reachable only under R3/R4 (priority-2 owner sentinel) or a last-writer raw duplicate (R4). Impact is limited: this table is only a fallback for `normal_fad` (`6672-6688`) and `find_depth_row`. |
| 3 | `3314-3346` (loop `3323`) `emitted_somatic_fragment_depth` | **published** FORMAT/FAD for the emitted allele list (tumor and every normal/second tumour sample) | `has_all` = `candidate_index(part, emitted) != max` for every emitted allele, then max-`priority` owner | **A** | `candidate_index` consumes no likelihood row, and the local (`3291-3312` → `grouped_somatic_fragment_depth`) returns a **well-formed all-zero table** for a sentinel owner: it never requires a finite row and never returns `nullopt` once the identity is found (`3184-3190`). The callers at `6207-6208` / `6634-6636` therefore publish `FAD=0,0,…` through the normal `has_value()` path. Falls under R3/R4; note `has_all` can in principle fall back to a halo owner if the core part lacks an emitted allele. |
| 4 | `3446-3485` (loop `3460`) `coverage_depth_for_candidate` | INFO/DP (tumour core contribution and the normal-view contribution at `6249-6250`) | `candidate_index(*part, event)` + priority; `selected_priority` starts `-1`, so among **equal** priorities the **first** owner wins | **A** | Two independent gaps: (i) the probe is non-row-consuming; (ii) the accepted value is `part.candidates[local].variant_depth`, a scalar that is only refreshed from `pairhmm.annotation_depth` *when that is non-zero* (`calling_pipeline.cpp:17037-17041`), so for a sentinel entry it keeps its construction/import value (typically `0`, or a merged `depth` from `16244`) instead of this owner's evidence count. First-wins-among-halo makes the result depend on the owner vector order, i.e. on where `-L` starts. |
| 5 | `4100-4181` (loop `4116`) `somatic_annotation_counts` (emitted branch) | **published** FORMAT/F1R2, F2R1, SB and INFO MBQ/MFRL/MMQ/MPOS, for tumour and every normal/second tumour sample | `has_all` = `candidate_index(part, emitted) != max` for every emitted allele, then max-`priority` owner | **A (highest impact)** | This is the direct Mutect2 analogue of HC's D2 and the clearest (A). `candidate_index` consumes no row. The local (`3581-4098`) *does* consult the context ordinal — but only to **exclude** sentinel candidates from a locus group (`3771-3790`); when **every** entry of the group is sentinel it `continue`s and the function still `return output;` (`4097`) with an all-zero `SomaticAnnotationCounts`. There is no emptiness or finiteness gate, so the wrapper's `has_value()` at `6205`/`6213`/`6222` accepts zeros as authoritative and the writer publishes `F1R2=0,0`, `F2R1=0,0`, `SB=0,0,0,0`. A sentinel owner that is *event-map-owned* passes both the probe and the ownership-preferred `candidate_index` (R2). |
| 6 | `4183-4253` (loop `4213`) `somatic_annotation_counts` (no-emitted-candidates branch) | per-merged-candidate reduction over every owner | `somatic_annotation_counts(part, reads)` per part (recursion), then `find_merged_candidate` identity lookup + priority (`if (priority < ownership[...]) continue;`) | **B — dead code in the tool flow** | Verified by grep: the only call sites are `6205`, `6213`, `6222`, and **all three pass a non-null `emitted_candidates`**; the recursion at `4216` lives inside this same branch. The branch is therefore not reachable from `Mutect2`, and the `Mutect2` tool is its only caller. (If a future caller passed `nullptr`, this would be (A) for the same reason as #5 — the local's zero-table behaviour is identical.) |
| 7 | `4375-4414` (loop `4388`) `emitted_somatic_read_depth` | **published** FORMAT/AD, DP and INFO/DP | `has_all` = `candidate_index(part, emitted) != max`, then max-`priority` owner | **B (CONSEQUENCE-1)** | Same probe shape as #3/#5, but the local (`4267-4373`) refuses to produce a value without evidence: it builds PairHMM requests only from **finite** rows (`4326-4341`) and returns `std::nullopt` when `row_count == 0 || requests.empty()` (`4344`). A sentinel owner has no finite row, so it returns `nullopt`, the caller's `has_value()` test at `6226-6233`/`6631-6636` fails, and the previous (candidate-scalar/`find_depth_row`) value is kept. This is exactly the HC `derive_multiallelic_depths`/`joint_candidate_pl` safety argument. |
| 8 | `4416-4523` (loop `4463`) `calculate_somatic` flatten reducer | INFO TLOD/NALOD/NLOD/AF and the **emission gate** for every merged candidate | none — *every* owner is consumed; acceptance is by identity key (`merged_candidate_index`) + priority, `if (priority < ownership[...]) continue;` (equal priority ⇒ last wins) | **A** | The copy `output.tlod[merged] = local->tlod[local_index]` (`4494`) and its siblings (`4495-4514`) are **unconditioned on row finiteness** — the only guard in the loop is priority (`4474-4493`). A sentinel entry contributes the sparse kernel's initialised `0.0` (`fastgatk-kernels/src/somatic.cpp:795`, dense path `658`) and can overwrite another entry's real TLOD for the same identity (R4). Because `tumor_somatic->tlod` drives emission (`effective_tumor_lod_emit`), the observable symptom is a **dropped or invented record**, not just a wrong field. |
| 9 | `4903-5030` (loop `4912`) `orientation_counts(Result, ReadBatch, size_t)` | F1R2/F2R1 orientation counts for a candidate | identity scan over every part's candidates (no row consumption, **not even an index-bounds check**); max-`priority` owner | **B — unreachable in the production tree** | The chosen owner's counts *do* require a finite row per read, but with an **OR** gate (`4997-5001`: skip only when REF **and** ALT are both missing), so a sentinel owner would silently return `{0,0}` instead of failing. That would make this (A)-shaped. It is (B) only because of reachability: the owner loop is entered solely from `5395`, i.e. from the branch `calculate_posterior` reaches when the result has **no** sub-owners (`5277` enters
the flattening reducer and returns first for a flattened result) or from its own recursion at `4936` on a part, and parts are leaves (no nesting). The two other call sites (`7602`, `9131`) pass a `ReadBatch` and take the flat overload at `4876`, which has no owner selection. |
| 10 | `5298-5390` (loop `5316`) `calculate_posterior` flatten reducer | somatic / germline / artifact probabilities, contamination-adjusted AF, orientation-bias probability, somatic-germline-artifact log10 evidence | none — every owner consumed; identity key + priority, `if (priority < ownership[...]) continue;` (equal ⇒ last wins); per-part pairing with `normal->assembly_region_likelihood_results[part_index]` (index-paired, `5319-5323`) | **A** | Identical shape and identical gap to #8: the copies at `5351-5370` are unconditional, so a sentinel entry's initialised defaults (`somatic_probability = 0.0` at `5290`, `artifact_probability = 1.0` at `5292`) are published for a locus, changing the germline/artifact posterior that `FilterMutectCalls` consumes. Also note the index-pairing with the normal's owner vector, which relies on both samples partitioning identically. |
| 11 | `5617-5655` (loop `5646`) `first_eventmap_haplotype` | ALT ordering inside one merged VCF row (a `Number=A/R`-affecting serialization decision) | identity match (`matching_index`), first owner (in `result` then each part) that yields a non-missing rank | **B** | The consumed value is `somatic_candidate_haplotype_indices` — **EventMap/haplotype membership**, populated for every candidate from the graph state independently of requests (`calling_pipeline.cpp:10967-11040`). It is neither a likelihood row nor the context ordinal, so the sentinel mechanism does not apply. (Whether the *halo* owner's region-local haplotype numbering is the GATK-authoritative one is a separate, non-sentinel question.) |
| 12 | `5896-5930` (loop `5900`) `phase_owner_for_candidate` | PGT/PID/PS physical phasing | `candidate_index` (flattened first, priority 0) then max-`priority` owner; membership check `candidate >= owner->somatic_candidate_haplotype_indices.size()` | **B** | Same reason as #11: the owner-derived payload is EventMap haplotype membership, not a likelihood row and not the context ordinal. A sentinel entry with empty memberships merely drops the site from phasing (`phase_sites_by_owner` only receives sites with non-empty membership, `5934-5937`), which cannot be produced by the sentinel mechanism. |
| 13 | `6061-6095` (loop `6065`) `owner_for_candidate` | INFO/ECNT and ECNTH (potential-event counts) | `candidate_index` + max-`priority` owner | **B** | The consumed payload is the owner's `graph_haplotype_event_maps` / `calling_regions` (a graph property) plus the **flattened** `tumor_somatic->tlod` / `normal->likelihoods` scalars (`6100-6120`); no owner likelihood row is read, so a sentinel owner cannot degrade it. (Halo-vs-core EventMap counting is a separate family — see M3 in §4, where GATK agrees with native on the window that changes ECNT/ECNTH.) |
| 14 | `7250-7277` `variant_non_reference_likelihoods` | per-locus `<NON_REF>` TLOD/AF for `-ERC GVCF` | none — flat scan of `tumor.candidates` with `priority` from `likelihood_candidate_event_map_owned` / `graph_derived` | **B / not an owner-selection site** | Consumes only flattened vectors (`tumor_somatic->non_reference_tlod`), never a part; retained here for completeness because it shares the priority idiom. |

Remaining `assembly_region_likelihood_results` hits in the file are **not** owner-selection loops and
were checked line by line: `7920`, `8015`, `8736`, `8868` are option-flag assignments
(`partition_somatic_assembly_regions` / `retain_assembly_region_likelihood_results = true`), which is
itself worth recording — **both** the tiled and the aggregate Mutect2 path enable the flattened
owner-retention mode, so this whole census applies to *every* Mutect2 invocation, not just
`--stream-by-region`.

---

## 3. The `calculate_output_variant_annotations` claim, and what Mutect2 uses instead

**Claim: "mutect2_tool.cpp has NO `calculate_output_variant_annotations` call."  → HOLDS (verified).**
Whole-tree grep finds the declaration in `fastgatk-native/include/fastgatk/calling/pipeline.hpp:1314`,
the definition in `calling_pipeline.cpp:13762`, and exactly three call sites, **all in
`hc_call.cpp`** (`3523`, `3530`, `4573`, `4579`, `4903`, `4909`). Enumerating every
`fastgatk::calling::` symbol used by `mutect2_tool.cpp` confirms it: `Result`, `AssemblyCandidate`,
`GenotypeCall`, `run`, `ForcedAllele`, `Options`, `Likelihoods`, `AssemblyRegionAssembly`,
`SomaticGroupLikelihood`, `SomaticActivityFeatureMask`, `MultiallelicDepth`,
`CandidateReadRealignment` — no annotation entry point, and no `calculate_variant_annotations`
either (that function is reachable only through `calculate_output_variant_annotations`).

**What Mutect2 uses instead.** Its own Host annotation family, all of it in `mutect2_tool.cpp`, all of
it owner-selecting, and all of it built on the "active-core priority + `candidate_index`" idiom:

* `somatic_annotation_counts` → FORMAT/F1R2, F2R1, SB, and INFO/MBQ, MFRL, MMQ, MPOS (#5, **A**);
* `emitted_somatic_fragment_depth` → FORMAT/FAD (#3, **A**);
* `emitted_somatic_read_depth` → FORMAT/AD, DP (#7, **B** by CONSEQUENCE-1);
* `coverage_depth_for_candidate` → INFO/DP (#4, **A**);
* `calculate_somatic` / `calculate_posterior` flatten reducers → INFO/TLOD, NALOD, NLOD, POPAF and the
  emission gate (#8, #10, **A**);
* `apply_grouped_somatic_depth` → the FAD fallback scalars (#2, **A**).

**Can that path consume a sentinel owner?** Yes, in principle, but only under R3/R4, and *not* through
the same door as HC. Two differences matter:

1. Mutect2's probes are **not** pure index-bounds checks: `candidate_index` resolves a duplicated
   identity to the event-map-owned row, and every loop prefers the owner whose own active core
   contains the locus. Since (R1) the request for an allele is issued only by the region that owns its
   active core, and (R3) cores are disjoint, the preferred owner is normally the requested owner.
   Mutect2 therefore has, in code, the mitigation that HC's `annotations_from_owner` lacked.
2. But the mitigation is **ownership-based, not row-based** (R2): a sentinel entry can carry
   `likelihood_candidate_event_map_owned == 1`, so nothing in the probe *proves* the owner has rows.
   And the (A) sites all fail to reject the all-zero result once such an owner is chosen — the
   zero table is well-formed and `has_value()`-true. The defect family is therefore live in Mutect2,
   just harder to reach than in HC.

---

## 4. Empirical work

All commands below ran inside single `bash` invocations with `timeout` and scratch in a temporary
directory. Native invocations used `--native-pair-hmm-threads 4 --add-output-vcf-command-line false`;
the GATK side used
`third_party/jdk17/bin/java -Xmx1g -jar third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar Mutect2 … --add-output-vcf-command-line false`
(no `--native-pair-hmm-threads`, which is an HC-only flag). Two fixtures were used:

* **HCC1143 matched normal** — `testdata/real/cnv_somatic/chr20/HCC1143_{tumor,normal}.bam`,
  reference `testdata/real/cnv_somatic/human_g1k_v37.chr-20.truncated.fasta`, samples
  `HCC1143` / `HCC1143 BL`; contig `20` is truncated to 1,000,000 bp and yields 44 data rows.
* **DREAM synthetic** — `testdata/real/dream_synthetic/chr20/{tumor,normal}.bam`, reference
  `testdata/downloads/reference/hs37d5.fa.gz`, samples `synthetic.challenge.set1.tumor` /
  `synthetic.challenge.set1.normal`; 5 data rows in the audited window.

### M1 — whole-contig parity (the strongest measurement)

`-L 20` on HCC1143, native vs pinned GATK, full data-row comparison (all 10 columns, no
normalization): **44 GATK rows, 44 native rows, 0 rows differing in any field.** Times: native 72 s,
GATK 8.6 s.

### M2 — owner scale actually exercised

Same run with `FASTGATK_DEBUG_REGION_SCHEDULING=1`: the flattening reducer reports
`regions=590` (tumour owners) / `regions=1628` (the second, joint/sample pass) over
`candidates=1154` / `candidates=3366`. So the owner-selection machinery is exercised at **hundreds of
owners per invocation** on this fixture, and M1 shows it still matched GATK exactly.

### M3 — `-L` start invariance, and the one window that legitimately differs

| run | windows | result |
| --- | --- | --- |
| DREAM chr20 | `10000000`, `10010000`, `10016000`, `10029999`, `10080000` (all ending `-10500000`) | 5 record identities, **every overlapping row byte-identical** across the windows that contain it |
| HCC1143 | `20`, `20:50000-950000`, `20:150000-800000`, `20:300000-1000000` | 17 common rows, **0 differing** |
| HCC1143 cluster 159k | `150000-165000`, `155000-165000`, `158000-165000`, `159300-165000`, `159395-165000` | `159395` **does** differ from the others (drops the 159459 record; `ECNT/ECNTH 4→3`; `TLOD 6.46/6.38/4.63 → 4.44/4.33/4.13`; phasing `PGT/PID/PS` lost; `AD 17,2→16,2`; `AS_SB_TABLE 19,8→19,7`) |
| HCC1143 cluster 231k | `230000-240000`, `231000-240000`, `231300-240000` | 4 records, **0 differing** |

The 159k difference is **not** a defect: running pinned GATK on the same five windows gives
**exactly the same rows with exactly the same fields in every window** — GATK performs the same
record-set and TLOD change when `-L` starts at 159395 (it clips the phasing cluster). Native matched
GATK 0-differing on all five windows (`150000-165000`: 5/5 rows, `159300-165000`: 5/5,
`159395-165000`: 4/4, `230000-240000`: 4/4, `231300-240000`: 4/4).

### M4 — a structurally different region partition

`--force-active true` on `20:150000-165000` (all reads active ⇒ a very different AssemblyRegion
schedule): 5 GATK rows, 5 native rows, **0 differing**.

### M5 — direct probe of the (A) mechanism: was any published record annotated from a sentinel owner?

`somatic_annotation_counts_local` prints one `[FASTGATK_SOMATIC_READ]`/`[FASTGATK_SOMATIC_FRAGMENT]`
line per retained read/fragment **only if the locus group survived the context-ordinal filter**
(`3761-3790`) and the values are informative (`gatk_mutect_best_allele` returns a value). A sentinel
owner therefore produces **zero** lines for its locus. Running the whole HCC1143 contig with
`FASTGATK_DEBUG_SOMATIC_ANNOTATIONS=1`:

* published rows: 44; positions with **zero** annotation-evidence debug lines: **none** (each of the
  44 loci produced 49–82 lines; 1369 `FRAGMENT` + 1494 `READ` lines in total);
* positions where **any** debug line reported `ref <= -1.0e299` (the missing-value sentinel): **none**;
* published records with all-zero `SB` (`0,0,0,0`) or all-zero `F1R2`: **none**.

So on this fixture the priority-2 owner that the code selected was, for every published record, an
owner that *did* hold finite likelihood rows for the published allele.

### M6 — streaming path

`--stream-by-region` on the whole 1 Mb HCC1143 contig **aborts** for every tile size tried
(100000, 50000, 20000, 5000): `error: RESOURCE_EXHAUSTED: single-base Mutect2 tile exceeds safe
memory budget` (exit 2). On small windows where it does run (10 kb window, `--stream-by-region 1000`)
its rows were identical to the non-streaming run. This is an unrelated robustness limitation, not part
of the audited family; it is recorded because it restricts how far the tiled path could be probed.

### Clean negative — stated plainly

**I could not provoke the D2 family in Mutect2.** Across the two fixtures, 44 + 5 + 4 + 4 published
records, five `-L` starts, a `--force-active` schedule, a 590-owner flattened result, and a
per-locus debug audit of owner selection, native Mutect2 was byte-identical to pinned GATK and no
published locus showed the sentinel signature (zero evidence lines, missing `ref`, or all-zero
`SB`/`F1R2`). The (A) classifications in §2 are reasoned from code; the negative is measured.

**What this negative does not cover** (honesty limits): the fixtures are a 1 Mb truncated chr20 and a
5-locus DREAM window; only tumor/normal matched mode; `-ERC GVCF` (which additionally enables
`materialize_reference_confidence_eventmap_calls` and the `<NON_REF>` writer) was **not** exercised;
`--stream-by-region` could not be exercised at contig scale (M6); `/tmp` non-persistence prevented
multi-invocation reuse, so every comparison was recomputed in one command.

---

## 5. Ranked (A) list and the cheapest next test for each

Ranking is by published-field impact × likelihood of the priority-2 owner being sentinel. "Cheapest
next test" assumes no code edits and ≤ ~2 minutes per command.

| rank | site | published effect if hit | cheapest next test |
| --- | --- | --- | --- |
| 1 | **#5** `somatic_annotation_counts` owner loop `4116` | FORMAT/F1R2, F2R1, SB (+ MBQ/MFRL/MMQ/MPOS) become `0`/`0,0,0,0` for tumour, every normal and every extra tumour sample; propagates into `FilterMutectCalls`' orientation/strand models | **Falsifier run (already scripted above, ~75 s):** whole-contig HCC1143 with `FASTGATK_DEBUG_SOMATIC_ANNOTATIONS=1`, then assert every published position has ≥1 debug line, no line reports `ref<=-1e299`, and no record publishes `SB=0,0,0,0` or `F1R2=0,0`. Extend it to a **locus set with core-boundary-crossing indels** (the R4 case) and to `-ERC GVCF`, where halo alleles are materialized more often. |
| 2 | **#3** `emitted_somatic_fragment_depth` owner loop `3323` | FORMAT/FAD = `0,0,…` (tumour and normal) while AD/DP stay non-zero — an internally inconsistent row that no other check would catch | Same run; add a cross-check that `FAD` is never all-zero on a record whose `F1R2+F2R1` (or AD alt) is non-zero, and that native `FAD` equals GATK's. |
| 3 | **#4** `coverage_depth_for_candidate` loop `3460` | INFO/DP too low (`0` or a stale scalar), and `-L`-start dependent because equal-priority owners are first-wins | Compare INFO/DP for the same record across two `-L` starts (**done for 17 rows: identical**) — repeat on windows whose core boundary cuts an indel's reference span, and on `--force-active` windows. |
| 4 | **#8** `calculate_somatic` flatten reducer `4463` | INFO/TLOD (and NALOD/NLOD/AF) overwritten with the kernel's initialised `0.0` for a sentinel entry; since TLOD gates emission the symptom is a **missing or invented record** | Row-**set** parity against GATK (missing/extra rows, not just fields) over a locus set with the documented raw/EventMap duplicate identities — a multi-ALT locus plus a low-quality-recovery call. Row-set parity is the cheapest canary: **done on 44 rows + 5 windows + force-active: zero missing/extra rows.** |
| 5 | **#10** `calculate_posterior` flatten reducer `5316` | somatic/artifact/germline posterior defaults (0 / 1) published for a locus → `FilterMutectCalls` germline/artifact decisions change | Same row-set parity, then `FilterMutectCalls` (native vs pinned GATK) on the same VCF and compare `FILTER`/`AS_FilterStatus`. |
| 6 | **#2** `apply_grouped_somatic_depth` loop `3244` | the FAD *fallback* scalars only (used when the emitted table is unavailable) — smallest blast radius | Cover it for free with the rank-2 FAD cross-check on a **normal** view where the primary table is absent. |

Note on `filter-mutect-calls`: `fastgatk-native/src/filter_mutect_tool.cpp` contains **no**
`assembly_region_likelihood_results` reference at all (grep over `fastgatk-native/src`: only
`mutect2_tool.cpp`, `calling_pipeline.cpp`, `hc_call.cpp`). `FilterMutectCalls` is a VCF-in tool and
has no owner-selection sites of this family; it can only inherit the defect through bad
`F1R2/F2R1/SB/TLOD` values written by `Mutect2` (which is why #5 and #8 rank first).

---

## 6. Reasoned vs measured (explicit)

**Reasoned (code reading only, no execution):**

* R1 (single requesting region per candidate), R2 (`event_map_owned` ≠ requested), R3 (priority rule
  selects the requested owner because cores are disjoint and a part's `calling_regions` are its own
  unpadded core), R4 (residual hole; the second variant for the merge-style reducers).
* The class (A)/(B) of every row of the §2 table, and every "reason" column.
* The `calculate_output_variant_annotations` claim (§3) — a grep/enumeration fact, but obtained by
  reading, not by a runtime assertion.
* That `fastgatk-kernels/src/somatic.cpp:795` initialises the sparse TLOD to `0.0` (so a
  never-written entry surfaces as `0`, not `-inf`).
* Dead-code verdict for #6 and unreachability verdict for #9 (derived from the complete call-site
  list, obtained by grep; no dynamic coverage run was performed).

**Measured (executed in this round, exact results in §4):**

* M1 native==GATK 44/44 rows, all fields, whole 1 Mb HCC1143 contig (0 differing).
* M2 `regions=590` / `candidates=1154` (and `1628`/`3366`) from
  `FASTGATK_DEBUG_REGION_SCHEDULING` — the owner machinery runs at scale.
* M3 `-L`-start sweeps: DREAM 5 windows and HCC1143 4 windows identical; the 159395 window's
  legitimate change reproduced identically by pinned GATK (0 field differences on 5/5, 4/4, 5/5, 4/4,
  4/4 rows).
* M4 `--force-active true` parity (5/5 rows).
* M5 per-locus owner audit with `FASTGATK_DEBUG_SOMATIC_ANNOTATIONS=1`: 44/44 published loci
  produced annotation evidence; 0 loci with a missing (`-1.0e299`) REF row; 0 records with all-zero
  `SB`/`F1R2`.
* M6 `--stream-by-region` aborts with `RESOURCE_EXHAUSTED` on the 1 Mb contig for sizes
  100000/50000/20000/5000; identical rows on a 10 kb window at size 1000.

**Not measured / not verified:** `-ERC GVCF` (BP_RESOLUTION/GVCF) annotation paths; any fixture with
a core-boundary-crossing indel at the locus of interest (the R4 trigger); multi-tumor-tag and
normal-only configurations; the `serial` backend (not used — the defect would be backend-independent,
and the point of this round was reachability, not backend equivalence). No claim in this report is
based on an unobserved run.
