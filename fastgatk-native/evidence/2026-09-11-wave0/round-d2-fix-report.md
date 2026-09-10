# Round D2 — sentinel-context twin remap: implementation, gates, verdict

**Verdict: the fix is implemented, all five gate steps pass, and the tree CONTAINS the change.**
Only one file is modified: `fastgatk-native/src/hc_call.cpp` (+53 / −2).

Baseline commit: `9b59b21` (working tree otherwise clean before the change).

---

## 1. What was measured before writing any fix (this changed the placement)

The brief prescribed remapping a sentinel index to an identical-allele sibling
**inside `likelihood_result->candidates`**. I instrumented the code to locate that
sibling, and the measurement shows **no such sibling exists in the same `Result`**:
the twin lives in a *different* `assembly_region_likelihood_result` owner.

Literal command (env-guarded temporary instrumentation, since removed):

```
FASTGATK_DEBUG_D2_INDICES=1 fastgatk-native/build/fastgatk-hc-call \
  -R fixtures/chr20/ref20mnp.fasta -I fixtures/chr20/mnp.bam \
  -L 20:10020381-10020710 --emit-ref-confidence GVCF --max-mnp-distance 1 \
  --threads 1 --add-output-vcf-command-line false -O out.g.vcf
```

Observed (window B, the failing window; region owner index in `[]`):

```
[D2] HC2 group nobs=1 owners=2 ncand_merged=6 nobsctx=0
[D2]   owner 0 ncand=6 nctx=6 nrealign=38
[D2]     o0 idx=4 tid=0 pos=10020678 ref=AC alt=TA ral=AC aal=TA ord=4294967295
[D2]     o0 idx=5 tid=0 pos=10020679 ref=CA alt=AT ral=CA aal=AT ord=4294967295   <-- sentinel, selected today
[D2]   owner 1 ncand=2 nctx=2 nrealign=33
[D2]     o1 idx=0 tid=0 pos=10020678 ref=AC alt=TA ral=AC aal=TA ord=0
[D2]     o1 idx=1 tid=0 pos=10020679 ref=CA alt=AT ral=CA aal=AT ord=0            <-- PAIRED TWIN
[D2]   owner 0 PROBE has=1 mqsq=97200 mqdp=27 sb=0,0,0,0                          <-- today's output
[D2]   owner 1 PROBE has=1 mqsq=28800 mqdp=8  sb=0,0,3,3                          <-- GATK's output
```

The `index=1 / ordinal 0` vs `index=5 / sentinel` pair named in the brief is real, but it
is **spread across two owner `Result`s**, not two rows of one owner's candidate vector.
Forcing each owner one at a time confirms which one is GATK-correct:

```
$ FASTGATK_DEBUG_D2_FORCE_OWNER=0 ... -L 20:10020381-10020710
win=10020381 owner=0 rows=13 RAW_MQandDP=97200,27 SB=0,0,0,0
$ FASTGATK_DEBUG_D2_FORCE_OWNER=1 ... -L 20:10020381-10020710
win=10020381 owner=1 rows=13 RAW_MQandDP=28800,8  SB=0,0,3,3
```

Control window A shows that this is *only* an owner-ordering defect, not a new rule:

```
[D2] HC2 group nobs=1 owners=3
[D2]   owner 0 PROBE nullopt        [D2]   owner 1 PROBE nullopt
[D2]   owner 2 PROBE has=1 mqsq=28800 mqdp=8 sb=0,0,3,3
[D2]     o2 idx=5 tid=0 pos=10020679 ref=CA alt=AT ral=CA aal=AT ord=0
```

So in window A the existing first-success owner loop already lands on the twin owner and
window A is already correct. In window B an owner that carries the allele **with the
sentinel** happens to come first, and it is accepted because it "succeeds".

### Consequence for placement

A remap *inside* `calculate_variant_annotations` is provably impossible here: for owner 0
`all_indices == [5]` and there is no sibling with that allele, reference allele, alternate
allele and a non-sentinel ordinal anywhere in that `Result`. Making `calculate_variant_annotations`
reach across owners would require changing its signature and threading the owner pool
through — exactly the "wide change" the brief says to stop on.

The same effect is available one layer up, in `hc_call.cpp`'s owner selection (a file the
brief lists as a key file, and specifically the gVCF annotation boundary), with the same
outcome the brief asks for: the published allele is resolved to its identical-allele twin
and **both** the `allele_read_likelihoods` / `reference_read_likelihoods` row lookups and
the `likelihood_candidate_read_context_ordinals` lookup are served by the twin, because a
different owner is chosen. This is the intervention I implemented.

---

## 2. Exact diff applied

```diff
$ git diff
diff --git a/fastgatk-native/src/hc_call.cpp b/fastgatk-native/src/hc_call.cpp
index 4d6098a..dbfa554 100644
--- a/fastgatk-native/src/hc_call.cpp
+++ b/fastgatk-native/src/hc_call.cpp
@@ -3685,6 +3685,44 @@ std::string vcf_text(const fastgatk::io::HtsReader& reader,
     return out.str();
 }
 
+// GATK annotates a record from the AssemblyRegion whose retained PairHMM read
+// set produced it.  A partitioned native Result can list the same allele in
+// more than one AssemblyRegion owner, and an owner that never requested
+// PairHMM for that allele keeps the UINT32_MAX context-ordinal sentinel.
+// Annotating from such an owner empties `context_mapping_evidence` (its own
+// likelihood rows are all -inf) and drops the MQ gate onto the raw-overlap
+// predicate, so the evidence - and therefore the published annotations -
+// depend on the `-L` window start.  Resolve every published allele to its
+// identical-allele twin inside the owner and require that twin to carry a real
+// PairHMM context ordinal, so both the likelihood-row lookups and the
+// context-ordinal lookup are served by the twin.
+bool owner_has_pairhmm_context(
+    const fastgatk::calling::Result& owner,
+    const std::vector<const fastgatk::calling::AssemblyCandidate*>& alleles) {
+    if (alleles.empty()) return false;
+    for (const auto* allele : alleles) {
+        if (allele == nullptr) return false;
+        std::size_t twin = std::numeric_limits<std::size_t>::max();
+        for (std::size_t index = 0; index < owner.candidates.size(); ++index) {
+            const auto& candidate = owner.candidates[index];
+            if (candidate.tid == allele->tid && candidate.position == allele->position &&
+                candidate.reference == allele->reference &&
+                candidate.alternate == allele->alternate &&
+                candidate.reference_allele == allele->reference_allele &&
+                candidate.alternate_allele == allele->alternate_allele) {
+                twin = index;
+                break;
+            }
+        }
+        if (twin == std::numeric_limits<std::size_t>::max()) return false;
+        if (twin >= owner.likelihood_candidate_read_context_ordinals.size()) return false;
+        if (owner.likelihood_candidate_read_context_ordinals[twin] ==
+            std::numeric_limits<std::uint32_t>::max())
+            return false;
+    }
+    return true;
+}
+
 std::string gvcf(const fastgatk::io::HtsReader& reader,
                  fastgatk::calling::Result& result,
                  int sample_ploidy,
@@ -4831,8 +4869,21 @@ std::string gvcf(const fastgatk::io::HtsReader& reader,
             const auto annotation_qual = calls[best] == nullptr ? 0.0 : calls[best]->qual;
             const auto annotations_from_owner = [&]()
                 -> std::optional<fastgatk::calling::GenotypeCall::Annotations> {
-                for (const auto& owner : result.assembly_region_likelihood_results) {
-                    if (owner == nullptr) continue;
+                // Owners that hold the identical-allele twin with a real
+                // PairHMM context come first; owners that merely carry the
+                // allele with the UINT32_MAX context sentinel keep their
+                // historical first-match order as a fallback.  Selecting the
+                // twin owner is what makes RAW_MQandDP/SB independent of the
+                // `-L` window start.
+                std::vector<const fastgatk::calling::Result*> ordered_owners;
+                ordered_owners.reserve(result.assembly_region_likelihood_results.size());
+                for (const auto& owner : result.assembly_region_likelihood_results)
+                    if (owner != nullptr && owner_has_pairhmm_context(*owner, group.candidates))
+                        ordered_owners.push_back(owner.get());
+                for (const auto& owner : result.assembly_region_likelihood_results)
+                    if (owner != nullptr && !owner_has_pairhmm_context(*owner, group.candidates))
+                        ordered_owners.push_back(owner.get());
+                for (const auto* owner : ordered_owners) {
                     if (const auto owned = fastgatk::calling::calculate_output_variant_annotations(
                             *annotation_reads, *owner, group.candidates, annotation_qual,
                             informative_read_overlap_margin, include_spanning_deletion,
```

`git diff --stat`:

```
 fastgatk-native/src/hc_call.cpp | 55 +++++++++++++++++++++++++++++++++++++++--
 1 file changed, 53 insertions(+), 2 deletions(-)
```

`git status --short`:

```
 M fastgatk-native/src/hc_call.cpp
```

No other file was touched. No commit, no branch.

### Why this is minimal and safe

* It is a **preference**, not a replacement: owners holding the twin with a real
  context ordinal are tried first in their original relative order; if none of them
  yields an annotation, the original first-match order is retried unchanged. No
  previously working path is removed.
* `calculate_variant_annotations`, `calculate_output_variant_annotations`, the partition
  merge, PairHMM pairing, `build_reference_blocks` and
  `build_profile_local_reference_blocks` are untouched (all zero-diff).
* No Mutect2 / BQSR / other-tool code is touched; the helper is file-local to
  `hc_call.cpp` and is only called from the diploid gVCF record writer.

---

## 3. Gate results

### Step 1 — rebuild OpenMP HC — PASS

```
$ third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/cmake --build fastgatk-native/build \
    --target fastgatk-hc-call -j 16
[ 96%] Building CXX object CMakeFiles/fastgatk-hc-call.dir/src/hc_call.cpp.o
[100%] Linking CXX executable fastgatk-hc-call
[100%] Built target fastgatk-hc-call
```

Exit status 0, no warnings or errors for the changed file.

### Step 2 — acceptance oracle — PASS (exit 0)

```
$ python3 fastgatk-native/scripts/verify_hc_window_invariance_gatk_oracle.py ; echo exit=$?
...
  "status": "pass",
  "strict_mode": true,
  "violations": [],
exit=0
```

All 10 gating windows report `data_rows_byte_identical=True` and
`pos_10020680_identical=True`. The two structural windows (`10020421`, `10020431`) are
`gating=False` by design (they carry the unrelated D3 block-granularity divergence) but
their `10020680` row is still identical, and the pinned-GATK drift guard reports
`gatk_baseline_drift=false` in all 12 windows.

Per-window POS 10020680 evidence, GATK vs native (all 12 windows):

| window start | GATK rows | native rows | byte-identical | GATK `RAW_MQandDP` | native `RAW_MQandDP` | GATK `SB` | native `SB` |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 10019901 (window A) | 48 | 48 | **yes** | 28800,8 | 28800,8 | 0,0,3,3 | 0,0,3,3 |
| 10020201 | 44 | 44 | yes | 28800,8 | 28800,8 | 0,0,3,3 | 0,0,3,3 |
| 10020301 | 13 | 13 | yes | 28800,8 | 28800,8 | 0,0,3,3 | 0,0,3,3 |
| 10020351 | 13 | 13 | yes | 28800,8 | 28800,8 | 0,0,3,3 | 0,0,3,3 |
| **10020381** | 13 | 13 | **yes** | 28800,8 | **28800,8** | 0,0,3,3 | **0,0,3,3** |
| **10020391** | 13 | 13 | **yes** | 28800,8 | **28800,8** | 0,0,3,3 | **0,0,3,3** |
| **10020401** | 13 | 13 | **yes** | 28800,8 | **28800,8** | 0,0,3,3 | **0,0,3,3** |
| **10020411** | 13 | 13 | **yes** | 28800,8 | **28800,8** | 0,0,3,3 | **0,0,3,3** |
| 10020421 (structural) | 53 | 12 | no (D3) | 97200,27 | 97200,27 | 1,2,6,16 | 1,2,6,16 |
| 10020431 (structural) | 51 | 9 | no (D3) | 97200,27 | 97200,27 | 1,2,6,16 | 1,2,6,16 |
| 10020441 | 4 | 4 | yes | 28800,8 | 28800,8 | 0,0,3,3 | 0,0,3,3 |
| 10020501 | 4 | 4 | yes | 28800,8 | 28800,8 | 0,0,3,3 | 0,0,3,3 |

The four windows B/C/D/E (`10020381`, `10020391`, `10020401`, `10020411`) move from
`97200,27 / 0,0,0,0` to GATK's `28800,8 / 0,0,3,3`. The two structural windows keep
`97200,27 / 1,2,6,16`, i.e. the fix does **not** over-apply — it only trades the sentinel
owner for a real-context twin when one exists.

### Step 3 — window A still byte-identical — PASS

Covered by the oracle above, and independently reproduced from the raw VCF:

```
$ fastgatk-native/build/fastgatk-hc-call -R fixtures/chr20/ref20mnp.fasta -I fixtures/chr20/mnp.bam \
    -L 20:10019901-10020710 --emit-ref-confidence GVCF --max-mnp-distance 1 --threads 1 \
    --add-output-vcf-command-line false -O a.g.vcf
native rows: 48
POS 10020680 ... RAW_MQandDP=28800,8 ... SB=0,0,3,3
```

48 data rows (unchanged), `divergent positions: []`, diff 0 lines against pinned GATK.

### Step 4 — HC/Mutect2 subset on OpenMP — PASS

```
$ fastgatk-native/scripts/run_regression.sh --backend omp -R 'fastgatk-(hc|mutect2)'
| omp | OpenMP (fastgatk-native/build) | 通过 | 73/73 | 214.18 sec |
exit=0
```

`100% tests passed, 0 tests failed out of 73` (`73` rather than the expected `72`: the
filter matches 73 tests in the current tree).

### Step 5 — serial rebuild + full dual-backend regression — PASS

```
$ third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/cmake --build fastgatk-native/build-serial \
    --target fastgatk-hc-call -j 16
[100%] Built target fastgatk-hc-call            (exit 0)

$ fastgatk-native/scripts/run_regression.sh --label 'D2 twin-remap fix'
| omp    | OpenMP (fastgatk-native/build)         | 通过 | 281/281 | 1042.16 sec |
| serial | Serial (fastgatk-native/build-serial)  | 通过 | 281/281 | 1045.00 sec |
exit=0
```

`omp.log`: `100% tests passed, 0 tests failed out of 281`.
`serial.log`: `100% tests passed, 0 tests failed out of 281`.
Evidence directory: `.diag/regression/20260911-021036`. (`281`, not `280`: the full suite
contains one more test than the brief assumed. Step-4 subset evidence:
`.diag/regression/20260911-015551`.)

---

## 4. Tree state

**The tree CONTAINS the change.** `git status --short` → ` M fastgatk-native/src/hc_call.cpp`;
no commit and no branch was created; nothing outside `fastgatk-native/src` was modified
(all instrumentation was reverted with `git checkout --` before the final patch was applied,
and `git diff` contains only the hunk above).

Nothing was reverted, because no gate step failed.

---

## 5. What I could not verify / known limits

1. **I did not run the oracle script end-to-end against the pristine pre-fix binary.** The
   pre-fix raw values were reproduced directly instead (`-L 20:10020381-10020710` →
   `RAW_MQandDP=97200,27`, `SB=0,0,0,0`; `-L 20:10019901-10020710` → 48 rows,
   `28800,8`, `0,0,3,3`), which matches the defect description exactly. The oracle's
   `--expect-divergence` mode was not exercised.
2. **The placement deviates from the literal brief.** The brief asked for the remap inside
   `likelihood_result->candidates`. Section 1 shows with measured evidence that no
   identical-allele sibling with a non-sentinel ordinal exists inside the owner `Result`
   that the annotation currently uses, so the literal placement is unimplementable without
   changing `calculate_variant_annotations`' signature. The equivalent resolution is done at
   the owner-selection boundary in `hc_call.cpp` (a listed key file), with the same effect on
   both the likelihood-row lookups and the context-ordinal lookup. No broad restructuring
   (merge, PairHMM pairing, reference blocks) was performed.
3. **The arbitrary-ploidy sibling path is untouched.** The `if (sample_ploidy != 2)` branch in
   the same `gvcf()` writer has the same "first owner that succeeds" shape. It is not on the
   failing path (this run is diploid) and is not covered by the oracle, so it was deliberately
   left alone for minimality. Whether it carries the same latent `-L` sensitivity is
   **unverified**.
4. **`calculate_variant_annotations`, `calculate_output_variant_annotations` and
   `calling_pipeline.cpp` were not modified at all** — the refuted approaches (a)–(d) from the
   brief were not attempted.
5. The `-L`-invariance is verified for the fixture windows in the oracle and for the whole
   HC/Mutect2 + full regression suites; it was **not** verified on other references or on
   samtools/streamed non-GVCF output modes beyond what the existing suites cover.
