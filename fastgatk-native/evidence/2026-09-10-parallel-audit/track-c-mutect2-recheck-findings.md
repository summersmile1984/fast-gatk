# Track C — Mutect2 independent re-check (findings)

Date: 2026-09-11. Scope: audit only; no production source touched.
Binaries: `fastgatk-native/build/fastgatk-mutect2` (Kokkos/OpenMP),
`fastgatk-native/build-serial/fastgatk-mutect2` (Serial). Both current with sources
(`calling_pipeline.cpp` mtime 2026-09-10 17:03, binaries 2026-09-10 18:17).
Oracle: pinned GATK 4.6.2.0 via `third_party/jdk17/bin/java -jar
third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar Mutect2`.

---

## (a) Exact commands and observed results

| # | Command | OMP `build/` | Serial `build-serial/` |
|---|---------|--------------|------------------------|
| 1 | `FASTGATK_REQUIRE_GATK_ORACLE=1 timeout 240 python3 fastgatk-native/scripts/verify_mutect2_recheck_normal_replay.py` | **PASS** (exit 0) | **PASS** (exit 0, `FASTGATK_MUTECT2_BINARY=…/build-serial/…`) |
| 2 | `FASTGATK_REQUIRE_GATK_ORACLE=1 timeout 240 python3 fastgatk-native/scripts/verify_mutect2_recheck_assembly_resultset_joint.py` | **PASS** (exit 0) — after the rewrite in (d) | **PASS** (exit 0) — after the rewrite in (d) |
| 2′ | same joint script, *original* revision, unmodified | **BROKEN** (exit 3) | not run (defect is backend-independent) |
| 3 | `third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest --test-dir fastgatk-native/build -R 'fastgatk-mutect2' -j 8` | **30/30 passed**, 144.8 s | — |
| 4 | same with `--test-dir fastgatk-native/build-serial` | — | **30/30 passed**, 142.2 s |
| 5 | `python3 fastgatk-native/scripts/verify_mutect2_tlod_formula.py` (both backends) | PASS, `maximum_tlod_formula_delta=4.6896e-13` | PASS, identical per-locus deltas |
| 6 | TLOD window-context probe (`-L 17:69368-69368`, `17:69200-69500`, `17:69000-70000`), GATK vs native | see (c) | see (c) |

Observed JSON summaries (abridged, both backends identical unless stated):

Script 1 (`normal_replay`, DREAM synthetic chr20 tumor/normal, whole contig `-L 20`):
```
{"status":"pass","gatk_call_count":5,"native_call_count":5,
 "gatk_native_allele_recall":1.0,"gatk_native_allele_precision":1.0,
 "missing_alleles":[],"extra_alleles":[],"row_difference_count":0,
 "tlod_drift_count":0,"schema_header_lines_match":true,"tbi_present":true,
 "manifest_tumor_sample":"synthetic.challenge.set1.tumor",
 "manifest_normal_sample":"synthetic.challenge.set1.normal"}
```
Run 1 on `build/` used the script's default binary; run on `build-serial/` reported
`"native_binary":"…/build-serial/fastgatk-mutect2"` with the same numbers.

Script 2 (`assembly_resultset_joint`, chr17 69k–70k + HCC1143 matched-normal leg):
`"native_tumor_only_matches_gatk":true`, `"native_joint_matches_gatk":true`,
`"manifest_ok"` satisfied, `"real_matched_normal_joint":{"available":true,
"interval":"20:67000-69000","native_allele_count":1,"gatk_allele_count":1,
"native_alleles":[["20",68037,"G","T"]],"gatk_alleles":[["20",68037,"G","T"]],
"allele_sets_equal":true,"field_differences":[]}`, `"status":"pass"`.
Cross-mode observation on the chr17 fixture: `observation_gatk_tumor_only_eq_joint=false`,
`observation_native_tumor_only_eq_joint=false` (both engines; see (c)).

ctest: `100% tests passed, 0 tests failed out of 30` on **both** trees.

---

## (b) Verdict on the regression question

**Not proven to have regressed; no regression observed.** Positive evidence, all measured:

1. 30/30 Mutect2 ctest tests pass on both backends in the current trees (~145 s each).
2. The DREAM synthetic tumor/normal full-contig chr20 run is byte-equivalent to pinned GATK
   on all 5 data rows: 0 missing alleles, 0 extra alleles, 0 row-field differences
   (QUAL/FILTER/INFO/FORMAT/sample all equal), and the schema headers match after stripping
   `##GATKCommandLine`.
3. chr17 69k–70k: native tumor-only emits exactly the 6 GATK alleles, and the joint
   (tumor + renamed mirrored normal) run emits exactly the same rows as GATK joint (both 0
   rows — see the caveat below).
4. Real matched-normal joint run (HCC1143 chr20, `20:67000-69000`): 1 allele each side,
   identical identity `20:68037 G>T`, and **zero** per-field data-row differences.
5. The two backends agree with each other bit-for-bit on every TLOD probed (below).

Honest limits of this verdict: (2)–(4) are equivalence checks on three small fixtures, not a
whole-genome sweep; they can detect a regression that perturbs these rows, not one confined to
loci absent from these fixtures. The chr17 joint leg is the weakest link — it compares 0 rows
to 0 rows, so the chr17 joint evidence is "no divergence in a degenerate case"; the real
matched-normal leg (HCC1143) is the substantive joint evidence.

---

## (c) Residual divergences

### c.1 Data-row divergences: none found

On the DREAM chr20 tumor/normal fixture the ordered data rows matched field-for-field
(0 diffs). On HCC1143 `20:67000-69000` joint mode likewise 0 field diffs. So on the audited
fixtures there is no "first differing field" to report — that is the finding.

### c.2 The previously recorded Δ3.5 TLOD residual at chr17:69368 is NOT reproducible now

Minimal repro (120 s cap each; `REF=gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta`,
`BAM=gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam`):

```
third_party/jdk17/bin/java -Xmx2g -jar third_party/gatk-package/gatk-4.6.2.0/…-local.jar \
   Mutect2 -R $REF -I $BAM --tumor-sample NA12878 -L 17:69000-70000 -O gatk.vcf.gz
fastgatk-native/build/fastgatk-mutect2 -R $REF -I $BAM --tumor-sample NA12878 \
   -L 17:69000-70000 -O native.vcf.gz
```

Observed (2-dp VCF TLOD, whole window):

| locus | GATK TLOD | native OMP | native Serial | Δ |
|-------|-----------|------------|---------------|---|
| 69067 | 6.68 | 6.68 | 6.68 | 0 |
| 69124 | 3.28 | 3.28 | 3.28 | 0 |
| **69368** | **70.51** | **70.51** | **70.51** | **0** |
| 69631 | 67.69 | 67.69 | 67.69 | 0 |
| 69803 | 3.72 | 3.72 | 3.72 | 0 |
| 69807 | 3.72 | 3.72 | 3.72 | 0 |

`AS_SB_TABLE` also matches GATK exactly at 69368 (`20,0|22,1`), i.e. the previously documented
native strand-bias defect (recorded as `21,0|27,1`, sum 49 ≠ DP 43) is **not present** in the
current binary.

Window-context sensitivity (the discriminator used in the earlier R15 round; previously
"GATK drifts with window, native constant 74.01"):

| interval | GATK @69368 | native @69368 |
|----------|-------------|---------------|
| `17:69368-69368` | 67.43 | **67.43** |
| `17:69200-69500` | 70.51 | **70.51** |
| `17:69000-70000` | 70.51 | **70.51** |

Native now reproduces GATK's window-dependent drift. Full-precision native echo
(`FASTGATK_DEBUG_TLOD=1`, lines tagged `[FASTGATK_TLOD_FULL]`) at chr17:69368 is
`70.514741250854556` on **both** backends (bit-identical).

Source-level interpretation: the previously hypothesised cause ("assembly/haplotype input
parity, A1/A-line boundary") is no longer exercised — whatever upstream input difference
produced 74.01 is not present in the current binaries. This audit cannot say *which* commit
changed it (the repo carries a single baseline commit plus doc commits, `git status` clean),
so the change is not attributable from history alone — only the current state is measured.

**Boundary status, stated precisely:** native ≈ GATK to within GATK's VCF rounding
(≤ 0.005 TLOD at every probed locus, including the context-sensitive 69368); machine-precision
GATK-side agreement is **not** established, because GATK's VCF serialises TLOD to 2 decimals
and no GATK-side full-precision channel was used here.

### c.3 Engine-parity figure does not reproduce as recorded — flag for the orchestrator

`A_B_EXECUTION.md` (R24) and `SESSION_HANDOFF.md` record engine parity **max Δ = 1.4e-14**
after the `[FASTGATK_TLOD_FULL]` channel was added. Re-running the same verifier now gives, on
**both** backends:

```
verify_mutect2_tlod_formula.py  →  {"loci_checked":6,
  "maximum_tlod_formula_delta":4.689582056016661e-13,
  "mean_tlod_formula_delta":8.37108160567368e-14,
  "per_locus_deltas":[1.78e-15,1.33e-15,2.84e-14,4.69e-13,8.88e-16,8.88e-16],
  "status":"pass"}
```

i.e. the max delta is today **4.69e-13 at locus 17:69631** (and 2.84e-14 at 69368), not the
recorded 1.4e-14 / 0.0. Per-locus deltas differ from the recorded list at exactly the loci the
recorded list reported as 0.0 and 1.4e-14, so this is a change in the numeric path, not merely
a display artefact. Magnitude is still float-noise and far under the verifier's `<1e-8`
full-precision guard, so `fastgatk-mutect2-tlod-formula` passes — but the **recorded 1.4e-14
headline figure is stale/inaccurate as of these binaries**.

Also worth correcting in the shared record: **`FASTGATK_TLOD_FULL` is not an environment
switch**. Nothing in `fastgatk-native/src` or `fastgatk-kernels/src` reads it
(`grep -rn 'getenv("FASTGATK_TLOD_FULL")'` → no hits). The gate is `FASTGATK_DEBUG_TLOD=1`;
`FASTGATK_TLOD_FULL` is merely the log tag that `verify_mutect2_tlod_formula.py` parses
(`verify_mutect2_tlod_formula.py:220`).

---

## (d) Scripts

### d.1 `fastgatk-native/scripts/verify_mutect2_recheck_normal_replay.py` — UNCHANGED, works

* Pins: DREAM synthetic chr20 tumor/normal (`testdata/real/dream_synthetic/chr20/{tumor,normal}.bam`,
  `testdata/downloads/reference/hs37d5.fa.gz`, samples `synthetic.challenge.set1.{tumor,normal}`,
  `-L 20`) against pinned GATK Mutect2. Full data-row equivalence (allele recall/precision,
  QUAL/FILTER/INFO/FORMAT/sample), schema header equality minus `##GATKCommandLine`, `.tbi`
  presence, and the `--output-manifest` contract (`schema_version==1`, tumor/normal selection).
* Verdict: **correct and passing as written on both backends**; no edit made.
* Caveat: it skips (exit 0, `status: skip`) unless `FASTGATK_REQUIRE_GATK_ORACLE=1`, so wire
  that env var in CMake for it to be a real gate.

CMake registration snippet (**NOT applied**), modelled on the existing
`fastgatk-mutect2-tlod-formula` block in `fastgatk-native/CMakeLists.txt`:

```cmake
    add_test(NAME fastgatk-mutect2-recheck-normal-replay
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_recheck_normal_replay.py")
    set_tests_properties(fastgatk-mutect2-recheck-normal-replay PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1"
        TIMEOUT 300)
```
(Runtime observed ≈ 100 s wall — the GATK pass dominates; a 300 s TIMEOUT is advised.)

### d.2 `fastgatk-native/scripts/verify_mutect2_recheck_assembly_resultset_joint.py` — REWRITTEN

Original revision, as delivered, was **unrunnable** (verified by running it):

```
$ FASTGATK_REQUIRE_GATK_ORACLE=1 timeout 240 python3 \
    fastgatk-native/scripts/verify_mutect2_recheck_assembly_resultset_joint.py
native joint failed (rc=2):
  error: BAD_INPUT: requested sample is absent from tumor input headers: NA12878_NORMAL_MATCH
EXIT=3
```

Two independent defects:

1. **The matched normal was never supplied.** The joint runs passed
   `--normal-sample NA12878_NORMAL_MATCH` but only one `-I` (the tumor BAM). Native rejects it
   (`mutect2_tool.cpp` sample validation); GATK would have failed identically.
2. **The central assertion was invalid by construction.** The script required
   "tumor-only allele set == joint allele set" and asserted GATK honours it too. Measured with
   the normal it describes (a renamed mirror of the tumor reads):

   ```
   GATK   tumor-only  6 rows      GATK   joint  0 rows
   native tumor-only  6 rows      native joint  0 rows
   ```

   GATK itself drops every call in joint mode there, so `gatk_tumor_only_joint_match` could
   never be true. The premise "a matched normal that mirrors the tumor reads … must not change
   which candidate events the tumor reports" is false under Mutect2 semantics.

Fix applied (rewrite, ~290 lines; the fixture, the intent and the JSON-summary style are kept):

* builds the mirrored matched normal in the `TemporaryDirectory` with `pysam`
  (rewrites `@RG/SM` → `NA12878_NORMAL_MATCH` + index) and wires it as a second `-I` for
  native and a second `-I -normal` for GATK;
* gates on pinned-GATK parity in **both** modes on chr17 (allele sets) and on the
  `--output-manifest` contract; cross-mode equality is reported as
  `observation_*_tumor_only_eq_joint` and is **not** a gate;
* adds a real matched-normal joint leg (HCC1143 chr20 `20:67000-69000`,
  `testdata/real/cnv_somatic/…`), which is non-degenerate (1 allele both sides, 0 field diffs)
  and is gated on allele-set equality with per-field diffs reported;
* reports `joint_rows_zero_by_construction` so the degenerate chr17 joint result cannot be
  mistaken for strong evidence.

Both gates passed on both backends (see (a)).

CMake registration snippet (**NOT applied**):

```cmake
    add_test(NAME fastgatk-mutect2-recheck-assembly-resultset-joint
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_recheck_assembly_resultset_joint.py")
    set_tests_properties(fastgatk-mutect2-recheck-assembly-resultset-joint PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1"
        TIMEOUT 300)
```
(Depends on `pysam` being importable; it self-skips with `status: skip` otherwise.
Runtime observed ≈ 60–90 s.)

---

## (e) Proven vs not proven

**Proven (measured this session, on these binaries):**

* The shared calling-pipeline change has **not** regressed any of the Mutect2 evidence checked:
  30/30 ctest green on both backends; DREAM chr20 tumor/normal identical to GATK on all 5 rows
  and all fields; chr17 tumor-only allele set identical to GATK; HCC1143 real matched-normal
  joint identical to GATK (allele identity + 0 field diffs).
* The previously recorded Δ3.5 TLOD residual at chr17:69368 is **absent** from the current
  binaries: native = GATK = 70.51 (whole window) and native now follows GATK's window-context
  drift (67.43 in the single-locus window) instead of being constant at 74.01. Native
  full-precision TLOD at that locus is `70.514741250854556` on both backends.
* The documented `AS_SB_TABLE` mismatch at 69368 (`21,0|27,1`) is absent; current value
  `20,0|22,1` equals GATK.
* The two backends are bit-identical on all probed TLOD values.
* `FASTGATK_TLOD_FULL` is a log tag, not an env channel; the real gate is `FASTGATK_DEBUG_TLOD`.
* The delivered joint script was broken in two independent ways (reproduced above).

**Not proven / open:**

* **Machine-precision GATK agreement of TLOD** anywhere: GATK's VCF writes TLOD at 2 decimals,
  so only "≤ 0.005 at every probed locus" is established. No GATK-side full-precision channel
  was used.
* The recorded engine-parity headline **1.4e-14 does not reproduce** (current max
  **4.69e-13** at 17:69631, both backends). Whether this is a genuine shift in the engine
  numeric path since R24, or a difference in how R24's run was configured/fixtured, is
  **unverified** — this audit ran the shipped verifier only and did not diff engine internals.
* No causal attribution of the Δ3.5 disappearance to a specific commit: the working tree is
  clean and `git log` shows a single baseline commit for `calling_pipeline.cpp`, so history
  alone cannot justify any claim about "what fixed it".
* Regression coverage is limited to three fixtures; loci outside them are not covered, and the
  chr17 joint leg is a 0-row/0-row comparison (degenerate).
* Neither script is registered in CMake (snippets above are provided but deliberately not
  applied); they are therefore not part of any automated gate yet.
