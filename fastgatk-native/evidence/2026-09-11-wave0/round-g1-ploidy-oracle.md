# Round G1 — ploidy-3 window-invariance oracle for the D2 ploidy-path fix

Deliverable: `fastgatk-native/scripts/verify_hc_ploidy_window_invariance_gatk_oracle.py`
Report: `.diag/round-g1-ploidy-oracle.md` (this file). No source, header, kernel, CMake or
`*.md`-at-root file was created, edited, rebuilt or committed in this round.

## 1. What the gate is

The diploid twin (`fastgatk-native/scripts/verify_hc_window_invariance_gatk_oracle.py`) gates the
diploid gVCF annotation path. The arbitrary-ploidy path (`sample_ploidy != 2`) has its **own**
annotation owner selection (`annotations_from_owner`, `hc_call.cpp` line 4553 in the
`sample_ploidy != 2` branch; the diploid twin's copy is at line 4886), and that second
copy carried a live instance of D2 until commit `a942fb1` (2026-09-11 03:24:49). This script is the
missing regression guard for that fix; it is written to mirror the diploid twin so the two gates
behave identically (`--native`/`FASTGATK_HC_BINARY`, `--serial-native`/`FASTGATK_HC_SERIAL_BINARY`,
`--expect-divergence`, `--gate-structural-windows`, `--windows`, `--threads`,
`tempfile.TemporaryDirectory` for all scratch, `main() -> int`, final JSON status line, exit code 1
on any violation, exit 0 with `status: skipped` when oracle assets are absent unless
`FASTGATK_REQUIRE_GATK_ORACLE` is set).

Per window in the pinned set it runs, with **`--sample-ploidy 3` on both sides**:

```
third_party/jdk17/bin/java -Xmx1g -jar third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
    HaplotypeCaller -R fixtures/chr20/ref20mnp.fasta -I fixtures/chr20/mnp.bam -L 20:<start>-10020710 \
    --emit-ref-confidence GVCF --max-mnp-distance 1 --sample-ploidy 3 --native-pair-hmm-threads 1 \
    --create-output-variant-index false --add-output-vcf-command-line false -O <tmp>/gatk.<start>.vcf

fastgatk-native/build/fastgatk-hc-call  -R ... -I ... -L 20:<start>-10020710 \
    --emit-ref-confidence GVCF --max-mnp-distance 1 --sample-ploidy 3 --threads 1 \
    --add-output-vcf-command-line false -O <tmp>/native.<start>.vcf
```

and then asserts:

1. **data-row byte parity** — the native data lines must equal **that window's live GATK data
   lines** (header/`#` provenance lines ignored); and
2. **the D2 row** — `POS 10020680` (`ALT=AT,*,<NON_REF>`) `RAW_MQandDP` + `SB` must be
   byte-identical to GATK in **every** window in the set, including the two only partially gated
   ones.

The expected GATK value is taken from the **live GATK run of that same invocation**, never from a
constant. The pinned table `GATK_D2_EXPECTATION` is a *drift guard on the GATK side only*: if live
GATK ever stops producing the recorded per-window value, the run fails with a "GATK baseline
drifted" violation instead of silently agreeing with a broken native build.

## 2. Windows, and why two of them are not fully gated

| start | class | full-row byte parity asserted | D2 row (10020680) asserted |
|---|---|---|---|
| 10019901 | gating | yes | yes |
| 10020381 | gating | yes | yes |
| 10020391 | gating | yes | yes |
| 10020401 | gating | yes | yes |
| 10020411 | gating | yes | yes |
| 10020421 | structural | **no** by default (`--gate-structural-windows` opts in) | yes |
| 10020431 | structural | **no** by default (`--gate-structural-windows` opts in) | yes |

Exclusion reason for `10020421` / `10020431`, observed this round: native does not reproduce GATK's
`<NON_REF>` **reference-block granularity** there — it collapses many GATK blocks into a single
degenerate block (`GQ=0`, `PL=0,0,0`, and a merged `END`), e.g. at `10020430` GATK
`GT:DP:GQ:MIN_DP:PL = 0/0/0:12:21:12:0,21,57,498` vs native `0/0/0:31:0:31:0,0,0,0` (98 GATK rows
vs 12 native rows in that window). That is the **D3 family, a different defect** (the diploid twin
already documents and excludes the same two starts for the same reason) and it is not what this gate
owns. Excluding them from *full-row* parity is what keeps the gate honest rather than permanently
red for an unrelated defect; their D2 row — the subject of this oracle — is still asserted in every
run, and their row counts and values are still printed. Requirement 2's byte-identical assertion is
therefore enforced in full for the five windows where it is not confounded by D3.

## 3. Per-window observed values (live, this round)

Both backends (OpenMP `fastgatk-native/build/fastgatk-hc-call` and Serial
`fastgatk-native/build-serial/fastgatk-hc-call`) produced identical results.

| window (`-L 20:<start>-10020710`) | rows GATK/native | full-row parity | GATK @10020680 `RAW_MQandDP` / `SB` | native @10020680 `RAW_MQandDP` / `SB` |
|---|---|---|---|---|
| 10019901 | 94 / 94 | **identical** | `28800,8` / `0,0,3,3` | `28800,8` / `0,0,3,3` |
| 10020381 | 13 / 13 | **identical** | `28800,8` / `0,0,3,3` | `28800,8` / `0,0,3,3` |
| 10020391 | 13 / 13 | **identical** | `28800,8` / `0,0,3,3` | `28800,8` / `0,0,3,3` |
| 10020401 | 13 / 13 | **identical** | `28800,8` / `0,0,3,3` | `28800,8` / `0,0,3,3` |
| 10020411 | 13 / 13 | **identical** | `28800,8` / `0,0,3,3` | `28800,8` / `0,0,3,3` |
| 10020421 | 98 / 12 | differs (D3 blocks) | `97200,27` / `1,2,6,16` | `97200,27` / `1,2,6,16` |
| 10020431 | 94 / 9 | differs (D3 blocks) | `97200,27` / `1,2,6,16` | `97200,27` / `1,2,6,16` |

`ALT` at `10020680` is `AT,*,<NON_REF>` in all seven windows on both sides.

Note on requirement 3: at **ploidy 3** GATK is *also* not globally window-invariant at this locus —
it emits `28800,8` / `0,0,3,3` for the first five starts but `97200,27` / `1,2,6,16` for the starts
`10020421` and `10020431`. Both of those legitimate GATK values are pinned *per window* (and
re-derived live on every run), not assumed away. (**Speculation**, not observed: the natural
explanation is that a later window start changes the retained read set / assembly-region boundary
around the CA→AT MNP and therefore the informative-overlap population GATK annotates from; I did
not instrument GATK to confirm this mechanism.)

## 4. Exit status in both modes

| command | observable | result |
|---|---|---|
| `python3 …/verify_hc_ploidy_window_invariance_gatk_oracle.py --serial-native fastgatk-native/build-serial/fastgatk-hc-call` | stdout JSON `"status": "pass"`, `"strict_mode": true`, `"violations": []` | **exit 0** (ran twice: 32.4 s and 33.8 s wall, both backends) |
| same `+ --expect-divergence` | stdout JSON `"status": "diagnostic"`, `"strict_mode": false` | **exit 0** |
| `--native /nonexistent` (assets missing, no `FASTGATK_REQUIRE_GATK_ORACLE`) | `{"status": "skipped", …}` | exit 0 |
| `--native /nonexistent` with `FASTGATK_REQUIRE_GATK_ORACLE=1` | `missing oracle assets: ['/nonexistent']` | **exit 1** |

Default (strict) runtime ≈ **33 s** for 7 windows × (GATK + native) with both native backends,
inside the ~3 min budget; the JSON status line is the last line of stdout and parses.

## 5. Negative case: was it demonstrated?

**Plainly stated: I could NOT demonstrate it against a real pre-fix binary**, and I am not claiming
otherwise.

* No pre-fix native binary exists to run: the only two `fastgatk-hc-call` artifacts
  (`fastgatk-native/build/`, `fastgatk-native/build-serial/`, mtimes 2026-09-11 03:07:12 and
  03:07:51) were built *after* the last source edit (`src/hc_call.cpp`, mtime 03:06:50) and
  therefore already contain the fix; `build/` is gitignored (`.gitignore:27`) and no binary was ever
  committed, so `git log --all --diff-filter=A -- '*fastgatk-hc-call'` is empty.
* Rebuilding a pre-fix binary or reverting the source is forbidden by this round's constraints
  (no source edits, no rebuild, no commit), so the true negative case was out of reach.

What I *did* do instead, and it is a simulation, not a pre-fix run: a throwaway harness (created
under `tempfile.TemporaryDirectory`, deleted with it, no repo file touched) loaded the oracle as a
module and monkeypatched its `run_native` so that, after running the real post-fix binary, it
rewrote the `POS 10020680` row of the native VCF with the **documented pre-fix values recorded by
the fix itself** in commit `a942fb1` (`RAW_MQandDP=97200,27`, `SB=0,0,0,0` for starts
`10020381/10020391/10020401/10020411`; window `10019901` left untouched because it was already
correct pre-fix). Results:

* strict mode over those five windows → **exit 1**, 8 violations, each D2-row violation naming the
  signature: `POS 10020680 annotation differs (GATK RAW_MQandDP=28800,8 SB=0,0,3,3 vs native
  RAW_MQandDP=97200,27 SB=0,0,0,0) -- matches the documented PRE-FIX ploidy-3 signature`; the
  injected window `10019901` stayed clean, reproducing the pre-fix window pattern exactly;
* the same injected run with `--expect-divergence` → **exit 0**, `"status": "diagnostic"`, 8
  violations reported.

So the gate provably **bites on the exact pre-fix value set** and its strict/diagnostic exit codes
invert correctly; what remains unproven is only that the *pre-fix binary itself* produced those
bytes on this machine — that claim rests on the fix commit's own recorded pre-fix observation
(`a942fb1` message), which also states the fix was in `src/hc_call.cpp` only.

## 6. CMake registration snippet — for a later round, NOT applied

Insert immediately after the existing diploid block (`fastgatk-native/CMakeLists.txt` lines
1414-1424, `add_test(NAME fastgatk-hc-window-invariance-gatk-oracle …)`). Do not add
`--serial-native` here: the diploid twin's CTest entry does not, and the two gates should behave
the same way (the serial backend was still exercised manually this round — see §3).

```cmake
    # STRICT gate for the second live instance of D2 (the arbitrary-ploidy path):
    # at --sample-ploidy 3 the published allele at 20:10020680 must be annotated
    # from the AssemblyRegion owner that carries its PairHMM context, so
    # RAW_MQandDP/SB must equal pinned GATK's live per-window value in all seven
    # windows (two of which are only checked on that row because the separate
    # reference-block-granularity defect D3 owns their full row set).  The script
    # exits non-zero on any mismatch and is green as of commit a942fb1
    # (7 windows x GATK + native, ~35 s).
    add_test(NAME fastgatk-hc-ploidy-window-invariance-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_ploidy_window_invariance_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-ploidy-window-invariance-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
```

## 7. Residual risk / what this guard does not cover

* It guards **one** ploidy value (3, the highest GATK accepts here). The commit message for
  `a942fb1` notes other `assembly_region_likelihood_results` traversals in `hc_call.cpp`
  (3134/3227/3240/3521/4021/4116/4535/4664/4705) whose first-success owner selection was not
  individually adjudicated; only ploidy 3 exercises the arbitrary-ploidy annotation branch, so a
  third instance in a *different* branch would still escape both oracles. (**Speculation:** whether
  such a third instance exists is untested here; this round did not survey those lines.)
* `10020421` / `10020431` full-row divergence is reported, not asserted, so D3 regressions in those
  windows are visible in the log but do not fail this gate (`--gate-structural-windows` opts in at
  the cost of a red test until D3 lands).
