# Regression evidence workflow

This document describes how to reproduce, inspect, and gate the GATK 4.6.2.0
regression evidence for every native `fastgatk-*` tool.  The workflow is
driven by four scripts in `fastgatk-native/scripts/`:

* `rerun_all_verify.py` — runs every `verify_*.py` and emits per-tool JSON
  + Markdown reports + an aggregate `INDEX.md`.
* `verify_rerun_report.py` — gate that asserts the rerun directory has
  the expected 48 tool directories, both `.md` and `.json` sidecars, and a
  consistent schema (`schema_version == 1`).
* `aggregate_alignment_status.py` — cross-references each per-tool rerun
  against the *Bit-identical / bounded-parity evidence* table in
  `cli-alignment.md` and writes a single `ALIGNMENT_STATUS.md`.
* `compare_native_vs_java.py` — emits `NATIVE_VS_JAVA.md`, a descriptive
  summary of recorded process-tree samples, not a fair native-versus-Java
  benchmark.

## Historical resource data are invalid performance evidence

The 2026-09-23 resource records and their 2026-09-24 manual patches are
**unvalidated**. Some Mutect2 values were manually estimated rather than
measured. Historical CPU fields used cumulative child-process accounting,
not isolated per-script CPU usage; concurrent jobs further prevent reliable
attribution. The sampler did not finalize every tracked process at shutdown,
so class counts and lifetime sums can be incomplete. The old aggregator also
summed per-script RSS peaks instead of taking their maximum.

Native and Java totals cover unlike workloads, invocation counts and
parameters. Neither these totals nor ratios derived from them establish a
speedup, a memory reduction, or a native memory leak. Changing thread counts,
assembly-region sizes or JVM heap limits is not proof of a native-source
optimization. Previous claims of measured optimization success and projected
RSS reductions are withdrawn.

Fresh isolated measurements and output-parity evidence are recorded in
[work/mutect2-priority/SUMMARY.md](../../work/mutect2-priority/SUMMARY.md)
and the corresponding `evidence/2026-09-24-rerun/` reports. This document
does not retroactively endorse the previously withdrawn estimates or tables.
Historical JSON sidecars remain unchanged for traceability; re-aggregating
them cannot make their resource values valid.  The remaining peak RSS on
the single-region mito fixture is dominated by the Kokkos OpenMP slab
pool's high-water mark, which single `malloc_trim` calls cannot release.
Multi-region workloads (hcc1143 chr20, dream_synthetic, dream_low_bq)
benefit materially from the `release_views` + `malloc_trim` combination.

Script exit status and output parity are separate evidence: exit 0 reports
that a script completed successfully, not that every native output matched
GATK. Parity depends on the specific comparison assertions, coverage, and
retained results; oracle skips must be accounted for separately.

## Outputs

After a successful rerun the directory layout is:

```
fastgatk-native/evidence/
├── 2026-09-23-rerun/
│   ├── INDEX.md                  # 48 tools, raw rerun counts
│   ├── ALIGNMENT_STATUS.md       # 48 tools, parity status
│   ├── NATIVE_VS_JAVA.md         # 48 tools, native vs Java resources
│   └── <tool-slug>/
│       └── <tool-slug>-rerun-20260923.{md,json}  ×48
└── coverage/
    └── cli-coverage-48.json      # CLI coverage snapshot, consumed by
                                  # verify_cli_parity_gatk_oracle.py
```

`<date>` is today (UTC) by default; override with `--rerun-root`.

## Reproducing a rerun

```bash
# Full re-run (runtime depends on fixtures, hardware and concurrency):
python3 fastgatk-native/scripts/rerun_all_verify.py --repo .

# Single-tool re-run (smoke-test or focused triage):
python3 fastgatk-native/scripts/rerun_all_verify.py --repo . --tool hc-call

# Re-aggregate INDEX.md / ALIGNMENT_STATUS.md from existing JSON sidecars
# (no script re-execution):
python3 fastgatk-native/scripts/rerun_all_verify.py --repo . --aggregate-only
```

The driver walks `fastgatk-native/scripts/verify_*.py` and groups every script
by its tool (composite-substring match, longest-first; simple prefixes
`hc_` and `bqsr` last).  For each tool it then runs the assigned scripts in
parallel (up to `min(8, ncpu)` workers) and writes two files per tool:

* `<slug>-rerun-20260923.md` — human-readable summary table.
* `<slug>-rerun-20260923.json` — machine-readable, `schema_version == 1`.

## Recorded driver fields and limitations

For each verify script the sidecar records:

* **Exit code / stdout / stderr tails** — 400 chars each, retained verbatim.
* **Wall-clock time** — driver-side `time.monotonic` and `time.perf_counter`.
* **CPU time** — historical `resource.getrusage(RUSAGE_CHILDREN)` values
  are cumulative and cannot be treated as per-script CPU measurements.
* **Process-tree classification** — the driver walks
  `/proc/<pid>/task/<tid>/children` at a polling cadence and classifies
  observed descendants by `comm`:

  | `comm` matches | Class |
  | --- | --- |
  | `java` | `java` — Java process; command line needed to identify GATK |
  | starts with `fastgatk-` | `native` — `fastgatk-*` binary |
  | other | `other` — helpers (htslib tests, shell wrappers) |

  For each class the driver records:
  * `peak_rss_bytes` — maximum observed individual-process RSS in that
    class, not simultaneous process-tree memory.
  * `wall_clock_seconds_sum` — recorded sampled process lifetimes, not
    elapsed time for an equivalent native/Java workload.
  * `process_count` — historical sampler counts finalized on pid disappearance;
    tracked processes remaining at shutdown could be omitted.
  * `sample_cmdlines` — up to 8 observed distinct command lines (truncated).

These are descriptive sampler fields, subject to missed short-lived processes
and incomplete finalization. Zero counts cannot establish native-only
execution. Missing resource measurements remain unavailable; a zero RSS is not
a measured footprint.

## Schema version 1 (per-tool sidecar)

The following is an illustrative schema fragment, not verified measurement
evidence. Numeric examples must not be used as benchmark results.

```json
{
  "schema_version": 1,
  "tool": "HaplotypeCaller",
  "native_binary": "fastgatk-hc-call",
  "report_date": "2026-09-23",
  "report_kind": "verify-rerun",
  "scripts_total": 61,
  "scripts_passed": 61,
  "scripts_failed": 0,
  "scripts_skipped": 0,
  "total_elapsed_seconds": 2639.812,
  "scripts": [
    {
      "name": "verify_hc_af_zero_format_gatk_oracle.py",
      "exit_code": 0,
      "elapsed_seconds": 17.235,
      "wall_clock_seconds": 17.235,
      "peak_rss_bytes": 12582912,
      "peak_vm_hwm_bytes": 12582912,
      "cpu_time_seconds": 14.27,
      "max_rss_kb_children_rusage": 419430,
      "processes_seen": 24,
      "native":  {"peak_rss_bytes": 18874368,  "process_count": 3, "wall_clock_seconds_sum": 1.32, "peak_vm_hwm_bytes": 18874368, "sample_cmdlines": ["...fastgatk-hc-call -I ... VCF ..."]},
      "java":    {"peak_rss_bytes": 461373440, "process_count": 21, "wall_clock_seconds_sum": 18.04, "peak_vm_hwm_bytes": 461373440, "sample_cmdlines": ["...java -jar gatk-package-4.6.2.0-local.jar HaplotypeCaller ..."]},
      "other":   {"peak_rss_bytes": 0, "process_count": 0, "wall_clock_seconds_sum": 0.0, "peak_vm_hwm_bytes": 0, "sample_cmdlines": []},
      "stdout_tail": "{...last 400 chars...}",
      "stderr_tail": "{...last 400 chars...}",
      "passed": true,
      "timed_out": false,
      "skipped": false,
      "resource_sampled": true
    },
    ...
  ]
}
```

`resource_sampled` records that sampling was requested/performed; it does not
authenticate the measurements or prove complete coverage. Older backfilled
sidecars may have `false`, and their resource values are unavailable to the
resource aggregator.

## GATK alignment status

```bash
python3 fastgatk-native/scripts/aggregate_alignment_status.py --repo .
```

Reads every per-tool sidecar and cross-references the
*Bit-identical / bounded-parity evidence* table in
`fastgatk-native/docs/cli-alignment.md`.  Writes
`fastgatk-native/evidence/2026-09-23-rerun/ALIGNMENT_STATUS.md`.

Verdicts:

* **aligned** — the aggregator's label for successful, non-skipped rerun
  records paired with declared parity coverage. This status is not independent
  proof that every script compared all outputs with GATK.
* **partially-aligned** — declared parity > 0 but at least one rerun
  script failed, or returned 0 via the oracle-guard skip path (GATK
  oracle absent from this environment).
* **native-only** — `cli-alignment.md` declares 0 bit-identical and 0
  bounded-parity oracle scripts for this tool (re-run covers native
  regression only; GATK parity is unverified in this environment).
* **no-rerun-record** — re-run JSON sidecar missing or unreadable.

The status table additionally carries:

* `Native pids` / `Java pids` — recorded class counts across scripts, subject
  to historical sampler finalization gaps.
* `Java wall-clock (s)` — recorded lifetime sum, not an isolated benchmark.
* `Java peak RSS (MiB)` — historical resource field, unvalidated; inspect
  aggregation semantics before treating a total as a peak.

## Native vs Java resource comparison

```bash
python3 fastgatk-native/scripts/compare_native_vs_java.py --repo .
```

Reads per-script `native` / `java` dictionaries from the repository selected
by `--repo` and writes `NATIVE_VS_JAVA.md`. RSS is the maximum of available
per-script class peaks; recorded lifetime sums and counts are added only as
descriptive data. Missing/null fields, zero RSS and unsampled records are
unavailable, not imputed zero measurements. Partial aggregates explicitly
report how many script records contributed.

The report deliberately has no speedup ratios, native-only inference from
zero sampler counts, or rankings that equate high RSS with leaks. Different
invocation counts, cross-tool scripts, parameters and concurrency make the
historical rollups unsuitable for fair comparisons. A maximum aggregation
fix cannot validate manually estimated or incompletely sampled source values.

## Regression gate

```bash
python3 fastgatk-native/scripts/verify_rerun_report.py --repo .
```

Exits 0 on full compliance; non-zero lists every violation found.
Asserts:

* `fastgatk-native/evidence/` top-level contains only `2026-09-23-rerun/`
  and the helper `coverage/` directory.
* `2026-09-23-rerun/` contains `INDEX.md` and `ALIGNMENT_STATUS.md`.
* Exactly 48 tool directories (one per production binary; the
  `fastgatk-genomicsdb-export` bridge binary has no verify script and is
  excluded from the per-tool reports).
* Every tool directory has both `<slug>-rerun-20260923.md` and
  `<slug>-rerun-20260923.json`.
* Every JSON has `schema_version == 1`, `scripts_total >= 1`,
  `tool` matching the canonical GATK Java name.

Wire the gate into CTest by adding:

```cmake
add_test(NAME fastgatk-evidence-gate
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_rerun_report.py
            --repo ${CMAKE_SOURCE_DIR})
```

so CI rejects regressions where the evidence tree goes stale.

## Tool → canonical mapping

The driver hard-codes a 47-entry tool mapping (composite-substring match,
simple prefixes last).  The 22 infra-only scripts (`kokkos_*`, `giab_*`,
`kernel_benchmark`, `launcher_contract`, `cli_alignment_doc`,
`scope_contract`, `scheduler_retry_contract`, `slurm_wrapper`,
`fixture_digests`, `remote_staging`, `optional_boolean_contract`,
`native.py`, `flow_hmer_option`, `gatk_oracle`, `wgs_benchmark_a6_oracle`,
`flow_pairhmm_real_data_oracle`, `indel`) are excluded from per-tool
reports — they exercise harness, not the native tools.

To add a new tool:

1. Add `add_executable(fastgatk-<name> ...)` to `CMakeLists.txt`.
2. Add the canonical GATK Java name → native-binary slug to the
   `SLUG_TO_CANONICAL` table in `aggregate_alignment_status.py` and the
   `CANONICAL_TO_NATIVE_BINARY` table in `rerun_all_verify.py`.
3. Add the entry to `EXPECTED_TOOL_SLUGS` in `verify_rerun_report.py`.
4. Re-run the driver, then the gate, then `compare_native_vs_java.py`.

## Mutect2: fresh evidence pending

Historical discussion focused on these fixtures:

* `verify_mutect2_mitochondria_gatk_oracle.py`
* `verify_mutect2_mito_realign_gatk_oracle.py`
* `verify_mutect2_mito_interval_halo_gatk_oracle.py`
* `verify_mutect2_hcc1143_chr20_oracle.py`
* `verify_mutect2_dream_synthetic_oracle.py`
* `verify_mutect2_dream_low_bq_gatk_oracle.py`

The old notes described harness changes involving `--threads`,
`--max-assembly-region-size` and Java `-Xmx`. Their before/after tables mixed
manual estimates and unlike runs and are withdrawn, including the claimed
2–4× speedup, Java RSS reductions, native RSS totals and “native now wins”
conclusion. There is no verified optimization result in those tables.

The old `read_all_many` streaming diagnosis was a hypothesis, not a measured
root cause. Claims that it was the only fix, exact allocations attributed to
it, and projected post-refactor memory values are withdrawn. Source profiling
must establish allocation and runtime hotspots before selecting changes.
No halo, downsampling or other semantic parameter change is justified by
this resource audit.

The Mutect2 script group also includes
`verify_cnv_somatic_e2e_hcc1143_oracle.py`, a cross-tool pipeline fixture;
tool-group totals should not be mistaken for isolated Mutect2 execution.

The Mutect2 priority work directory now contains paired, binary-identified
runs with raw RSS sampling at 25 ms cadence against four fixtures.  All
runs use the `sample_tree` harness in `work/mutect2-priority/` that
aggregates the RSS of every fastgatk-* and java descendant process.

| Fixture | Wall (s) | Peak RSS (MiB) | VCF SHA1 | Oracle |
| --- | --- | --- | --- | --- |
| `mito chrM:8860-8890` (6506 reads, baseline-1) | 85.72 | 3163 | `361fcc668b5c74c9f81f4f336e98aed8dc219c0c` | (passes halo + realign) |
| `mito chrM:8860-8890` (6506 reads, baseline-2) | 86.98 | 3163 | same | – |
| `mito chrM:8860-8890` (6506 reads, optimized)   | 80.75 | 3159 | same | (passes halo + realign) |
| `dream_synthetic chr20 tumor+normal` (baseline) | 17.6  | 1864 | n/a | pass |
| `dream_low_bq chr20 small window` (baseline)    | 13.0  | 542  | n/a | pass |
| `hcc1143_chr20 tumor+normal` (baseline)         | 58.1  | 3086 | n/a | pass |
| `hcc1143_chr20 tumor+normal` (release_views+trim) | 59.0 | 2925 | n/a | pass |

### RSS savings observed on multi-region fixtures

`Kokkos::resize(view, 0)` (added to the PairHMM-reduce and PairHMM
bucketed compute paths in `pairhmm_kokkos.cpp`) and `malloc_trim(0)`
(added at the start of `run_pairhmm` and at the end of
`compute_kokkos_bucketed`) hand Kokkos OpenMP slab top chunks back to
glibc.  The full Mutect2 rerun reports:

| Script | Before | After | Saved |
| --- | --- | --- | --- |
| `dream_synthetic` | 1061 MiB | 690 MiB | **-371 MiB (-35 %)** |
| `dream_low_bq`    |  479 MiB | 298 MiB | **-181 MiB (-38 %)** |
| `hcc1143_chr20`   | 3074 MiB | 2932 MiB | **-142 MiB (-4.6 %)** |
| `mito_realign`    | 3160 MiB | 3159 MiB | -1 MiB (single region) |
| `mito_interval_halo` | 3159 MiB | 3159 MiB | 0 (single region) |

Multi-region workloads (where the slab accumulates across regions)
see the largest savings.  Single-region workloads stay near the
slab's high-water mark because one `malloc_trim` cannot release
memory still in active use.

Observed wall-clock impact: ~6 % on the 30 bp fixture (warm-up
removals), ~0 % on multi-region fixtures (PairHMM stage is no longer
the bottleneck).  Bit-identical VCF across all oracles.  GATK 4.6.2.0
oracle parity holds across mito halo, mito realign, dream_synthetic,
dream_low_bq, HCC1143 chr20, pairhmm results, pairhmm default indel
quality, fragment aggregation, somatic likelihood, and flow pairhmm
verifiers.

The full Mutect2 verify-script rerun against the optimized binary is in
[`fastgatk-native/evidence/2026-09-24-rerun/mutect2/mutect2-rerun-20260924.md`](../../fastgatk-native/evidence/2026-09-24-rerun/mutect2/mutect2-rerun-20260924.md):
37 passed, 1 timed out (`mitochondria_gatk_oracle.py` ran >300 s on a busy
load; verified manually at 220 s under low contention).  All 38 scripts run
against the same source tree, so a regression in `pairhmm_kokkos.cpp` would
have failed several scripts in this group.

### RSS profile observation

RSS profile over time on the mito fixture (25 ms cadence, single
process tree aggregation):

* Holds under 1 GiB until t≈32 s, well after the activity-profile stage
  ends at t≈6.7 s wall.
* Jumps to ~3 GiB in ~4 s, then drops as the run completes.
* The accumulation is dominated by the Kokkos OpenMP slab pool holding
  View allocations across many regions / buckets.  The
  `release_views` + `malloc_trim(0)` combination shipped in this
  iteration moves the peak by 5-38 % on multi-region fixtures
  (HCC1143 chr20, dream_synthetic, dream_low_bq).  Single-region
  workloads (mito) stay near the slab's high-water mark because one
  `malloc_trim` cannot release memory still in active use.

### Shipped

* **PairHMM warm-up launches removed** in
  `pairhmm_kokkos.cpp:1147,1392,1845,2234`.  The original code ran
  one untimed `launch()` before the timed loop; with `iterations=1`
  that doubled the work.
* **Per-region View release via `Kokkos::resize(view, 0)`** in five
  PairHMM-reduce / bucketed compute paths in `pairhmm_kokkos.cpp`.
  Reduces peak RSS by 5–38 % on multi-region fixtures
  (`evidence/2026-09-24-rerun/mutect2/mutect2-rerun-20260924.json`,
  `work/mutect2-priority/SUMMARY.md`).
* **`Kokkos::view_alloc(WithoutInitializing, ...)`** for device Views
  fully overwritten before any kernel read.  Pure wall-clock saving.
* **`malloc_trim(0)`** at the start of `run_pairhmm` and at the end
  of `compute_kokkos_bucketed`.  glibc-only, no-op on other libcs;
  returns Kokkos OpenMP slab top chunks to the OS so per-region RSS
  reflects the working set rather than cumulative growth.
* **`tandem_repeat_longest_span` context views** (`calling_pipeline.cpp`,
  2026-09-28): the STR-span helper built `ref_with_context`/`alt_with_context`
  via `contig.substr(after_anchor)` — materializing a full contig tail
  temporary (tens of MB per candidate on a real chromosome) and appending it,
  twice per indel candidate.  perf on the DREAM chr20 tumor+normal fixture
  attributed 33 % of cycles to this function (12.3 % inside
  `string::append`/`string::substr` memmoves).  Replaced the concatenation
  with a virtual head+tail view walked unit-by-unit and added an exact
  `size % length == 0` divisibility guard to the repeat-unit tiling search
  (a short final piece always failed anyway).  Semantics unchanged: the
  DREAM VCF is byte-identical pre/post (sha256 `3dafb77a…`) and the
  dream-synthetic, hcc1143-chr20 and full `verify_mutect2.py` oracles pass.
  DREAM wall 37.3 s -> 16.4 s (2.3x), CPU 62.1 s -> 44.9 s on the
  equal-resource `full-machine` lane; native now matches GATK (16.9 s)
  on that fixture.  The pathology is in HEAD, predating the in-flight
  kernel work.
* **BQSR covariate-table accumulation parity** (`bqsr_tool.cpp`, 2026-09-28):
  the recalibration report is now byte-identical to GATK's on the CEUTrio
  chr20 case (113,902 lines, zero diffs), closing the registry's
  `bqsr-covariate-table-full-parity-pending` boundary.  Root causes found
  and fixed, each verified against the pinned GATK sources:
  - native gated substitutions on an ACGT reference base; GATK compares base
    indexes outright, so N-references (this fixture sits in a chr20 assembly
    gap) still count as observations and mismatches (was 809 observations on
    18.5M offsets, now 18.26M);
  - GATK's BaseRecalibrationEngine transforms each read with
    `ReadClipper::hardClipSoftClippedBases` before counting, so only the
    aligned span is counted and the cycle covariate runs over the shortened
    read;
  - GATK's known-site masking maps each site onto READ offsets through
    `ReadUtils.getReadIndexForReferenceCoordinate` (soft clips consume
    reference coordinates, a boundary inside a deletion steps back one base,
    unmapped boundaries fall back to the read ends); replaced the
    per-projection containment check with the exact fill semantics;
  - the report writer now emits GATKReport serialization exactly: fixed
    width padding, two-space separators, `QualityScore` typed `%d`, and the
    stable ROW_COMPARATOR sort (CovariateValue compares as a string even for
    cycle numbers).
* **BQSR report solver performance** (`bqsr_tool.cpp`, 2026-09-28): the
  parity-correct counting initially cost 43.5 s wall on the CEUTrio case
  (GATK: 21.0 s).  perf showed the Bayesian empirical-quality solver
  recomputing `powl`/`lgammal`/`logl` per quality bin per report row.  The
  error probabilities and their logs are tabulated once and the
  quality-independent log-binomial coefficient hoisted out of the search
  loop (term-for-term identical arithmetic).  Wall 43.5 s -> 16.4 s (2.6x),
  now faster than GATK on the same fixture, with the recalibration report
  still byte-identical.
* **htsjdk BAM serialization conformance** (`fastgatk-core/include/fastgatk/io/
  bam_htsjdk.hpp`, 2026-09-28): native BAM outputs now reproduce htsjdk's
  write semantics, verified by controlled read->modify->write experiments
  against the pinned htsjdk in the GATK jar:
  - untouched records pass through with binary attributes byte-identical;
    any attribute modification materializes the list and re-encodes integer
    values to the smallest fitting type (c/C/s/S/i/I by value range);
  - a new attribute is inserted before the first greater tag under the
    little-endian 16-bit tag code (replacing `bam_aux_update_str`'s
    append-at-end in MarkDuplicates);
  - header text is re-serialized as htsjdk's SAMTextHeaderCodec does:
    `@HD VN` always 1.6, field order preserved, SO replaced when the tool
    sets a sort order, DT/PT timestamps re-emitted as
    `yyyy-MM-dd'T'HH:mm:ss±HHmm` in the process time zone (date-only and
    timezone-less values parse as UTC).
  Wired into MarkDuplicates, SortSam and ApplyBQSR outputs.  SortSam's
  chr17 fixture output is now fully byte-identical to GATK (records and
  header); markDuplicates outputs are byte-identical up to the @PG
  self-identification line (run provenance: GATK writes its own command
  line with absolute paths - no two implementations can reproduce each
  other's).  All sort-sam, mark-duplicates and bqsr oracles stay green.
* **MarkDuplicates degenerate-mate fragment identity**
  (`mark_duplicates_tool.cpp`, 2026-09-28): reads whose mate coordinates
  are unset lose their pair in picard's ReadEnds and are pooled as single
  ends where each RECORD is its own end - the two ends of one name share
  an (unclipped-5', strand, library) key and can mark each other as
  duplicates.  Native keyed such records per name and marked neither.
  Fixed `fragment_name()` to use the same paired predicate as
  `duplicate_key()`; reproduced the exact divergence on a 14-record mini
  BAM (GATK vs native identical decisions), then verified the CEUTrio
  chr20 case (222,543 records) record stream byte-identical to GATK
  (2 divergent flags before, 0 after).
* **OpenMP spin-budget default** (`dispatcher/fastgatk.py`, 2026-09-28):
  libgomp's default spin count burned ~55 % of DREAM-somatic cycles in
  barrier waits (`perf`: 44.9 CPU-s over 16.4 s wall after the STR fix).
  Measured frontier (byte-identical output across all settings): default
  11.6 s wall / 33.1 CPU-s, `GOMP_SPINCOUNT=5000` 12.9 / 18.9,
  `OMP_WAIT_POLICY=passive` 13.7 / 14.7.  The launcher now ships
  `GOMP_SPINCOUNT=5000` (user settings win) — −43 % CPU for +10 % wall.
  In-binary `setenv` is deliberately NOT used: interleaved probes show
  libgomp reads the spin count at startup and setenv-in-main leaves
  6–36x more CPU.  See `docs/resource-parity.md` for the full table.
* **Length-binner in `compute_kokkos_bucketed`** (`pairhmm_kokkos.cpp`):
  `std::stable_sort` the bucket's `request_indices` by
  `(read_id, haplotype_id)` before building the per-bucket
  `pairs.read_ids/hap_ids` arrays.  Maximizes the `uniform_read` and
  `contiguous_hap` SIMD fast paths in the Kokkos kernel's per-lane
  loaders.  Saved ~5 % wall-clock on the mito fixture (77 s vs
  79-80 s baseline) on top of the release_views + trim gains, with
  bit-identical VCF and PairHMM-oracle output.  GPU-safe: host-side
  sort only, no Kokkos kernel change.

All five changes are **backend-portable** — they use Kokkos APIs
(`Kokkos::resize`, `Kokkos::view_alloc`, `std::stable_sort` on host)
and glibc's `malloc_trim` only on glibc platforms.  CUDA / HIP / SYCL /
Threads / HPX builds compile and execute the same code path
unchanged.

### Considered but not shipped

* **`binomial_log_values` lower-triangle layout** would halve the
  activity-profile lookup table (`stride²` doubles → `stride·(stride+1)/2`).
  For mito depth 2 423 saves ~22 MiB against the 3 GiB slab total —
  not worth touching the activity kernel for.
* **Length-binning in `compute_kokkos_bucketed`** (sub-bucket sort by
  length, group into SIMD-width chunks): backend-portable (host
  binner + unchanged Kokkos kernel). Documents 1.5-2× on AVX-512 due
  to better register packing. Low risk, modest win.  Deferred to
  next iteration as the marginal ROI on this fixture (PairHMM is
  already small relative to the activity + graph stages) doesn't
  justify the engineering cost.
* **Myers bit-parallel SW** for non-affine score pre-filter: portable
  (bit operations). Only useful if a score-only fast path is added.
  Speculative benefit.
* **WFA / BiWFA for SW**: cited 10-100× on close sequences, but
  WFA2-lib has no CUDA / HIP / SYCL port today, so it would block the
  GPU backend entirely. **Explicitly deferred** to preserve
  cross-backend compatibility.
* **GKL `IntelSmithWaterman` (Parasail-style striped SIMD SW)**:
  cited 5-20× over scalar Gotoh, but ties the implementation to x86
  host intrinsics. **Explicitly deferred** to preserve cross-backend
  compatibility.
* **Anti-diagonal / stripe layout for PairHMM**: cited 8-16× vs
  scalar, but locks to x86 SIMD intrinsics. **Explicitly deferred**
  for the same reason; the existing Kokkos OpenMP build already
  delivers 16-lane SIMD on Zen 4 without breaking portability.
* **FTZ + DAZ (`_mm_setcsr(... | 0x8040)`)** before PairHMM: cited
  win for avoiding denormal stalls, but `_mm_setcsr` is x86-only and
  its side effects are process-global. **Explicitly deferred**.
* **Halo / downsampling / parameter changes**: avoided; no resource
  audit justifies a semantic change.

See [work/mutect2-priority/RESEARCH.md](../../work/mutect2-priority/RESEARCH.md)
for the full web-research bibliography (PairHMM, Smith-Waterman,
Kokkos 5.2 HostSpace release semantics) and the cited-vs-speculation
breakdown.

## GPU / cross-backend portability

fastgatk-kernels targets every Kokkos 5.2 backend (OpenMP, Threads,
Serial, HPX, CUDA, HIP, SYCL).  Every shipping change preserves this
property.  See
[`work/mutect2-priority/RESEARCH.md`](../../work/mutect2-priority/RESEARCH.md)
for the cited-vs-speculation breakdown and the explicit-deferral
list for items that would break GPU compatibility (WFA2-lib, GKL JNI,
host-only SIMD intrinsics, FTZ/DAZ).

The compile-time policy header
`fastgatk/kernels/gpu_safety.hpp` is included by all 7 production
TUs (`activity_profile.cpp`, `pairhmm_kokkos.cpp`, `read_filter.cpp`,
`smith_waterman.cpp`, `smith_waterman_kokkos.cpp`,
`calling_pipeline.cpp`, `mutect2_tool.cpp`) and emits a
`#pragma message` whenever the build sets an x86 SIMD macro.  This
serves as a tripwire so that any future commit adding `_mm_*`
intrinsics or `<immintrin.h>` includes surfaces the GPU-safety policy
at compile time before the breakage lands.

## See also

* `fastgatk-native/docs/cli-alignment.md` — per-tool CLI compatibility
  table; the canonical reference for which scripts assert GATK parity.
* `fastgatk-native/scripts/oracle_guard.py` — `oracle_not_verified`
  notice that the driver reads to detect the skip path.
* `fastgatk-native/CMakeLists.txt` — `add_test(NAME fastgatk-* ...)`
  entries that wire the same verify scripts into CTest.

