# Round: SILENT DATA LOSS on a malformed FORMAT record

Scope: the side finding of `.diag/round-crosssample-merge.md` §0 — on a record
whose FORMAT declares 4 keys but carries 5 sample values, pinned GATK exits 3
while native exits 0 and silently drops the record.  Gate first, fix second.

Probes: `.diag/malformed_probe.py` (`.diag/malformed-probe.log`, 18 shapes),
`.diag/malformed_reach_probe.py` (`.diag/malformed-reach.log` pre-fix /
`.diag/malformed-reach-after.log` post-fix).
New gate: `fastgatk-native/scripts/verify_genotype_gvcf_malformed_gatk_oracle.py`
(`.diag/malformed-gate-before.log`, `.diag/malformed-gate-after.log`).

## 0. Verdict

**(i) A genuine defect — and worse than the previous round reported.**

The drop is not record-local.  Native's linear traversal read its input with
`while (bcf_read(...) == 0)`, which cannot tell htslib's *end of input* (-1)
from htslib's *record did not parse* (-2): the first malformed record ended the
traversal as if the file had ended, so **that record and every record after it
in that input were discarded**, with exit 0 and no message beyond htslib's own
`[E::vcf_parse_format_max3]` log line.  Measured pre-fix: a malformed record at
chr1:2 followed by a well-formed record at chr1:5 produced **0 output rows,
exit 0**; with the malformed record *last*, the chr1:2 row was written and
everything from the malformed record on was lost, still exit 0.

Native is not entitled to handle this differently, because native *already has*
the correct behaviour in the same file: the indexed traversal
(`genotype_gvcf_tool.cpp:6505`) and the indexed streaming traversal (`:4943`)
check `status < -1` and throw.  Measured side by side on one input
(`.diag/malformed-reach*.log`):

| traversal | pre-fix | post-fix |
| --- | --- | --- |
| default, unindexed (`-V x.g.vcf`) | **exit 0, 0 rows** | exit 2, `BAD_INPUT: GenotypeGVCFs input record parse failed` |
| `--stream-by-locus`, unindexed | **exit 0, 0 rows** | exit 2, `BAD_INPUT: streaming GenotypeGVCFs VCF record parse failed` |
| `-V x.g.vcf.gz -L chr1:1-100` (tabix) | exit 2, reported | exit 2, reported (unchanged) |
| `-L chr1:1-100 --stream-by-locus` | exit 2, reported | exit 2, reported (unchanged) |

So the same input and the same tool already said "malformed input" on the
indexed path; only the unindexed path — the *default* one — swallowed it.

Reachability: the malformed record itself is hand-made, but it is ordinary
input-side corruption (a truncated/interrupted VCF line, a record written by a
third-party tool with a FORMAT/value count mismatch).  The severity is in the
*consequence*: any such record silently truncates the run's output while the
process reports success, so a downstream pipeline cannot detect it.  No warning
is emitted by native; the only hint is htslib's stderr line, which native's own
exit status contradicts.

## 1. STEP 1 — measured shapes (pre-fix)

Fixture per shape: the malformed record at chr1:2, a well-formed
`GT:DP:AD:PL`/`0/1` record at chr1:5; reference = 100 bp of `A`; both tools get
identical `-R/-V`, native additionally `--gatk-compatible-annotations`.
`rows` counts data rows written.  Full log: `.diag/malformed-probe.log`.

| shape | GATK exit | GATK rows | native exit | native rows | native stderr (htslib) |
| --- | --- | --- | --- | --- | --- |
| `too-many-format-values` (4 keys, 5 values) | 3 | 0 | **0** | **0** | `[E::vcf_parse_format_max3] Incorrect number of FORMAT fields at chr1:2` |
| `too-many-format-values-alone` | 3 | 0 | **0** | 0 | same |
| `empty-sample-column` | 0 | 1 | **0** | **0** | `[E::vcf_parse_format_check7] Number of columns at chr1:2 does not match the number of samples (0 vs 1)` |
| `dp-non-numeric` | 3 | 0 | **0** | **0** | `[W::vcf_parse_format_fill5]` + `[E::vcf_parse_format_fill5] Invalid character 'a' in 'DP'` |
| `ad-non-numeric` | 0 | 2 | **0** | **0** | `[W::vcf_parse_format_fill5]` + `[E::vcf_parse_format_fill5] Invalid character 'x' in 'AD'` |
| `pl-non-numeric` | 0 | 1 | **0** | **0** | `[W::vcf_parse_format_fill5]` + `[E::vcf_parse_format_fill5] Invalid character 'x' in 'PL'` |
| `gt-non-numeric` | 3 | 0 | **0** | **0** | `[E::vcf_parse_format_fill5] Couldn't read GT data: value not a number or '.'` |
| `missing-sample-column` (no sample column) | 3 | 0 | 0 | 2 | `[E::vcf_parse_format_empty1] FORMAT column with no sample columns` |
| `pl-width-short` (2 of 6 PLs) | 3 | 0 | 0 | 2 | (none) |
| `gt-allele-index-out-of-range` (`0/5`) | 3 | 0 | 0 | 2 | (none) |
| `extra-format-key-with-value` (undeclared `ZZ`) | 3 | 0 | 0 | 2 | `[W::vcf_parse_format_dict2] FORMAT 'ZZ' ... not defined in the header` |
| `too-few-format-values` (4 keys, 3 values) | 0 | 1 | 0 | 1 | (none) |
| `pl-width-long` (8 PLs) | 0 | 2 | 0 | 2 | (none) |
| `pl-width-tetraploid` (10 PLs, diploid GT) | 0 | 2 | 0 | 2 | (none) |
| `pl-empty` (`PL=.`) | 0 | 1 | 0 | 1 | (none) |
| `gt-single-allele-tetraploid-format` (`0/1/0/1`) | 0 | 2 | 0 | 2 | (none) |
| `ad-width-short` (2 of 3 ADs) | 0 | 2 | 0 | 2 | (none) |
| `format-key-not-in-header` (`ZZ` key, no value) | 0 | 2 | 0 | 2 | `[W::vcf_parse_format_dict2]` |

Answers to the STEP 1 questions:

* **Which drop silently**: the seven shapes in **bold** above (the first/third
  through seventh rows) exit 0 having written no row at all — the malformed
  record *and the well-formed record behind it*.  All seven are htslib
  `vcf_parse` failures (return -2).
* **How far the drop reaches**: to the end of that input.  Measured: 0 of 2
  records emitted with the bad record first; 1 of 2 emitted with the bad record
  last (`.diag/malformed-gate-before.log`, case
  `too-many-format-values-bad-last`: *"native wrote 1 data row(s) for a failed
  run"*).  Records *before* the malformed one survive; the malformed one and
  everything after are lost.
* **Realistic reachability**: not from a well-formed gVCF.  It is reachable from
  a corrupted or partially written input (a truncated line at any position), and
  from any producer that emits a FORMAT/values mismatch — i.e. it is the
  standard failure mode a pipeline must notice, not an exotic hand-crafted VCF.
* **Not fixed here** (measured, reported, no data loss): `missing-sample-column`
  (native emits a record with no sample data where GATK exits 3),
  `pl-width-short` (native writes a padded `PL=100,0,.`), and
  `gt-allele-index-out-of-range` / `extra-format-key-with-value` (native exit 0
  where GATK exits 3).  See §6.

## 2. GATK's rule, with file:line

GATK's reader rejects the record before any tool code sees it; the exception is
raised in htsjdk:

```
htsjdk.tribble.TribbleException: The provided VCF file is malformed at approximately line
number 11: there are 1 genotypes while the header requires that 9 genotypes be present for
all records at chr1:5, for input source: .../x.g.vcf
  at htsjdk.variant.vcf.AbstractVCFCodec.generateException(AbstractVCFCodec.java:887)
  at htsjdk.variant.vcf.AbstractVCFCodec.createGenotypeMap(AbstractVCFCodec.java:759)
  at htsjdk.variant.vcf.AbstractVCFCodec$LazyVCFGenotypesParser.parse(AbstractVCFCodec.java:121)
  at htsjdk.variant.vcf.AbstractVCFCodec.parseVCFLine(AbstractVCFCodec.java:454)
  at htsjdk.tribble.TribbleIndexedFeatureReader$WFIterator.readNextRecord(TribbleIndexedFeatureReader.java:377)
  at org.broadinstitute.hellbender.engine.VariantLocusWalker.traverse(VariantLocusWalker.java:131)
```

* The genotypes/keys count check is `AbstractVCFCodec.createGenotypeMap`
  (htsjdk 4.x, bundled in `gatk-package-4.6.2.0-local.jar`, class verified
  present with `javap`); the rejection is thrown as `TribbleException` and
  propagates out of the record iterator, so GATK exits **3** and writes no
  further row (rows already written stay written).
* `gatk-source/` contains **no** PL/format-width validation on this path: the
  only `getPL().length` comparison in the whole GATK tree is
  `gatk-source/src/main/java/org/broadinstitute/hellbender/utils/variant/writers/GVCFBlockCombiner.java:145`
  (block coalescing, unrelated), and `GenotypeGVCFsEngine.java` touches PL only
  at `:377-378`, `:490-491` (assignment method, `builder.noPL()`).  The rule is
  a *reader* rule: GATK rejects malformed input instead of genotyping it.
* htslib, native's reader, applies the same rule and *says so*:
  `third_party/htslib-build/htslib-src/vcf.c` — `vcf_parse_format_max3` (:2970),
  `vcf_parse_format_fill5` (:3135), `vcf_parse_format_check7` (:3424),
  `vcf_parse_format_empty1` (:2897); every one of them logs
  `[E::vcf_parse_format_*]` and makes `vcf_parse` (:3747, `int ret = -2` at
  :3749) return **-2**.  Native's error is therefore not "GATK is stricter";
  native had the verdict in hand and dropped it on the floor.

## 3. Native's skip path, with file:line

htslib's contract, read from the vendored source: `bcf_read` (:2019) returns
`0` = record, `-1` = **end of input** (`vcf_read` :3930 → `hts_getline` EOF),
`-2` = **parse failure** (`vcf_parse1` is a macro for `vcf_parse`,
`htslib/vcf.h:274`; `vcf_parse` returns -2 from its `err:` label).

Pre-fix (`git show HEAD:fastgatk-native/src/genotype_gvcf_tool.cpp`):

| line | code | effect |
| --- | --- | --- |
| `:6515` | `while (bcf_read(input, header, record) == 0) process_record(record);` | **the default traversal**: -1 and -2 both end the loop |
| `:5595` | `while (bcf_read(input, header, record) == 0)` (streaming span probe) | same |
| `:4923` | `return bcf_read(cursor.input, cursor.header, cursor.raw) == 0;` (streaming cursor, unindexed branch) | same |
| `:6505` | `if (status < -1) { ... throw "BAD_INPUT: indexed GenotypeGVCFs traversal failed"; }` | **already correct** |
| `:4943` | `if (status < -1) throw "BAD_INPUT: indexed streaming GenotypeGVCFs traversal failed";` | **already correct** |
| `:6490-6495` | `if (vcf_parse(&line, header, record) < 0) throw "BAD_INPUT: indexed GenotypeGVCFs VCF record parse failed";` | **already correct** |

That is the whole defect: three unindexed read loops did not apply the check the
indexed loops in the same file already applied.  The tools' error convention is
`throw std::runtime_error("BAD_INPUT: ...")` → `fastgatk-genotype-gvcf: <what>`
on stderr, exit **2** (`genotype_gvcf_tool.cpp:7097-7098`).

## 4. STEP 2 — the gate

New file: `fastgatk-native/scripts/verify_genotype_gvcf_malformed_gatk_oracle.py`
(exit status, not rows, so it is separate from the data-row oracles, as the task
suggested).  Contract asserted per case:

> on a malformed record native must either (a) exit non-zero, print a diagnostic
> and write **no** data row, or (b) exit 0 and still emit every well-formed
> record.  Exit 0 with a well-formed record missing is a violation.

14 cases: 8 gated fail-loud (the seven silent-drop shapes plus
`too-many-format-values-bad-last`, which also pins that no partial file is
written), 1 gated control (`control-valid-input`, which fails if the fix
over-rejects valid input), and 5 REPORTED-ONLY cases with `expect` constants
that fail loudly if pinned GATK moves.  Each case also pins GATK's measured exit
and rows, so the gate cannot silently re-pin itself to a changed GATK.

Exit status before the fix: **1** — `"status": "divergence"`, 25 violations,
all of them the silent-drop contract
(`.diag/malformed-gate-before.log`).  After the fix: **0** —
`"status": "pass"`, `"violations": []` (`.diag/malformed-gate-after.log`).

**Registration is NOT done** (hard constraint 1 forbids editing
`fastgatk-native/CMakeLists.txt`).  The orchestrator should add, next to the
spandel oracle at `CMakeLists.txt:1556-1560`:

```cmake
    # GenotypeGVCFs must not silently truncate its input on a malformed record
    # (bcf_read returns -1 for EOF and -2 for a parse failure).
    add_test(NAME fastgatk-genotype-gvcf-malformed-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_malformed_gatk_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-malformed-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
```

## 5. STEP 3 — the fix (applied, minimal)

One file, +34/-4, three call sites — no engine restructure, no Kokkos ABI
change, no other tool touched:

```diff
--- a/fastgatk-native/src/genotype_gvcf_tool.cpp
+++ b/fastgatk-native/src/genotype_gvcf_tool.cpp
@@ -4919,8 +4919,17 @@ bool read_next_genotype_cursor_record(...) {
-    if (cursor.traversal_index == nullptr)
-        return bcf_read(cursor.input, cursor.header, cursor.raw) == 0;
+    if (cursor.traversal_index == nullptr) {
+        // htslib signals end of input with -1 and an unparseable record with
+        // -2 (vcf_parse's error return).  Collapsing the two into "no more
+        // records" would silently truncate the traversal at the first
+        // malformed record, so the parse failure is reported instead.
+        const int status = bcf_read(cursor.input, cursor.header, cursor.raw);
+        if (status < -1)
+            throw std::runtime_error(
+                "BAD_INPUT: streaming GenotypeGVCFs VCF record parse failed");
+        return status == 0;
+    }
@@ -5592,7 +5601,8 @@ (streaming span probe)
-            while (bcf_read(input, header, record) == 0) {
+            int probe_status = 0;
+            while ((probe_status = bcf_read(input, header, record)) == 0) {
@@ -5620,6 +5630,15 @@
             }
+            // A record htslib could not parse (-2) must not read as end of
+            // input (-1); the probe would otherwise index a truncated file.
+            if (probe_status < -1) {
+                bcf_destroy(record);
+                bcf_hdr_destroy(header);
+                bcf_close(input);
+                throw std::runtime_error(
+                    "BAD_INPUT: streaming GenotypeGVCFs VCF record parse failed");
+            }
             bcf_destroy(record);
@@ -6512,8 +6531,19 @@ (default traversal)
             if (!used_index) {
-                while (bcf_read(input, header, record) == 0)
+                int read_status = 0;
+                while ((read_status = bcf_read(input, header, record)) == 0)
                     process_record(record);
+                // -1 is end of input; -2 is a record htslib could not parse.
+                // Treating the parse failure as end of input would silently
+                // discard this record and everything after it.
+                if (read_status < -1) {
+                    bcf_destroy(record);
+                    bcf_hdr_destroy(header);
+                    bcf_close(input);
+                    throw std::runtime_error(
+                        "BAD_INPUT: GenotypeGVCFs input record parse failed");
+                }
             }
```

Post-fix behaviour (`.diag/malformed-gate-after.log`): the eight gated shapes
exit **2** with `fastgatk-genotype-gvcf: BAD_INPUT: ... record parse failed` and
`native_output_exists=False` (no partial file, no partial rows), while
`control-valid-input` still matches pinned GATK byte for byte.  Three of the eight
(`empty-sample-column`, `ad-non-numeric`, `pl-non-numeric`) now exit non-zero
where GATK exits 0 and merely skips the site: stricter, but exactly the
project's rule that an unreadable record must not be discarded silently.

Not fixed, deliberately (see §6): the four shapes where native exits 0 with
*rows present*, i.e. where data is not lost.

## 6. STEP 4 — stale assertions

**None; no test line was corrected.**  Two independent checks:

1. Structural scan: all 147 VCF-ish data rows built inside
   `fastgatk-native/scripts/*.py` were checked for a FORMAT-key/sample-value
   count mismatch.  Six apparent hits, all verified by reading to be multi-line
   Python string continuations (`verify_genotype_gvcf.py:634,692`,
   `verify_gatk_genotype_gvcf_inbreeding.py:117`,
   `verify_select_variants_refonly_gatk_oracle.py:174,307`,
   `benchmark_variant_filtration.py:39`) — no fixture is malformed this way.
2. Empirical: the full 300-test suite passes on both backends after the change
   (§7), so no registered test depended on the truncation.

Adjacent divergences found but **not** repaired this round (recorded as
REPORTED-ONLY gate cases so they cannot drift unnoticed; each needs its own
decision, none of them loses data):

* `missing-sample-column`: htslib logs `[E::vcf_parse_format_empty1]` but its
  caller turns that -1 into success (`vcf_parse_format`, `vcf.c:3446`,
  `:3457-3458`:
  `if ((ret = vcf_parse_format_empty1(...))) return ret ? 0 : -1;`), so the
  record parses with zero samples and native writes
  `chr1 2 . A G . . DP=20` — a row with no FORMAT and no genotype — where GATK
  exits 3.  Follow-up: reject when
  `bcf_hdr_nsamples(header) > 0 && record->n_sample != bcf_hdr_nsamples(header)`
  (htslib sets `BCF_ERR_NCOLS` on that record).
* `pl-width-short`: native pads the PL vector (`PL=100,0,.`) where GATK exits 3.
* `gt-allele-index-out-of-range`, `extra-format-key-with-value`: native exit 0
  vs GATK exit 3.
* Systemic, out of scope (constraint 3): the same unchecked
  `while (bcf_read(...) == 0)` idiom still exists in 19 other
  `fastgatk-native/src/*.cpp` files (e.g. `apply_vqsr_tool.cpp:424,779`,
  `mutect2_tool.cpp:319`, `variant_recalibrator_tool.cpp:651,2395,2707`).
  Some tools already apply the `status < -1` check (`gather_vcfs_tool.cpp:670`,
  `genomicsdb_import_tool.cpp:864,1211`, `variant_eval_tool.cpp:901`,
  `combine_gvcf_tool.cpp:1178`, `get_pileup_summaries_tool.cpp:833`,
  `collect_f1r2_counts_tool.cpp:489`, `read_metrics_tool.cpp:1280,1289`), so the
  fix above follows the house convention — but the idiom should be audited
  tool by tool, since each site is an independent silent truncation.

## 7. STEP 5 — gate results (mandatory order)

| Step | Command | Result |
| --- | --- | --- |
| 0 | `python3 .../verify_genotype_gvcf_malformed_gatk_oracle.py` **before the fix** | **exit 1** — `status: divergence`, 25 violations (`.diag/malformed-gate-before.log`) |
| a | same gate **after the fix** | **exit 0** — `status: pass`, `violations: []`, 14 cases (`native_exit=2`, `native_rows=0` for all 8 gated malformed shapes; control still byte-identical) (`.diag/malformed-gate-after.log`) |
| b1 | `python3 fastgatk-native/scripts/verify_genotype_gvcf.py` | **exit 0** — `{"status": "pass", "output_records": 1}` (`.diag/malformed-verify-gvcf.log`) |
| b2 | `ctest --test-dir fastgatk-native/build -R 'genotype-gvcf' -V` | **exit 0** — **15/15 passed**, 561.06 s (`.diag/malformed-ctest-genotype-gvcf.log`) |
| c | the 17-name strict-gate filter of the task | **exit 0** — **19/19 passed**, 1466.12 s (`.diag/malformed-strict-gates.log`) |
| d | `fastgatk-native/scripts/run_regression.sh --label malformed-report` | **exit 0** — **300/300 omp** (1401.48 s) and **300/300 serial** (1326.52 s) (`.diag/regression/20260912-142328/{omp,serial}.log`, `summary.txt`) |

Steps were run in order a → b1 → b2 → c → d by one chained job
(`.diag/malformed-steps.log`: `STEP_B1=0 STEP_B2=0 STEP_C=0 STEP_D=0`).

**Nothing was edited after the suite run**, so the binaries need no re-proof:
source `genotype_gvcf_tool.cpp` mtime 13:46:23 and the gate script 13:45:39 are
both older than the binaries (build 13:46:55, build-serial 13:47:06), which are
older than the suite start (14:23).  No rebuild or re-run was needed after the
tree last changed; the only later action was this report (a `.diag/*.md`
document no test executes).

## 8. Tree state

* The tree **CONTAINS a change**: `git status --short` =
  ` M fastgatk-native/src/genotype_gvcf_tool.cpp` and
  `?? fastgatk-native/scripts/verify_genotype_gvcf_malformed_gatk_oracle.py`;
  `git diff --numstat` = `34  4`.  No commit, no branch.
* Production code changed, so `run_regression.sh` was required and run (§7d).
  Binary hashes after the suite:
  `857060bd88fdeb738bbe8702c6a59500  fastgatk-native/build/fastgatk-genotype-gvcf`,
  `c14feaf54bf2a42554b9002a171842b6  fastgatk-native/build-serial/fastgatk-genotype-gvcf`
  (pre-fix they were `4f9943c6…` / `b179c11d…` — both trees really were rebuilt).
  Gate script md5 `c66f3a6bc6fa832f189cdc52379a1ef6`.
* NOT touched: `fastgatk-native/CMakeLists.txt`, every root `*.md`, Mutect2,
  every other tool, every other test script.
* Artifacts: probes `.diag/malformed_probe.py` (`-probe.log`),
  `.diag/malformed_reach_probe.py` (`malformed-reach.log`,
  `malformed-reach-after.log`), comparison table `.diag/malformed-table.md`;
  gate logs `.diag/malformed-gate-before.log` / `-after.log`; step logs
  `.diag/malformed-verify-gvcf.log`, `malformed-ctest-genotype-gvcf.log`,
  `malformed-strict-gates.log`, `malformed-regression.log`,
  `malformed-build.log`, `malformed-steps.log`.

## 9. What remains unproven

* The gate is **not registered** in CTest (CMakeLists.txt is off-limits); it was
  run by hand for §7a.  Until the orchestrator registers it (§4), the fix is
  not protected by the suite.
* The eight gated shapes were measured with `-V` pointing at a plain `.vcf`
  and no `-L`; the `--include-non-variant-sites` dense path shares the same
  aggregate read loop but was not measured separately.
* The BCF half of the contract is read, not measured: `bcf_read1_core` returns
  -1 (EOF) / -2 (corruption) for BCF, so `status < -1` is right there too, but
  no corrupted `.bcf` fixture was built.
* Multi-input (`-V a -V b`): pre-fix, the rest of the affected input was lost
  while the other inputs continued; post-fix the run aborts entirely.  Only the
  single-input case was measured.
* The GenomicsDB bridge path (`gendb://`) reads through a separate process and
  was not exercised.
* `control-valid-input` proves the fix does not reject a well-formed fixture;
  it is one fixture, and the evidence that it does not over-reject generally is
  the 300/300 dual-backend suite, not an exhaustive argument.
* Whether *other* tools silently truncate is reported from a grep of the idiom
  (§6), not measured — deliberately out of scope.
