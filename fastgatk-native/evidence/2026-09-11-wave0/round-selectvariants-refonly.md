# Round: SelectVariants reference-only record retention (native drop vs pinned GATK 4.6.2.0)

Date: 2026-09-11. Working tree: `22f4f51` + uncommitted changes (no commit, no branches created).
Scope: `fastgatk-select-variants` only. Mutect2 and every other tool untouched.

## 1. The divergence

`fastgatk-native/scripts/verify_select_variants.py` (pre-round line 820) asserted

```python
assert ref_only_records == []
```

for this fixture (that file's own fixture, declared with **no `##INFO` lines at all**) and these arguments:

```
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1
chr1   20  .  A   C,G 50   PASS   .    GT:AD:PL:GQ  0/0:30,0,0:0,30,60,40,70,80:30

gatk SelectVariants -V ref-only-unused.vcf.gz -O out.vcf.gz --remove-unused-alternates
```

Native dropped the record; GATK keeps it. The audit row already recorded this
(`IMPLEMENTATION_STATUS.md:586`: *"`verify_select_variants.py:820` | `ref_only_records == []`（丢弃全 hom-ref 记录） | **保留**，ALT='.'"*) —
this round re-measured it independently, localised the native code, fixed it and gated it.

Measurement note: pinned GATK refuses the test's input as written, because Python's
`gzip` module does not produce BGZF: `IndexFeatureFile` fails with *"Input file is not in
valid block compressed format"* and `SelectVariants` then fails with *"An index is
required but was not found"*. All measurements below therefore use an **uncompressed
`.vcf`** plus GATK's own `IndexFeatureFile` (Tribble index), which both engines read;
the byte content of the records is identical to the test's fixture.

## 2. Measured GATK truth (STEP 1) — option by option

Pinned jar: `third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar`,
`third_party/jdk17/bin/java`. Same fixture, four argument sets (all exit 0):

| arguments | GATK data rows |
| --- | --- |
| *(none)* | `chr1 20 . A C,G 50 PASS . GT:AD:GQ:PL 0/0:30,0,0:30:0,30,60,40,70,80` |
| `--remove-unused-alternates` | `chr1 20 . A . 50 PASS AN=2 GT:AD:GQ:PL 0/0:30:30:0` |
| `--exclude-non-variants` | *(none — record dropped)* |
| `--remove-unused-alternates --exclude-non-variants` | *(none — record dropped)* |

**The retention behaviour is option-dependent, and the deciding option is
`--exclude-non-variants`, not `--remove-unused-alternates`.** GATK keeps the all-hom-ref
record and writes `ALT='.'` under `--remove-unused-alternates`; it drops the same record
when `--exclude-non-variants` is added (the flag defaults to false).

GATK source for the rule (local `gatk-source/`):

* `SelectVariants.java:265-266` — `@Argument(fullName="exclude-non-variants", …) private boolean excludeNonVariants = false;`
* `SelectVariants.java:696` — `VariantContext result = subsetGenotypesBySampleNames(vc, preserveAlleles, removeUnusedAlternates);`
* `SelectVariants.java:1233` — `return preserveAlleles ? subset : GATKVariantContextUtils.trimAlleles(subset,true,true);`
* `GATKVariantContextUtils.java:1455-1460` — `trimAlleles(...)`: `if ( inputVC.getNAlleles() <= 1 … ) return inputVC;` — the trimmed context keeps the **reference allele alone**, it is not a dropped record.
* `SelectVariants.java:708-715` — the non-variant removal, which runs *after* subsetting/trimming: `if (excludeNonVariants) { final boolean nonVariant = ! result.isPolymorphicInSamples() || GATKVariantContextUtils.isSpanningDeletionOnly(result); if (nonVariant) return; }`

Deeper rules, all measured on the same fixture family (record `A C,G`, `GT:AD:GQ:PL`):

| probe | GATK row under `--remove-unused-alternates` | rule |
| --- | --- | --- |
| base (`GQ=30`, PL `0,30,60,40,70,80`) | `… AN=2 GT:AD:GQ:PL 0/0:30:30:0` | PL/AD collapse to one allele; AN refreshed |
| `GQ=99` in the input | `… AN=2 GT:AD:GQ:PL 0/0:30:99:0` | **original GQ is preserved**, not re-derived: `AlleleSubsettingUtils.java:94-100` *"if we subset to just ref allele, keep the GQ"* |
| `GT=./.` (all no-call) | `… AN=0 GT:AD:GQ:PL ./.:30:30:0` | `AN=0` is written explicitly |
| input `INFO=AC=0,0;AN=2;AF=0,0` | `… AN=2 …` | AC/AF are cleared (no ALT left), AN survives |
| input ALREADY `ALT='.'`, `INFO='.'` | `chr1 20 . A . 50 PASS . GT:AD:GQ:PL 0/0:30:30:0` | no trimming happened → GATK adds **nothing**, not even AN |
| GQ declared but absent in the record | `… AN=2 GT:AD:PL 0/0:30:0` | GQ emitted only when the input genotype had one (`AlleleSubsettingUtils.java:128`) |
| haploid `GT=0` | `… AN=1 GT:AD:GQ:PL 0:30:30:0` | the PL collapse is ploidy-aware |
| two records (POS 10 variant + POS 20 ref-only) | POS 20 row identical to the single-record case | interaction with the write path |

Ordering fact (measured, driven by the select-variants ctest subset, see §4b): GATK's order
is **trim (`:696`) → filtered genotypes to no-call (`:698`) → non-variant test (`:708-715`)**.
A record whose only ALT carrier is a filtered genotype therefore still carries its ALT at
trim time.

Fixture/format note: with the test's `##INFO`-less header, GATK *adds*
`##INFO=<ID=AC|AF|AN|DP>` to the output header and writes `AN=2`; both tools' data rows are
byte-identical for the canonical FORMAT order (see §6 for what is *not* identical).

## 3. STRICT GATE (STEP 2)

New file: `fastgatk-native/scripts/verify_select_variants_refonly_gatk_oracle.py`
(house conventions: `main() -> int`, `argparse --native` + `FASTGATK_SELECT_VARIANTS_BINARY`
default, `TemporaryDirectory` scratch, final JSON status line, `oracle_guard`, modelled on
`verify_variant_recalibrator_culprit_gatk_oracle.py`).

Design: every case builds its fixture in a fresh temp dir, runs pinned GATK and native on the
**same** input file and arguments, and requires the **data rows to be byte-identical**, plus a
fixture-validity gate that the GATK side still reproduces the rows measured here (so the gate
cannot pass vacuously if GATK changes). 13 gated cases, 1 reported-only case
(`--expect-divergence` mode available).

Exit status, before vs after the fix — same script, same fixtures:

* **before** (pre-change binary, rebuilt from `git show HEAD:…/select_variants_tool.cpp`):
  `EXIT=1`, log `.diag/refonly-gate-prefix-binary.log`:

```
[rua-keeps-ref-only] gated args=['--remove-unused-alternates']
    GATK   rows: ['chr1\t20\t.\tA\t.\t50\tPASS\tAN=2\tGT:AD:GQ:PL\t0/0:30:30:0']
    NATIVE rows: []
    VIOLATION: rua-keeps-ref-only: record count differs: GATK=1 NATIVE=0
[rua-preserves-original-gq] GATK ['chr1\t20\t.\tA\t.\t50\tPASS\tAN=2\tGT:AD:GQ:PL\t0/0:30:99:0'] NATIVE []
[rua-mixed-file] GATK POS20 ['chr1\t20\t.\tA\t.\t50\tPASS\tAN=2\tGT:AD:GQ:PL\t0/0:30:30:0'] NATIVE []
```

* **after**: `EXIT=0`, 13/13 gated cases, no violations, log `.diag/refonly-gate-final.log`
  (`{"status": "pass", … "gated_cases": 13, "reported_only_cases": ["diagnostic-nocall-before-trim"], "violations": []}`).

Registration: `fastgatk-native/CMakeLists.txt` was **not** touched (hard constraint 1); the
orchestrator registers gates. The script is runnable standalone and exits non-zero on mismatch.

## 4. FIX (STEP 3) — minimal, localised

### 4a. Localisation (what dropped the record)

Pre-change `fastgatk-native/src/select_variants_tool.cpp:2390-2400` (that is `HEAD`, quoted
verbatim):

```cpp
                if (kept.size() < static_cast<std::size_t>(record->n_allele)) {
                    // GATK drops records whose samples do not use any ALT
                    // when --remove-unused-alternates is enabled.  A
                    // ref-only record has no valid Number=G/Number=R target
                    // space, so do this contract check before invoking the
                    // Kokkos remap kernels.
                    if (kept.size() == 1) {
                        ++filtered_records;
                        bcf_clear(record);
                        continue;
                    }
```

The comment was a **native-only assumption**: `GATKVariantContextUtils.java:1455` returns a
single-allele context rather than nothing, and `SelectVariants.java:708-715` is where the
record actually disappears — only under `--exclude-non-variants`.

### 4b. What was changed (three localised edits, no restructuring)

1. **`src/select_variants_tool.cpp:1956-2010` (new helper `trim_record_to_reference_only`)** —
   writes the trimmed record as reference-only (`bcf_update_alleles_str` with the REF allele,
   so the writer emits `ALT='.'`), compacts Number=R/G via the existing `compact_fields`, then
   applies the measured GATK payload contract: original GQ preserved when the input genotype
   had one and dropped when it did not (`AlleleSubsettingUtils.java:94-100`, `:128`), AC/AF
   cleared and AN written (including `AN=0`) for the no-ALT allele set.
2. **`src/select_variants_tool.cpp:1885-1897` (`compact_fields` PL branch)** — a new
   `new_alleles == 1` arm collapses PL to the single all-ref cell (old PL index 0, per
   `AlleleSubsettingUtils.subsettedPLIndices`, `:426-441`) and leaves GQ to the caller. The
   Kokkos kernel `fastgatk-kernels/src/genotype.cpp:797-800` rejects
   `target_allele_count < 2` (*"invalid genotype PL remap dimensions"*) — the shared kernel was
   deliberately **not** changed, since other tools use it.
3. **`src/select_variants_tool.cpp:2469-2507`** — the drop branch is replaced by: if
   `--exclude-non-variants` is set, count `non_variant_skipped` and drop (the post-trim
   non-variant test GATK performs at `:708-715`); otherwise emit the reference-only record and
   fall through to the normal write path (so `--sites-only-vcf-output`,
   `--keep-original-ac/--keep-original-dp` and the writer keep their existing behaviour). The
   multi-ALT path is unchanged, just moved into an `else`.

Edit 3 was forced by the ctest subset: with only edits 1+2, `fastgatk-select-variants-filtered-nocall-gatk-oracle`
failed, because native evaluated `--exclude-non-variants` **before** `--set-filtered-gt-to-nocall`
and trimming, so a record that becomes non-variant only through the no-call conversion was
resurrected. GATK evaluates it after both. A gate case
(`rua-env-drops-after-nocall-conversion`) now pins that.

Full diff: `.diag/refonly-source.diff` (one file changed: +142/−48, of which most of the
deletions are the re-indented `else` block).

## 5. STALE ASSERTION CORRECTED (STEP 4)

File `fastgatk-native/scripts/verify_select_variants.py`, diff in `.diag/refonly-test.diff`.
Every line touched (post-edit numbering):

| lines | change |
| --- | --- |
| 798-815 | comment rewritten: states the measured GATK contract + `SelectVariants.java:1195-1233`, `GATKVariantContextUtils.java:1455`, `SelectVariants.java:708-715`, `:265-266`; records that the previous expectation was native-only; records the literal measured GATK row `chr1 20 . A . 50 PASS AN=2 GT:AD:GQ:PL 0/0:30:30:0` |
| 836-837 | `assert ref_only_records == []` → `assert ref_only_records == ["chr1\t20\t.\tA\t.\t50\tPASS\tAN=2\tGT:AD:PL:GQ\t0/0:30:0:30"]` (exact record, not a looser check) |
| 838-842 | added field-level asserts on the same row: `ALT == "."`, `INFO == "AN=2"`, `FORMAT == "GT:AD:PL:GQ"`, sample `"0/0:30:0:30"` |
| 844-852 | added the option-dependence check: the same input with `--remove-unused-alternates --exclude-non-variants` must emit **no** data rows |

Only that assertion block depended on the dropped record; nothing else in the file needed a
change (the rest of the file passed untouched, and `non_variant_skipped == 4` at line 231 is a
different invocation without `--remove-unused-alternates`, so it is unaffected). The other
`verify_select_variants_*_oracle.py` scripts needed no correction: the only other
`--remove-unused-alternates` user, `verify_select_variants_filtered_nocall_oracle.py`, always
combines it with `--exclude-non-variants` and its expectation (rows 1 and 3 removed) is what
GATK does — it now passes again after edit 3. No assertion in this file was weakened: the
record-count increase `[]` → one exact row is strictly stronger.

Known-value caveat recorded in the test comment: GATK re-emits FORMAT keys in its own rebuild
order (`GT:AD:GQ:PL`) while native preserves the input header order (`GT:AD:PL:GQ`); the
asserted values are identical (`AD=30`, `PL=0`, `GQ=30`). That key-order difference is
pre-existing, unrelated to the drop, and is measured + reported (not gated) by the new oracle
(`diagnostic-input-format-order`).

## 6. STEP 5 — gate results (mandatory order)

| step | command | exit | result |
| --- | --- | --- | --- |
| 5a | `python3 fastgatk-native/scripts/verify_select_variants_refonly_gatk_oracle.py --native build/fastgatk-select-variants` | **0** | 13/13 gated cases byte-identical; 1 reported-only case (`.diag/refonly-gate-final.log`) |
| 5b | `python3 fastgatk-native/scripts/verify_select_variants.py` (omp binary, and serial binary with `FASTGATK_EXPECTED_EXECUTION_SPACE=Serial`) | **0** / **0** | `{"status": "pass", "output_records": 1, "gatk_oracle_exact": true}` |
| 5b | `ctest --test-dir fastgatk-native/build -R 'select-variants' -V` | **0** | 5/5 passed (`.diag/step5b-select-variants.log`) |
| 5c | `ctest -R 'window-invariance-gatk-oracle\|alleles-overlap-gate-oracle\|span-del-qual-gatk-oracle\|gvcf-symbolic-prior-gatk-oracle\|af-zero-format-gatk-oracle\|arbitrary-ploidy-span-del-prior-oracle\|polyploid-gvcf-span-del-prior-oracle\|spanning-prior-genotype-gq-oracle\|multialt-owner-annotation-oracle\|mutect2-recheck\|gvcf-indel-end-gatk-oracle\|culprit-gatk-oracle' -V` | **0** | 14/14 passed (the filter matches 14 registered tests) (`.diag/step5c-strict-gates.log`) |
| 5d | `fastgatk-native/scripts/run_regression.sh --label select-variants-refonly` (after rebuilding the changed tool in **both** trees) | **0** | omp **295/295**, serial **295/295**; evidence `.diag/regression/20260911-165816/`, log `.diag/step5d-regression.log` |

No criterion failed, so nothing had to be reverted. The full suite was run exactly once.

Two environment notes on 5d, both independent of this change:
`run_regression.sh` printed *"源码比二进制新（…/variant_recalibrator_tool.cpp）"* for both
backends — that file was already newer than the binaries before this round (another agent's
in-flight work, not touched here), and it is a warning, not a failure; the script's snapshot
also reports the tree as `22f4f51` + 3 uncommitted changes.

## 7. Change state

The tree **CONTAINS** the change (no commit, no branch, nothing reverted):

* `fastgatk-native/src/select_variants_tool.cpp` — modified (fix, `M`)
* `fastgatk-native/scripts/verify_select_variants.py` — modified (corrected assertion, `M`)
* `fastgatk-native/scripts/verify_select_variants_refonly_gatk_oracle.py` — new gate (`??`)
* rebuilt binaries: `fastgatk-native/build/fastgatk-select-variants` (OpenMP) and
  `fastgatk-native/build-serial/fastgatk-select-variants` (Serial), both built from the current
  source after the final edit.
* `.diag/` artifacts: `round-selectvariants-refonly.md` (this file),
  `refonly-source.diff`, `refonly-test.diff`, `refonly-gate-prefix-binary.log`,
  `refonly-gate-final.log`, `step5b-select-variants.log`, `step5c-strict-gates.log`,
  `step5d-regression.log`, `prefix-evidence.sh`, `select_variants_tool.cpp.fixed-backup`.
* Untouched, as required: `fastgatk-native/CMakeLists.txt`, every `*.md` at the repository
  root, Mutect2 and all other tools.

## 8. What remains unproven / deliberately out of scope

Measured but **not** fixed (reported by the gate, never gated):

1. **FORMAT key order after allele subsetting.** GATK emits `GT:AD:GQ:PL`; native preserves
   the input header order. Visible on the triploid surviving record too
   (`0/1/1:30,8:10:0,10,20,30` vs `0/1/1:30,8:0,10,20,30:10`), which this round did not touch.
   Reproducing GATK's order means reordering FORMAT fields per record — not a minimal edit,
   and it would invalidate assertions in `verify_select_variants.py:790-795` that pin native's
   order deliberately.
2. **INFO key order / AF precision on records that still carry ALTs.** GATK writes
   `AC=1;AF=0.5;AN=2` (and `AF=0.667` where native writes `0.666667`); native keeps the input
   order. Case `rua-mixed-file` gates the reference-only row byte-for-byte and reports the
   sibling row.
3. **`--set-filtered-gt-to-nocall` before/after trimming order.** Case
   `diagnostic-nocall-before-trim` (reported-only): without `--exclude-non-variants`, GATK keeps
   `ALT=C` because it trims while the filtered genotype still calls it, while native converts to
   no-call first and trims the allele away. Both now emit one record. Fixing the order is a
   separate change to the `--set-filtered-gt-to-nocall` contract and was not attempted here.

Not measured at all (no claim made):

* `--keep-original-ac` / `--keep-original-dp` on a reference-only record: the new path falls
  through to the existing restore logic; no GATK comparison was run. *Speculation:* GATK would
  restore the original AC/AN values it stored, which for a ref-only record has not been probed.
* `--sites-only-vcf-output` combined with the new reference-only path (it goes through the same
  common write path, so `bcf_subset_format` still applies, but no oracle case exercises it).
* Records with symbolic ALTs (`<NON_REF>`, `*`) trimmed to reference-only: `is_non_variant`
  treats a spanning-deletion-only record as non-variant (`GATKVariantContextUtils.isSpanningDeletionOnly`
  in `SelectVariants.java:711`); native's post-trim check uses `kept.size() == 1` only, so a
  `*`-only survivor under `--exclude-non-variants` is **unverified**.
* GVCF/NON_REF inputs feeding `--remove-unused-alternates`, and multi-sample ref-only records
  (the gated fixtures use one or two samples).
* The new gate is not registered in CTest (constraint 1: `CMakeLists.txt` untouched), so its
  result is currently only as durable as this round's logs.
