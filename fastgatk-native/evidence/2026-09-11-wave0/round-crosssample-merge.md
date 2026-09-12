# Round: MULTIPLE input records at ONE position (cross-sample ALT sets)

Scope: the divergence the previous round measured but did not repair
(`.diag/round-nda.md` §6, first bullet): for two records at the same position
whose ALT sets differ, pinned GATK emits **two rows** (each with `NDA=1`) while
native emits **one** grouped row.  Gate first, fix second.

Probes: `.diag/crosssample_probe.py` (`.diag/crosssample-probe.log`),
`.diag/crosssample_realistic_probe.py` (`.diag/crosssample-realistic.log`).
Gate: the four new cases in
`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`.

## 0. Verdict

**Unsupported input shape — not a genotyping defect.  No production change.**

The divergence is reproducible (GATK 2 rows vs. native 1 row), but it exists
**only** for an input GATK does not support: two records starting at the same
locus.  GATK takes exactly one input track (§2, `GenotypeGVCFs.java:65-66`), the
only documented producer of a multi-sample GVCF collapses that shape into one
record (measured, §1.3), and in GATK's group-by-locus traversal the shape is a
**hard error** (`IllegalStateException`, exit 3, §1.2).  In its default
by-variant traversal GATK does accept it — but there it does not merge records
at all: it genotypes **each record on its own** and writes one row per record
(§2), so "cross-sample ALT union at a locus" is not a concept GenotypeGVCFs has.
Native's coalescing is the CombineGVCFs model applied to GenotypeGVCFs' input;
on every supported input the two models coincide and the rows are byte-identical
(§1.3, gated).

Consequently Step 3 was **not** attempted: matching GATK here means changing
native's record-grouping/emit policy, which the task defines as structural
(§6).  The measured GATK rows are recorded as REPORTED-ONLY gate cases (§5).

Side finding (not this round's subject, no case added): on malformed input
(a record whose FORMAT declares 4 fields but carries 5 sample values) GATK exits
3 with `TribbleException: ... There are too many keys for the sample ...` while
native exits 0, emits a warning and silently drops the record.

## 1. STEP 1 — measurements

Reference: 100 bp `chr1` of `A`s.  Fixtures are plain VCFs indexed with
`IndexFeatureFile`; both tools get identical `-R/-V/args`, native additionally
`--gatk-compatible-annotations` (the registered convention).  Fixture shorthand:
`CALL` = `0/1:20:0,20,0:100,0,100,100,100,100` (`GT:DP:AD:PL`), `MISS` = `./.:.:.:.`.

### 1.1 The diverging shape (default args)

`two-samples-two-records-different-alt`: two samples, two records at `chr1:2`,
one per sample, different ALT sets.

```
chr1	2	.	A	G,<NON_REF>	.	PASS	DP=20	GT:DP:AD:PL	CALL	MISS
chr1	2	.	A	T,<NON_REF>	.	PASS	DP=20	GT:DP:AD:PL	MISS	CALL
```

pinned GATK 4.6.2.0 (**2 rows**, exit 0):

```
chr1	2	.	A	G	92.64	.	AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63	GT:AD:DP:GQ:PL	0/1:0,20:20:99:100,0,100	./.
chr1	2	.	A	T	92.64	.	AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63	GT:AD:DP:GQ:PL	./.	0/1:0,20:20:99:100,0,100
```

native (**1 row**, exit 0):

```
chr1	2	.	A	G	92.64	.	AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63	GT:AD:DP:GQ:PL	0/1:0,20:20:99:100,0,100	./.:.:.:.:.
```

Two things to read out of that pair:

* GATK wrote **one row per input record**, in input order, and the sample the
  record carries no call for is rendered as a bare `./.` (not `./.:.:.:.:.`) —
  the row is exactly the genotyping of *that one record*.
* native kept the **first** record of the group, dropped the second (its sample
  set overlaps the first's: `genotype_gvcf_tool.cpp:6639-6646`), which loses the
  `T` allele **and** S2's call.  Native's ALT union did contain `T` — with
  `--annotate-with-num-discovered-alleles` native writes `NDA=2`, i.e. the allele
  was in the merged set — but the allele was then pruned because the only
  evidence for it (record 2) had been discarded.

### 1.2 Input-order and mode dependence of the GATK rule

| fixture | GATK | native |
| --- | --- | --- |
| §1.1, reversed input order (T record first) | **2 rows** in input order: `T` row first | 1 row, `ALT=T` (first record wins) |
| §1.1 + `--annotate-with-num-discovered-alleles` | 2 rows, each `NDA=1` | 1 row, `NDA=2` |
| §1.1 + `--include-non-variant-sites` | **exit 3**, no row | 1 row |
| one sample, two records, different ALTs | **2 rows** (both `0/1`) | 1 row, `ALT=G` |
| one sample, two BYTE-IDENTICAL records | **2 identical rows** | 1 row |
| `<NON_REF>`-only ref block + variant record at one position | 1 row (`ALT=T`, S1=`./.`, S2=`0/1`) | 1 row, `ALT=T` but S1 carries the **ref block's** call and S2 is `./.` — the two records are effectively swapped and record 2's call is lost |
| same, `--include-non-variant-sites` | **exit 3** | 1 row |

The dense-mode failure is GATK's own assertion, not a crash:

```
Caused by: java.lang.IllegalStateException: Variant input contains more than one variant starting at location: chr1:2-2
	at ...GenotypeGVCFsEngine.getVariantSubsetToProcess(GenotypeGVCFsEngine.java:362)
```

`reverse order` and `byte-identical duplicates` show the rule is *purely*
per-record: GATK's output is a function of each record alone, so the same locus
appears once per input record.

### 1.3 The supported shape, and the pipeline that produces it

`two-samples-one-record-union`: one record at `chr1:2`, ALT `G,T,<NON_REF>`,
S1 = `0/1`, S2 = `2/2`.  GATK and native rows are **byte-identical**:

```
chr1	2	.	A	G,T	277.88	.	AC=1,2;AF=0.250,0.500;AN=4;DP=40;ExcessHet=0.0000;MLEAC=1,2;MLEAF=0.250,0.500;QD=6.95	GT:AD:DP:GQ:PL	0/1:0,20,0:20:99:100,0,100,100,100,100	2/2:0,0,20:20:99:200,200,200,200,200,0
```

End-to-end through the documented pipeline, starting from two per-sample GVCFs
(`-V s1.g.vcf -V s2.g.vcf`, one record each at `chr1:2`, `G` and `T`):

```
$ gatk CombineGVCFs -R ref.fa -V s1.g.vcf -V s2.g.vcf -O combined.g.vcf      # exit 0
chr1	2	.	A	T,G,<NON_REF>	.	.	DP=40	GT:AD:DP:PL	./.:0,0,20,0:20:100,100,100,0,100,100,100,100,100,100	./.:0,20,0,0:20:100,0,100,100,100,100,100,100,100,100

$ gatk GenotypeGVCFs -R ref.fa -V combined.g.vcf --annotate-with-num-discovered-alleles ...
$ native GenotypeGVCFs -R ref.fa -V combined.g.vcf --annotate-with-num-discovered-alleles ...
chr1	2	.	A	T,G	184.52	.	AC=1,1;AF=0.250,0.250;AN=4;DP=40;ExcessHet=1.7609;MLEAC=1,1;MLEAF=0.250,0.250;NDA=2;QD=4.61	GT:AD:DP:GQ:PL	0/2:0,0,20:20:99:100,100,100,100,100,100	0/1:0,20,0:20:99:100,0,100,100,100,100
                                          ^ identical on both sides, byte for byte
```

So the two-record shape is **not reachable** from the documented workflow:
CombineGVCFs emits exactly one record for the locus, with the ALT union already
computed — which is precisely what native's coalescing assumes.

Two further measured facts about the input contract:

* GATK refuses a second `-V`: `A USER ERROR has occurred: Illegal argument
  value: Argument 'V/variant' cannot be specified more than once.` (exit 1).
  Native accepts `-V a -V b` (multi-input is its own extension, the
  CombineGVCFs role) and merges them — for that input there is no GATK oracle at
  all.
* the ref-block-plus-variant row in §1.2 shows the practical cost of native's
  rule on this unsupported shape: the surviving record's samples do not match
  GATK's, and one record's data is silently dropped.

## 2. The GATK rule, with file:line

GATK's traversal mode decides everything, and the mode is **not** a per-input
choice — it is implied by the arguments:

```java
// GenotypeGVCFs.java:284-285  (onTraversalStart)
if (!(includeNonVariants || forceOutputIntervalsPresent)) {
    changeTraversalModeToByVariant();
}
```

* **default mode (by-variant)** — one `apply()` call per input record:

  ```java
  // VariantLocusWalker.java:132-142
  apply(variant,
        Collections.singletonList(variant),      // :138
        new ReadsContext(reads, variantInterval, readFilter),
        new ReferenceContext(reference, variantInterval),
        new FeatureContext(features, variantInterval));
  ```

  `GenotypeGVCFsEngine.callRegion()` merges that one-element list
  (`GenotypeGVCFsEngine.java:128` -> `:136`) and `GenotypeGVCFs.apply()` writes
  the single resulting record (`GenotypeGVCFs.java:320-331`).  **The merge key
  is therefore the record itself**; there is no position- or allele-level
  grouping to speak of.  Two records at one position = two `apply()` calls =
  two rows, in input order.  Per-row `NDA` is the ALT count of that one record
  (`GenotypingEngine.java:464-465`), hence `NDA=1`.

* **group-by-locus mode** (`--include-non-variant-sites`, or any
  `--force-output-intervals`) — the locus stream is one base at a time
  (`VariantLocusWalker.java:188-190`, `IntervalLocusIterator`) and `apply()` gets
  **every variant overlapping** that base (`:154-172`):

  ```java
  // GenotypeGVCFsEngine.java:349-368  (getVariantSubsetToProcess)
  final List<VariantContext> matchingStart =
          preProcessedVariants.stream().filter(vc -> vc.getStart() == loc.getStart()).collect(Collectors.toList());
  if (matchingStart.size() == 0)      { return preProcessedVariants; }
  else if (matchingStart.size() == 1) { return matchingStart; }
  // since this tool only accepts a single input source, there should never be
  // more than one variant at a given starting locus            // :359-360
  throw new IllegalStateException(String.format(
          "Variant input contains more than one variant starting at location: %s", ...));   // :361-364
  ```

  Here the key **is** position — but only one record may *start* at the locus;
  every other member of the list must be a record that started earlier and spans
  it.  The actual cross-record merge
  (`ReferenceConfidenceVariantContextMerger.merge()`, called at
  `GenotypeGVCFsEngine.java:136`) is a locus-level union of the alleles:
  `isSpanningEvent = loc.getStart() != vc.getStart()` (`:150`) and the union is
  built by `collectTargetAlleles()` (`:325-348`).  **It is unreachable for two
  records that start at the locus**, because that case throws first.

The documented input settles which mode matters:

```
// GenotypeGVCFs.java:65-66
* The GATK4 GenotypeGVCFs tool can take only one input track.  Options are 1) a single single-sample GVCF 2) a single
* multi-sample GVCF created by CombineGVCFs or 3) a GenomicsDB workspace created by GenomicsDBImport.
```

and the supported multi-sample producer emits one record per locus
(`CombineGVCFs extends MultiVariantWalkerGroupedOnStart`, `CombineGVCFs.java:82`;
one `vcfWriter.add(mergedVC)` per locus at `:412-420`) — measured in §1.3.

## 3. Native's rule, with file:line

```cpp
// fastgatk-native/src/genotype_gvcf_tool.cpp:1180-1187
std::string record_key(const bcf1_t* record) {
    // Records at one locus must be coalesced even when different input
    // shards carried different concrete ALT subsets.  The union/remapping is
    // performed after grouping, so ALT names do not belong in this key.
    const char* ref = record->n_allele > 0 ? record->d.allele[0] : "";
    return std::to_string(record->rid) + ":" + std::to_string(record->pos) + ":" + ref;
}
```

* merge key = **`rid:pos:REF`** (ALT names deliberately excluded);
* the aggregate path groups on that key (`:6541`), builds the concrete ALT union
  (`:6546-6560`) and emits **one record per group**; within a group it keeps the
  first record and **discards any later record whose sample set overlaps the
  accepted ones** (`:6639-6646`);
* the streaming path (`--stream-by-locus`) uses the identical rule
  (`:5869-5880` group by `record.key`, `:5969-5972` overlap drop) — verified by
  running both paths on §1.1: one row each, byte-identical between the paths;
* multiple `-V` inputs are accepted (`:592-593`) and merged by the same key.

**How the two differ.**  GATK (default) groups by *input record* and emits one
row per record; native groups by *locus* and emits one row per locus with the
sample sets unioned.  They agree exactly when the input has one record per
locus (the documented input shape) — and there GATK's own group-by-locus mode
would also merge.  They can only disagree when ≥2 records share
`(chrom, pos, REF)`, which is the shape GATK's engine asserts cannot occur.

## 4. STEP 2 — the gate

Changed file: `fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`
(+204/-1).  Nothing was added to the repo other than `.diag/*` artifacts.

Added one **gated** case and three **REPORTED-ONLY** cases, all on the existing
byte-exact row comparison:

| case | gated | what it pins |
| --- | --- | --- |
| `cross-sample-alt-union-single-record` | **yes** | the supported shape: one record, ALT union `G,T`, two samples — GATK row pinned and required byte-identical to native's |
| `same-position-two-records-cross-sample-alt` | no | the two measured GATK rows for §1.1 (via `expect`), native's 1 row printed |
| `same-position-two-records-cross-sample-alt-nda` | no | the same, with `--annotate-with-num-discovered-alleles`: GATK 2 rows `NDA=1` each, native 1 row `NDA=2` |
| `same-position-two-records-dense-unsupported` | no | `--include-non-variant-sites`: GATK exit 3 and no row (pinned via the new per-case `gatk_expected_exit`), native 1 row |

The `expect` constants make each reported-only case fail loudly (and print the
moved rows) if pinned GATK ever changes, without affecting the exit status; the
gated case fails the gate if native regresses on the supported shape.  The only
harness change is honouring an optional per-case `gatk_expected_exit` so that
the deliberately-failing GATK run in the last case is an observation rather than
a violation (default 0: no behaviour change for the 37 pre-existing cases).

Exit status: **0** — 41 cases, `"status": "pass"`, `"violations": []`
(`.diag/cs-gate-full.log`).  Before this round the shape was not gated at all,
so there is no "failing gate" to show; the pre-fix evidence is the measurement
itself (§1) and the runtime probe logs.

## 5. STEP 3 — the fix: not attempted, deliberately

The only way to reproduce GATK's rows here is to make native stop coalescing
records that start at the same locus, i.e. to emit one output row per **input
record** in input order with non-carried samples rendered `./.`.  That is a
change to the record-grouping/emit policy (`record_key` at
`genotype_gvcf_tool.cpp:1180-1187`, the group loops at `:6538-6650` and
`:5859-5980`) and it would remove the cross-sample merge that native offers for
`-V a -V b` input — exactly the "structural change to native's
record-grouping/merge key" the task says to stop and report instead of
attempting.  It would also be a behaviour change for every input that *benefits*
from the coalescing, i.e. a large blast radius for a shape GATK declares
invalid.

**Recommendation (not implemented).**

1. Keep the coalescing; it is the CombineGVCFs model and it is exact on the
   documented input (gated by `cross-sample-alt-union-single-record`).
2. If this shape ever needs to be supported, do it as an explicit
   *by-variant/1:1 traversal mode* mirroring `GenotypeGVCFs.java:284-285` (one
   record in, one row out, input order), never as a silent change of the default
   grouping — and consider whether `--include-non-variant-sites` should mirror
   GATK's hard error on >1 record starting at a locus.
3. Independent of parity: a group that discards a record because its samples
   overlap is silent data loss (measured above: `T` and S2's call vanish).  A
   diagnostic (stderr warning + a counter in the telemetry JSON) would make the
   shape visible without changing any byte of the output.  Not done here: the
   task forbids touching the merge policy, and this would need its own round.

## 6. STEP 4 — stale assertions

**None; no test line was corrected.**  No registered test builds a fixture with
two records at the same position (an AST scan of every
`fastgatk-native/scripts/*genotype*.py` for a single fixture string containing
two data lines with the same `(chrom, pos)`: 0 hits), and no registered test
invokes native GenotypeGVCFs with two `-V` (AST scan for ≥2 `-V` arguments in
one command: 0 hits).  So nothing asserted the old behaviour and nothing had to
be updated; the only test-side change is the four cases in §4.

## 7. STEP 5 — gate results (mandatory order)

| Step | Command | Result |
| --- | --- | --- |
| a | `python3 fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` | **exit 0** — 41 cases, 0 violations (`status: pass`); 5 of them new (§4) (`.diag/cs-gate-full.log`) |
| b1 | `python3 fastgatk-native/scripts/verify_genotype_gvcf.py` | **exit 0** — `{"status": "pass", "output_records": 1}` (`.diag/cs-verify-gvcf.log`) |
| b2 | `ctest --test-dir fastgatk-native/build -R 'genotype-gvcf' -V` | **exit 0** — **15/15 passed**, 579.18 s (`.diag/cs-ctest-genotype-gvcf.log`); includes `fastgatk-genotype-gvcf-contract`, the native HC → GenotypeGVCFs end-to-end path |
| c | the 17-name strict-gate filter of the task | **exit 0** — **19/19 passed**, 1419.48 s (`.diag/cs-strict-gates.log`) |
| d | `run_regression.sh` | **not run — not required**: no production code and no compiled byte changed (§8), so the task's condition ("if you changed production code") is not met |

### 7.1 Dual-backend regression (d)

Not applicable; §8 gives the binary hashes proving nothing was rebuilt.

### 7.2 The one edit made after the suite, and why it did not invalidate it

After the suite finished I corrected **one misleading sentence** in the docstring
I had added (it said native "keeps the first record of each group whose samples
overlap an already-accepted record"; the rule is the opposite — it keeps the
first record and *drops* later overlapping ones).  It is comment text inside the
Python gate script; no compiled byte exists for it.  Proof by re-run instead of
by rebuild:

```
python3 fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py   -> exit 0, 41 cases, 0 violations
ctest --test-dir fastgatk-native/build -R 'spandel-gatk-oracle' -V           -> exit 0, 1/1 passed, 371.72 s
```

The post-edit gate payload is byte-identical to the pre-edit one once the random
temporary-directory name is normalised away (all 41 cases' rows, exits and
violations; verified in `.diag/cs-gate-full-after.log` vs `.diag/cs-gate-full.log`).
Steps 5(b1), 5(b2) and 5(c) were **not** re-run: no script those tests execute
other than the spandel oracle was touched, and the binaries are unchanged.

## 8. Tree state

* The tree **CONTAINS a change**: one modified file, no commit, no branch.
  `git status --short` = ` M fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`;
  `git diff --numstat` = `204  1`; the file's md5 is
  `4887f3cc7c40a640e2dff85f5ef8b226` (the hash the post-suite re-run of §7.2 used).
* **No production code was touched** — `fastgatk-native/src/genotype_gvcf_tool.cpp`
  is byte-identical to the previous round's commit, no rebuild happened, and the
  binaries are unchanged:
  `4f9943c601227c1f8c1e49513f55e762  fastgatk-native/build/fastgatk-genotype-gvcf`,
  `b179c11d01fd6b154e18d250eb90f0ec  fastgatk-native/build-serial/fastgatk-genotype-gvcf`.
* NOT touched: `fastgatk-native/CMakeLists.txt`, every root `*.md`, Mutect2,
  every other tool, every other test script.
* Artifacts: probes `.diag/crosssample_probe.py` (`.diag/crosssample-probe.log`),
  `.diag/crosssample_realistic_probe.py` (`.diag/crosssample-realistic.log`);
  gate logs `.diag/cs-gate-full.log` (+ `.diag/cs-gate-full-after.log`, §7.2),
  `.diag/cs-verify-gvcf.log`, `.diag/cs-ctest-genotype-gvcf.log`,
  `.diag/cs-strict-gates.log`, `.diag/cs-ctest-spandel-after.log`.

## 9. What remains unproven

* The GenomicsDB half of the input contract (`GenotypeGVCFs.java:65-66`) is read,
  not measured: no GenomicsDB workspace was built.  Everything else in §1-§3 is
  measured.
* "The shape is unreachable from the documented pipeline" is proven for
  CombineGVCFs (measured) and argued from the docstring for GenomicsDBImport.
* The scan behind §6 is structural (AST), not exhaustive: a test could build
  such a fixture dynamically at run time.  The empirical check is step 5.
* Native's coalescing is a no-op only if no *supported* input contains two
  records sharing `(chrom, pos, REF)`.  That is true for the shape GATK
  documents, and step 5(b2) exercises the native HC → GenotypeGVCFs end-to-end
  path (`fastgatk-genotype-gvcf-contract`), but no systematic sweep of every
  fixture in the tree was run.
* `--force-output-intervals` (the second way to enter group-by-locus mode) was
  not exercised on this shape; only `--include-non-variant-sites` was.
* The malformed-FORMAT divergence noted in §0 (GATK exit 3 vs. native exit 0,
  record dropped) is measured but not gated; its own input is malformed VCF, so
  it was left alone this round.
