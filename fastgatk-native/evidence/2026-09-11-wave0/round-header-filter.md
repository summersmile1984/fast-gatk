# Round: the output header's `##FILTER` line (GATK always declares `LowQual`, native did not)

Scope: the divergence the previous round measured and deliberately left ungated
(`.diag/round-filter-rgq.md` §7, last bullet): on the fixture of
`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`, pinned GATK
4.6.2.0 GenotypeGVCFs writes `##FILTER=<ID=LowQual,Description="Low quality">` in
its output header while native's `--gatk-compatible-annotations` header had **no**
`##FILTER` line unless the input declared one.  The existing oracles compare data
rows only, so no gate could see it.  Gate first, fix second.

## 0. Verdict

* **GATK truth (measured, pinned 4.6.2.0):** GenotypeGVCFs declares
  `##FILTER=<ID=LowQual,Description="Low quality">` in its output header for
  **every** run, whether or not the input declares any filter.  Measured at
  header index 2 of 31 (plain input) / index 2 of 32 (input also declares
  `q10`) / index 3 of 32 (input declares `AAA`, which sorts first).  The text,
  the id and the description are byte-identical to the input-less case in all
  fixtures measured.
* **GATK rule with source:** `headerLines.add(GATKVCFHeaderLines.getFilterLine(
  GATKVCFConstants.LOW_QUAL_FILTER_NAME));` is the **last** `add()` of
  `GenotypeGVCFsEngine.setupVCFWriter()`
  (`gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/GenotypeGVCFsEngine.java:416`),
  immediately before `new VCFHeader(headerLines, ...)` (`:418-419`) and
  `vcfWriter.writeHeader(outputHeader)` (`:420`).  Nothing guards it: it is not
  conditional on the input header, on `--include-non-variant-sites`, on
  `--keep-combined`, on `--dbsnp` or on the output mode, and `setupVCFWriter()`
  has exactly one caller (`GenotypeGVCFs.java:305`).  The text comes from
  `GATKVCFHeaderLines.java:89`
  (`addFilterLine(new VCFFilterHeaderLine(LOW_QUAL_FILTER_NAME, "Low quality"))`)
  with `LOW_QUAL_FILTER_NAME = "LowQual"` (`GATKVCFConstants.java:179`).  The
  sibling tool documents the same intent explicitly:
  *"FILTER fields are added unconditionally as it's not always 100% certain the
  circumstances where the filters are used. For example, in emitting all sites
  the lowQual field is used"* (`HaplotypeCallerEngine.java:557-559`).
* **Fix kept in the tree:** one local change in
  `fastgatk-native/src/genotype_gvcf_tool.cpp`
  (`gatk_compatible_header_text`): the `##FILTER` group is now collected, given
  GATK's unconditional `LowQual` declaration when absent, sorted the way htsjdk
  sorts it, and re-inserted where the input's filter lines were — the same
  collect/sort/re-insert shape the function already used for `##FORMAT`.
  46 changed lines in that file (`+45/-1`), of which 25 are code (including
  braces) and 20 are comment; no restructuring, no kernel ABI, no
  `CMakeLists.txt`, no Mutect2, no other tool.
* Oracle: **exit 1 before** the fix (22 `##FILTER`-only violations, literal
  header lines in §2.1), **exit 0 after** (25 cases: 23 gated + 2 reported-only,
  0 violations).  One new **GATED** case pins the group *order* as well as the
  line's presence.
* Step 5: (a) exit 0; (b1) exit 0, (b2) 15/15 exit 0; (c) 19/19 exit 0;
  (d) **300/300 on both backends**, `REG_EXIT=0`.  Change **kept in the tree**.

## 1. STEP 1 — measured truth, byte for byte

Reproducer: `.diag/header_filter_probe.py` — log and JSON
`.diag/header-filter-probe.log` / `.diag/header-filter-probe.json` before the fix
(6 fixtures) and `.diag/header-filter-probe-after.log` /
`.diag/header-filter-probe-after.json` after it (8 fixtures).  It builds each
fixture with the
oracle's own 100 bp `chr1` reference, indexes it with `gatk IndexFeatureFile`,
and runs

```
GATK:   java -jar gatk-package-4.6.2.0-local.jar GenotypeGVCFs -R ref.fa -V in.g.vcf \
          [-all-sites] -O out.vcf --create-output-variant-index false
native: fastgatk-genotype-gvcf -R ref.fa -V in.g.vcf [-all-sites] \
          --gatk-compatible-annotations -O out.vcf
```

dumping **every** header line of both outputs in order.

### 1.1 The `##FILTER` lines, both sides, all fixtures measured

| fixture (input header) | mode | GATK `##FILTER` (index) | native BEFORE (index) | native AFTER (index) |
| --- | --- | --- | --- | --- |
| no `##FILTER` line at all | default | `LowQual` (2 of 31) | *(none)* | `LowQual` (3 of 25) |
| no `##FILTER` line at all | `-all-sites` | `LowQual` (2 of 31) | *(none)* | `LowQual` (3 of 25) |
| `##FILTER=<ID=LowQual,Description="Low quality">` | default | `LowQual` (2 of 31) | `LowQual` (3 of 25) | `LowQual` (3 of 25) |
| `##FILTER=<ID=q10,Description="Quality below 10">` | default | `LowQual`, `q10` (2,3 of 32) | `q10` (3) | `LowQual`, `q10` (3,4 of 26) |
| `##FILTER=<ID=AAA,Description="A filter">` | default | `AAA`, `LowQual` (2,3 of 32) | `AAA` (3) | `AAA`, `LowQual` (3,4 of 26) |
| `##FILTER=<ID=PASS,Description="All filters passed">` | default | `LowQual`, **`PASS`** (2,3 of 32) | `LowQual` (3) | `LowQual` (3 of 25) — **unchanged, still divergent** (§7.3) |
| no `##FILTER` line, ordinary variant path | default | `LowQual` (2 of 31) | *(none)* | `LowQual` (3 of 25) |
| no `##FILTER` line, ordinary variant path | `-all-sites` | `LowQual` (2 of 31) | *(none)* | `LowQual` (3 of 25) |

The line is byte-identical on both sides where it exists, including the
description text and the quoting:

```
##FILTER=<ID=LowQual,Description="Low quality">
```

### 1.2 Full header order (the fixture of the gate, `-all-sites`, no input filter)

GATK, 31 lines — `##fileformat`, then **sorted by key** (`ALT` < `FILTER` <
`FORMAT` < `GATKCommandLine` < `INFO` < `contig` < `source`), values sorted
inside each key:

```
[ 0] ##fileformat=VCFv4.2
[ 1] ##ALT=<ID=NON_REF,Description="Represents any possible alternate allele">
[ 2] ##FILTER=<ID=LowQual,Description="Low quality">          <-- the divergence
[ 3] ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
[ 4] ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
[ 5] ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
[ 6] ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
[ 7] ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
[ 8] ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
[ 9] ##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, ...">
[10] ##GATKCommandLine=<ID=GenotypeGVCFs,CommandLine="...",Version="4.6.2.0",Date="...">
[11] ##INFO=<ID=AC,Number=A,...>            ... [12] AD ... [13] AF ... [14] AN ...
[15] ##INFO=<ID=BaseQRankSum,...>           [16] ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
[17] ##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
[18] ExcessHet  [19] FS  [20] InbreedingCoeff  [21] MLEAC  [22] MLEAF  [23] MQ
[24] MQRankSum  [25] QD  [26] ReadPosRankSum  [27] SOR
[28] ##contig=<ID=chr1,length=100>
[29] ##source=GenotypeGVCFs
[30] #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	STAR
```

native AFTER, 25 lines — input order preserved, `##contig` left where the input
put it, the tool's own INFO rank order, no `##GATKCommandLine`:

```
[ 0] ##fileformat=VCFv4.2
[ 1] ##contig=<ID=chr1,length=100>
[ 2] ##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
[ 3] ##FILTER=<ID=LowQual,Description="Low quality">          <-- added by this round
[ 4] ##FORMAT=<ID=AD,...>  [5] DP  [6] GQ  [7] GT  [8] PL  [9] RGQ
[10] ##INFO=<ID=AC,...>  [11] AF  [12] AN  [13] DP  [14] ExcessHet  [15] FS
[16] InbreedingCoeff  [17] MLEAC  [18] MLEAF  [19] MQ  [20] QD  [21] SOR  [22] AD
[23] ##source=GenotypeGVCFs
[24] #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	STAR
```

(Line 3 is where the input's own ordering would have put a `##FILTER` line:
immediately after `##ALT`/before the first `##FORMAT`; see §3.)

### 1.3 Position, precisely

* **Inside the group:** GATK's `##FILTER` lines are ordered by the lines' own
  text (htsjdk sorts a `VCFHeader`), measured with three spellings: `AAA` before
  `LowQual` before `q10`.  Native now reproduces exactly that order (§1.1).
* **Of the group in the header:** GATK puts it at index 2 (right after
  `##fileformat`/`##ALT`); native puts it at index 3 (right after its retained
  `##contig`, which htsjdk sorts to index 28).  **The absolute index therefore
  still differs by one.**  That is a property of the *pre-existing* whole-header
  ordering divergence (native preserves input order; GATK emits htsjdk's sorted
  header), not of the missing line — it is reported, not gated, not fixed
  (§7.2).

### 1.4 Other header lines: what differs (observations, NOT gated)

Measured for `covered-star-only-record-dense-non-pass`: **16** lines only in GATK
before the fix, **15** after (the removed one is the `##FILTER` line — see
§2.2), and **10** lines only in native in both runs.

| Kind | GATK | native | gated? |
| --- | --- | --- | --- |
| `##FILTER` | `LowQual` always | now the same line | **YES** (this round) |
| `##fileformat` | `VCFv4.2` | `VCFv4.2` | identical |
| `##contig` | sorted to position 28 | kept at position 1 | no (pre-existing) |
| `##GATKCommandLine` | written (the gate's GATK invocation does not pass `--add-output-vcf-command-line false`) | never written | no (pre-existing) |
| `##source` | `GenotypeGVCFs` | `GenotypeGVCFs` | identical |
| standard annotation lines | GATK adds `##INFO=BaseQRankSum`, `MQRankSum`, `ReadPosRankSum`, a second `##INFO=AD`, a second `##INFO=DP` ("Approximate read depth; some reads may have been filtered") and a second `##FORMAT=AD` ("Allelic depths for the ref and alt alleles in the order listed") | not added | no (pre-existing) |
| `Description=` quoting | htsjdk re-quotes every description (`Description="Allele depths"`) and rewrites MLEAC/MLEAF wording | native keeps the input's unquoted text (`Description=Allele depths`, `Description=Maximum likelihood allele count`) | no (pre-existing) |
| total | 31-32 lines | 25-26 lines | — |

So the `##FILTER` gap was the *only* header difference this round fixes, and no
new header difference was introduced (§2.2).  None of the rows marked "no" is
new, and none is in scope for this round.

### 1.5 Does the same gap exist for other tools?

* **HaplotypeCaller — no gap.** Native already emits the line verbatim, in both
  output profiles: `fastgatk-native/src/hc_call.cpp:3074` (plain VCF) and
  `:3795` (gVCF), both `##FILTER=<ID=LowQual,Description="Low quality">`, which
  matches GATK's `HaplotypeCallerEngine.java:559`.
* **GnarlyGenotyper — no native tool exists** (`ls fastgatk-native/src` has no
  gnarly implementation), so there is nothing to fix.
* Every other GATK tool that mentions `LOW_QUAL_FILTER_NAME` either *applies*
  the filter (`GenotypingEngine.java:185`, `GnarlyGenotyperEngine.java:106`,
  `GnarlyGenotyper.java:267`) or *declares* it (`GATKVCFHeaderLines.java:89`).
  No other native tool needs the line for this reason.

**Conclusion for step 1's last question: the gap was GenotypeGVCFs-only**, and
the fix is confined to that tool.

## 2. STEP 2 — the gate, written and run BEFORE the fix

**Extended** the already-registered
`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` rather
than adding a new script, for the same reason as the previous round: it already
runs pinned GATK and native on the **same fixture with the same arguments**, and
it is already registered in `fastgatk-native/CMakeLists.txt:1556-1560` under the
strict-gate filter of step 5(c).  `CMakeLists.txt` must not be edited this round,
so a brand-new file could not have been registered as a ctest gate at all — the
gate would have existed but never run in CI.

Added:

* a docstring section *"The output header's `##FILTER` lines"* quoting
  `GenotypeGVCFsEngine.java:416`, `:418-420`, `GenotypeGVCFs.java:305`,
  `GATKVCFHeaderLines.java:89`, `GATKVCFConstants.java:179` and the measured
  bytes/positions;
* `GATK_LOWQUAL_FILTER_LINE` (the measured declaration, byte for byte);
* `HEADER_WITH_Q10_FILTER` — the same fixture header with
  `##FILTER=<ID=q10,Description="Quality below 10">` declared instead of LowQual;
* `header_lines()` / `filter_header_lines()` helpers;
* in `run_case()`: capture both outputs' headers, record
  `gatk_filter_header_lines` / `native_filter_header_lines` and a
  **non-gating** `header_observations` block (`only_in_gatk`, `only_in_native`,
  `order_differs`; `##GATKCommandLine` excluded from the diff because it embeds
  the run's absolute paths and timestamp and native never writes it), and add
  the gated assertion

  ```python
  if gatk_filters != native_filters:
      result["violations"].append(
          "##FILTER header lines are not byte-identical: "
          f"GATK={gatk_filters} NATIVE={native_filters}")
  ```

  for **every** gated case (not just the new one), so the whole gate now covers
  the header's filter contract as well as the data rows;
* one new **GATED** case `header-filter-line-with-extra-declared-filter`
  (the `q10` fixture), which fails both when the `LowQual` declaration is
  missing *and* when it is appended at the wrong place in the group;
* the JSON payload now reports `header_compared` and `header_observations`.

### 2.1 Exit status and the literal header lines

**BEFORE the fix** — **exit 1**, `"status": "divergence"`, **22 violations**,
every one of them the `##FILTER` group and nothing else
(`.diag/header-filter-gate-before.log`); 23 gated cases, 22 of them with an
input header that declared no filter:

```
[default] VIOLATION: ##FILTER header lines are not byte-identical:
  GATK=['##FILTER=<ID=LowQual,Description="Low quality">'] NATIVE=[]
[include-non-variant-sites] ... GATK=['##FILTER=<ID=LowQual,Description="Low quality">'] NATIVE=[]
... (the same for all 22 cases with a filter-less input header) ...
[header-filter-line-with-extra-declared-filter] VIOLATION: ##FILTER header lines are not byte-identical:
  GATK=['##FILTER=<ID=LowQual,Description="Low quality">', '##FILTER=<ID=q10,Description="Quality below 10">']
  NATIVE=['##FILTER=<ID=q10,Description="Quality below 10">']
```

`non-pass-variant-input-declared-filter` (the one gated case whose input header
already declares `LowQual`) was **green before the fix** — the two halves of the
contract separate cleanly.  The data-row assertions stayed green throughout:
the 22 violations were header-only.

**AFTER the fix: exit 0**, `"status": "pass"`, 0 violations, 25 cases
(23 gated + 2 reported-only), every case's `##FILTER` group byte-identical to
GATK's, including the ordered pair for the `q10` case
(`.diag/header-filter-gate-after.log`).

### 2.2 The observations before vs after (proof the fix is exactly one line)

For `covered-star-only-record-dense-non-pass`:

```
only_in_gatk :  before 16  after 15     removed by fix: ['##FILTER=<ID=LowQual,Description="Low quality">']
                                        added   by fix: []
only_in_native: before 10  after 10     removed by fix: []   added by fix: []
```

i.e. the fix removed exactly the one divergent header line and introduced no
new header difference — every other header observation is unchanged and
pre-existing.

## 3. STEP 3 — the fix

`git diff --stat` = **2 files changed, 173 insertions(+), 1 deletion(-)**
(oracle +128, tool +45/-1).  The whole behaviour change in
`fastgatk-native/src/genotype_gvcf_tool.cpp`, inside
`gatk_compatible_header_text()` (the GATK-compatibility header builder used by
**both** writers: `:5602` streaming, `:6422` aggregate):

```diff
+constexpr const char* kGatkLowQualFilterLine =
+    "##FILTER=<ID=LowQual,Description=\"Low quality\">";
+
 std::string gatk_compatible_header_text(const std::string& formatted) {
 ...
     std::vector<std::string> output;
-    output.reserve(retained.size() + info_lines.size() + 1);
+    output.reserve(retained.size() + info_lines.size() + 2);
+    std::vector<std::string> filter_lines;
+    for (const auto& line : retained)
+        if (line.rfind("##FILTER=", 0) == 0) filter_lines.push_back(line);
+    if (std::find(filter_lines.begin(), filter_lines.end(), kGatkLowQualFilterLine) ==
+        filter_lines.end())
+        filter_lines.emplace_back(kGatkLowQualFilterLine);
+    std::sort(filter_lines.begin(), filter_lines.end());
     std::vector<std::string> format_lines;
 ...
     bool inserted_filter = false;
     bool inserted_format = false;
     for (const auto& line : retained) {
+        if (line.rfind("##FILTER=", 0) == 0) {
+            if (!inserted_filter) {
+                output.insert(output.end(), filter_lines.begin(), filter_lines.end());
+                inserted_filter = true;
+            }
+            continue;
+        }
         if (line.rfind("##FORMAT=", 0) == 0) { ... unchanged ... }
 ...
+    if (!inserted_filter) {
+        const auto anchor = std::find_if(output.begin(), output.end(), [](const auto& line) {
+            return line.rfind("##FORMAT=", 0) == 0 || line.rfind("#CHROM", 0) == 0;
+        });
+        output.insert(anchor, filter_lines.begin(), filter_lines.end());
+    }
```

Why this is the correct and minimal fix:

* **Text** is GATK's, byte for byte, quoted from
  `GATKVCFHeaderLines.java:89` + `GATKVCFConstants.java:179` and written into a
  named constant with the `file:line` provenance next to it.
* **Position inside the group** follows htsjdk's rule (the group is sorted by
  line text), which is what makes the ordered pair `LowQual, q10` come out
  right and `AAA, LowQual` come out right.
* **Position of the group** is inherited from the input (where the input's own
  `##FILTER` lines were); when the input declared none, the group is created
  where a `##FILTER` line would have sat in the input's own ordering — before
  the first `##FORMAT=` line.  This is the same collect/sort/re-insert shape the
  function already used for `##FORMAT`, so no new mechanism is introduced.
* **Unconditional is correct, and it was checked rather than assumed.**  The
  requirement was: if the line were added in a configuration where GATK does not
  declare it, that would be a new divergence.  Measured across default mode,
  `--include-non-variant-sites`, single- and two-sample inputs, the ordinary
  variant path, and four different input filter declarations, GATK declared it
  every time (§1.1), and the code path is unguarded with a single caller
  (`GenotypeGVCFsEngine.java:416`, `GenotypeGVCFs.java:305`).  So there is no
  configuration in which adding it unconditionally over-declares.
* **It adds nothing when the input already declares the line** (the
  `std::find` guard), which is why the fixtures whose gVCF comes from GATK
  HaplotypeCaller — and therefore already carries the line — are byte-unchanged
  in the header (§4, §6).
* It is local to the `--gatk-compatible-annotations` boundary: both call sites
  guard on `options.gatk_annotation_compatibility`, so the non-compatibility
  diagnostic profile (which has no `##FILTER` line and is explicitly not a GATK
  parity surface) is untouched, as are the Kokkos kernels, the merge/union
  stages and every other tool.

Not restructured: the genotyping engine, the FILTER-column handling of the
previous round, the AF/PL kernels and their ABI, the `##INFO`/`##FORMAT`
ordering logic, the GVCFBlock/PASS stripping, the `##source` insertion, and
every other tool.  `fastgatk-native/CMakeLists.txt` and every root `*.md` are
unmodified.

### 3.1 Collateral, explicitly

* The changed function is in the anonymous namespace of
  `genotype_gvcf_tool.cpp` and is called from exactly two places
  (`:5602`, `:6422`), both inside this tool's two encode stages.
* `output.reserve(...)` was bumped by one slot (a hint only).
* The `PASS` strip at `:4528` is *not* touched, and `filter_lines` is built from
  the already-stripped `retained`, so an input `##FILTER=<ID=PASS,...>` line is
  still dropped exactly as before (§7.3).

## 4. STEP 4 — stale assertions

**No registered test pinned the old behaviour; no test line needed correcting.**

| File / check | Why it could have been stale | Outcome |
| --- | --- | --- |
| `verify_gatk_genotype_gvcf.py:234-236` | asserts set equality between GATK's and native's *semantic* headers **including `##FILTER=`** | **not stale, and not affected**: its input is a GATK HaplotypeCaller gVCF, which already declares `LowQual`, so native retained it before this round and still does.  The `std::find` guard means the header is byte-identical to before for that input.  (This is also why the pre-existing set-equality assertion could never have caught this gap.) |
| `verify_gatk_genotype_gvcf_legacy_qual.py:105`, `verify_gatk_genotype_gvcf_multisample.py:116` | `assert native_header == gatk_header` — **order-sensitive** list equality | **not stale**: both build the list from `sorted({...})`, i.e. an order-insensitive set, and their inputs are HaplotypeCaller gVCFs that already declare `LowQual` |
| `verify_hc_dense_gvcf_genotype_gatk_oracle.py:128` | the strictest one: `normalized_joint_text(gatk) == normalized_joint_text(native)` compares **every** non-ignored header line in order | **not stale**: the input gVCF is native/GATK HaplotypeCaller output, whose header already has exactly one `##FILTER` line, and `--add-output-vcf-command-line false` removes GATK's command-line line — so the FILTER group is one element and the sort is a no-op.  Ran green in step 5(b2)/5(d) |
| `verify_gatk_genotype_gvcf_inbreeding.py:137-149`, `verify_gatk_genotype_gvcf_multiallelic.py:123-134` | compare header **ID sets** (`required <= ids`) | not stale: subset assertions, and adding GATK's own line cannot remove an id |
| `verify_genotype_gvcf*.py`, `verify_hc_dense_gvcf_genotype*`, `benchmark_genotype_gvcf.py` | could pin a header line *count* or a header *list* | `grep -n "len(header\|header_lines\|header_count"` over `fastgatk-native/scripts/*.py`: no genotype-output header count/list assertion |
| every `fastgatk-native/scripts/*.py` with `##FILTER` | could pin a header that lacks the line | 66 hits, all **fixture inputs** for combine/reblock/select-variants/variant-filtration/VQSR/variant-eval tools or already-correct expectations; none is a GenotypeGVCFs output-header expectation |
| `fastgatk-native/src/hc_call.cpp:3074`, `:3795` | already declared `LowQual`; the "other tool" half of step 1 | verified correct against `HaplotypeCallerEngine.java:559`; untouched |

The only test-script change in this round is the oracle itself (§2).

## 5. Diff summary

```
 fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py | 128 +++++++++++++++++++++
 fastgatk-native/src/genotype_gvcf_tool.cpp                          |  46 +++++++-
 2 files changed, 173 insertions(+), 1 deletion(-)
```

`git status --short` (nothing else in the tree):

```
 M fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py
 M fastgatk-native/src/genotype_gvcf_tool.cpp
```

## 6. STEP 5 — gate results (mandatory order)

| Step | Command | Result |
| --- | --- | --- |
| a | `python3 fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` | **exit 0** (25 cases: 23 gated + 2 reported-only, 0 violations); **exit 1 before** the fix (22 header violations, §2.1) |
| b1 | `python3 fastgatk-native/scripts/verify_genotype_gvcf.py` | `B1_EXIT=0` — `{"status": "pass", "output_records": 1}` (`.diag/header-filter-verify-gvcf.log`) |
| b2 | `ctest --test-dir fastgatk-native/build -R 'genotype-gvcf' -V` | `B2_EXIT=0` — **15/15 passed**, incl. `fastgatk-genotype-gvcf-spandel-gatk-oracle` (`.diag/header-filter-ctest-genotype-gvcf.log`) |
| c | the 17-name strict-gate filter from the task | `C_EXIT=0` — **19/19 matched tests passed**, incl. `fastgatk-genotype-gvcf-spandel-gatk-oracle` (212.14 s) (`.diag/header-filter-strict-gates.log`) |
| d | `fastgatk-native/scripts/run_regression.sh --label header-filter` | see below |

`fastgatk-native/scripts/run_regression.sh --label header-filter`
(omp + serial, `FASTGATK_REQUIRE_GATK_ORACLE=1` by default, ctest parallelism 8,
both backends in parallel):

```
REG_EXIT=0
| 后端   | 构建目录                              | 结果 | 通过/总数 | 耗时        |
| omp    | OpenMP (fastgatk-native/build)        | 通过 | 300/300   | 1284.07 sec |
| serial | Serial (fastgatk-native/build-serial) | 通过 | 300/300   | 1284.62 sec |
```

Evidence block `.diag/regression/20260912-072155/{omp.log,serial.log,evidence.md}`
(+ `.diag/header-filter-regression.out` for `REG_EXIT`), run window
07:21:55 → 07:43:20, label `header-filter`, git `d36c134` + the two uncommitted
files of §5.  The runner's staleness check reported no warning (`grep -c 警告`
= 0), i.e. no source file is newer than the binaries it tested.

**Nothing was edited after the suite run.**  `md5sum -c
.diag/header-filter-binaries.md5` re-confirms

```
09024bf17c6b82f1b3eb42d1c37b1e7e  fastgatk-native/build/fastgatk-genotype-gvcf        (omp)
b8a939ee9086610855910303e58a6155  fastgatk-native/build-serial/fastgatk-genotype-gvcf (serial)
```

and both binaries were built after the last source edit and before every gate
above (logs `.diag/header-filter-build-{omp,serial}.log`).  The two changed test
inputs are frozen at

```
deb31671db6ebf8917c67afa86416d52  fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py
eff2c633ffc720ec8ba9a8de95590c43  fastgatk-native/src/genotype_gvcf_tool.cpp
```

Timestamps, all
before the suite window 07:21:55 → 07:43:20: oracle
`verify_genotype_gvcf_spandel_gatk_oracle.py` 06:34:37 < tool
`genotype_gvcf_tool.cpp` 06:38:51 < omp binary 06:39:12 < serial binary
06:39:25.  The only file touched afterwards is this report, which is not a test
input (no source, test script, header or build input changed after the suite).
No rebuild and no re-run were therefore needed.

## 7. How far the divergence reached, and what is deliberately left

Measured, with the tree at the fixed revision:

| Surface | Before | After |
| --- | --- | --- |
| header, input with no `##FILTER` line, default mode | no `##FILTER` | `##FILTER=<ID=LowQual,Description="Low quality">` ✅ |
| header, input with no `##FILTER` line, `--include-non-variant-sites` | no `##FILTER` | same ✅ |
| header, ordinary variant path (non-dense) | no `##FILTER` | same ✅ |
| header, input declares an unrelated filter (`q10`) | `q10` only | `LowQual, q10` — GATK's order ✅ |
| header, input declares a filter sorting first (`AAA`) | `AAA` only | `AAA, LowQual` ✅ |
| header, input already declares `LowQual` | `LowQual` | `LowQual` (unchanged, byte-identical) ✅ |
| both writers (aggregate and `--stream-by-locus`) | shared function | both fixed ✅ |
| non-compatibility native diagnostic profile | no `##FILTER` | **no `##FILTER` — unchanged, by construction, NOT fixed** |
| input declaring `##FILTER=<ID=PASS,Description="All filters passed">` | native drops the line, GATK keeps it | **unchanged, NOT fixed** (§7.3) |
| absolute index of the group in the header | 2 (GATK) vs 3 (native) | **unchanged, NOT fixed** (§7.2) |

### 7.1 Non-compat profile (residual, by construction)

`gatk_compatible_header_text()` is only reached when
`options.gatk_annotation_compatibility` is set (`:5599`, `:6419`), so the native
diagnostic profile's header is byte-identical to before this round.  That profile
is explicitly not a GATK parity surface (it also publishes RCQ/RCP and skips the
AF-based allele pruning), so no change is wanted there.

### 7.2 Absolute position of the group (residual, reported)

native's compatibility header preserves the input's line order while GATK emits
htsjdk's sorted header, so the `##FILTER` group sits at index 3 rather than 2
(displaced by `##contig`, which htsjdk sorts to index 28).  Reproducing GATK's
absolute index would mean reproducing htsjdk's whole-header sort — including
`##contig`, `##source` and the second standard `##INFO=DP`/`##FORMAT=AD` lines
(§1.4) — i.e. rewriting the header builder.  *Speculation (not measured this
round):* that rewrite would also break
`verify_hc_dense_gvcf_genotype_gatk_oracle.py:128`, whose strict ordered header
comparison passes today only because native's order happens to coincide with
GATK's for HaplotypeCaller-derived inputs.
That is far outside "small and local", so it stays reported.  The gate asserts
the group's *content and internal order*, which is the part this round's
divergence was about; it deliberately does not assert the group's absolute index
(documented in the oracle docstring).

### 7.3 `##FILTER=<ID=PASS,...>` is stripped by native (residual, measured)

`genotype_gvcf_tool.cpp:4528` drops a `##FILTER=<ID=PASS,Description="All
filters passed">` line unconditionally
(`if (line == "##FILTER=<ID=PASS,...>") continue;`).  Measured: with such a line
in the input, GATK's output header keeps it (index 3, right after `LowQual`) and
native's does not.  This is a **pre-existing, separate** divergence — not the one
this round was dispatched for, and not gated (no fixture in the oracle declares
`PASS`, so the new gate never sees it).  It is left unfixed deliberately: the
line is an intentional native normalization, removing it changes behaviour for
every input that declares `PASS`, and the task bounds this round to the missing
`LowQual` declaration.  *Speculation:* deleting that one `continue` would likely
align the two, but its collateral is unmeasured.

## 8. What remains unproven

* The gate asserts the `##FILTER` group byte-for-byte on the fixtures it pins
  (23 gated cases; §1.1 covers 8 fixture/mode combinations and 5 distinct input
  filter spellings).  It is a point
  measurement, not a proof that native's filter-group handling matches GATK for
  every input shape.  Untested shapes: an input `##FILTER` line whose
  `Description` text differs from GATK's for the *same* id, `.bcf`/`.vcf.gz`
  inputs, `--keep-combined-raw-annotations`, `--dbsnp`, a GenomicsDB source, and
  inputs whose `##FILTER` lines are not contiguous in the header.
  *Speculation:* the code is a collect/sort/re-insert over the whole retained
  header, so line position within the input should not matter, but that is not
  measured.
* The ordering rule for the group is inferred as "htsjdk sorts the whole
  `VCFHeader` by key then value" from three measured spellings
  (`AAA`/`LowQual`/`q10`) and from the observed sorted order of the `##INFO` and
  `##FORMAT` groups in GATK's output.  I did not read htsjdk's comparator from
  the jar this round.
* §1.4's "other header differences are pre-existing" is established by the
  before/after observation diff (§2.2) for one fixture, not for all of them.
* §4's "no test pins the old behaviour" is a grep plus the step-5 results, not a
  proof: a test could construct such an expectation dynamically.  Step 5(d) is
  the empirical check (300/300 on both backends with
  `FASTGATK_REQUIRE_GATK_ORACLE=1`), which is evidence, not proof.
* The residuals in §7.2 and §7.3 are measured but not repaired; whether they can
  be repaired without disturbing the HaplotypeCaller-derived oracles is
  **unproven**.

## 9. Tree state

* Change **kept in the tree** (no commit, no branch, no branch switch).
* Modified — exactly two files (§5).  Scratch lives under `.diag/` only, which
  is git-ignored.
* Not touched: `fastgatk-native/CMakeLists.txt`, every root `*.md`, Mutect2,
  every other tool, and every registered test script other than the oracle.
* Artifacts: probe `.diag/header_filter_probe.py` (+ logs/JSON
  `.diag/header-filter-probe*.log|json`), gate logs
  `.diag/header-filter-gate-{before,after}.log`, build logs
  `.diag/header-filter-build-{omp,serial}.log`, step-5 logs
  `.diag/header-filter-{verify-gvcf,ctest-genotype-gvcf,strict-gates}.log`,
  binary hashes `.diag/header-filter-binaries.md5`.
