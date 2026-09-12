# Round: the input's `##FILTER=<ID=PASS,...>` line (GATK propagates it, native stripped it)

Scope: the residual divergence the previous round measured and deliberately left
unfixed (`.diag/round-header-filter.md` §7.3): on an input whose header declares
`##FILTER=<ID=PASS,Description="All filters passed">`, pinned GATK 4.6.2.0
GenotypeGVCFs keeps the line in its output header while native dropped it at
`fastgatk-native/src/genotype_gvcf_tool.cpp:4528`.  Gate first, fix second.

## 0. Verdict

* **GATK truth (measured, pinned 4.6.2.0):** GenotypeGVCFs writes back the
  input header's own `##FILTER` lines **verbatim**, the reserved `PASS` line
  included, with the input's own `Description` text, and it never synthesizes a
  `PASS` line of its own.  Source:
  `final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());`
  (`gatk-source/.../GenotypeGVCFsEngine.java:395`), the only removal being
  `headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().startsWith(GVCF_BLOCK));`
  (`:398-399`); the only `##FILTER` line the tool adds is LowQual (`:416`).  The
  group is emitted in htsjdk's sorted order, so LowQual precedes PASS precedes
  `q10`.
* **The retention is a GATK choice, not an htsjdk behaviour.**  Neither the
  pinned `gatk-package-4.6.2.0-local.jar` nor `htsjdk-4.2.0.jar` contains the
  byte string `All filters passed`; htsjdk's `VCFHeader` never references
  `VCFConstants.PASSES_FILTERS_v4`; and GATK's own PASS line
  (`GATKVCFHeaderLines.java:90`, *"Site contains at least one allele that passes
  filters"*) is used only by VQSR.  Decisive control: an input with **no**
  `##FILTER` line at all produces **no** PASS line in GATK's output.  So the
  text can only come from the input header.
* **The special-casing is on the native side, and it is HTSlib's.**
  `bcf_hdr_parse()` injects a synthetic
  `##FILTER=<ID=PASS,Description="All filters passed">` record as the **second**
  record of *every* header it parses — *"The filter PASS must appear first in the
  dictionary"* (`third_party/htslib-build/htslib-src/vcf.c:1271-1276`) — and
  `bcf_hdr_add_hrec()` then destroys any later `##FILTER` record whose ID is
  already registered (`vcf.c:1092-1111`, the destroy at `:1107-1110`).  So the
  parsed header carries HTSlib's canonical text whether or not the input
  declared PASS, and **loses the input's own text** when it did.  Reproduced
  with the shipped binary: `.diag/filter-pass-htslib-injection.log` shows the raw
  HTSlib header (`bcf_hdr_write`, non-compatibility profile) of a PASS-less input
  containing that line, while the GATK-compatible header of the same input
  contains only LowQual.
* Consequently the native strip at `:4528` was **half right**: it removed
  HTSlib's synthetic line, but it also removed a genuine input declaration.
  Simply deleting it is **wrong** — measured, native then emits PASS for inputs
  that never declared it (§3.1).
* **Fix kept in the tree:** one local change in
  `fastgatk-native/src/genotype_gvcf_tool.cpp` — a new
  `input_pass_filter_line()` helper that reads the input file's own
  `##FILTER=<ID=PASS...>` line back out of the raw input text (HTSlib has already
  destroyed it), plus one `continue` retargeted at HTSlib's synthetic line and a
  re-insertion of the input's line into the `##FILTER` group the previous round
  built.  `+82/-4` in that file; no restructuring, no kernel ABI, no
  `CMakeLists.txt`, no Mutect2, no other tool.
* Oracle: **exit 1 before** the fix (3 `##FILTER`-only violations, §2.1),
  **exit 0 after** (28 cases: 26 gated + 2 reported-only, 0 violations).
* Step 5: (a) exit 0; (b1) exit 0, (b2) 15/15 exit 0; (c) 19/19 exit 0;
  (d) **300/300 on both backends**, `REG_EXIT=0`.  Change **kept in the tree**.

## 1. STEP 1 — measured truth, byte for byte

Reproducer: `.diag/filter_pass_probe.py` — log/JSON
`.diag/filter-pass-probe-before.{log,json}` (before) and
`.diag/filter-pass-probe-after.{log,json}` (after, final fix).  It builds each
fixture on the oracle's 100 bp `chr1` reference, indexes it with
`gatk IndexFeatureFile`, and runs

```
GATK:   java -jar gatk-package-4.6.2.0-local.jar GenotypeGVCFs -R ref.fa -V in.g.vcf [-all-sites] \
          -O out.vcf --create-output-variant-index false
native: fastgatk-genotype-gvcf -R ref.fa -V in.g.vcf [-all-sites] \
          --gatk-compatible-annotations -O out.vcf
```

dumping every `##FILTER` line of both outputs with its whole-header index.

### 1.1 The `##FILTER` lines, both sides, every fixture

| # | input header declares | mode | GATK `##FILTER` (index) | native BEFORE (index) | native AFTER (index) |
| --- | --- | --- | --- | --- | --- |
| 1 | *(no `##FILTER` line at all)* | default | `LowQual` (2) | `LowQual` (3) | `LowQual` (3) |
| 2 | `PASS` / `"All filters passed"` | default | `LowQual` (2), `PASS` (3) | `LowQual` (3) — **PASS dropped** | `LowQual` (3), `PASS` (4) |
| 3 | `PASS` / `"All filters passed"` | `-all-sites` | `LowQual` (2), `PASS` (3) | `LowQual` (3) — **PASS dropped** | `LowQual` (3), `PASS` (4) |
| 4 | `PASS` / `"Some other PASS text"` | default | `LowQual` (2), `PASS`*(that text)* (3) | `LowQual` (3) — **PASS dropped** | `LowQual` (3), `PASS`*(that text)* (4) |
| 5 | `PASS` / `"All filters passed"`, then `q10` | default | `LowQual` (2), `PASS` (3), `q10` (4) | `LowQual` (3), `q10` (4) — **PASS dropped** | `LowQual` (3), `PASS` (4), `q10` (5) |
| 6 | `q10` then `PASS` (reversed input order) | default | `LowQual` (2), `PASS` (3), `q10` (4) | `LowQual` (3), `q10` (4) — **PASS dropped** | `LowQual` (3), `PASS` (4), `q10` (5) |
| 7 | `LowQual` only | default | `LowQual` (2) | `LowQual` (3) | `LowQual` (3) |
| 8 | `q10` only | default | `LowQual` (2), `q10` (3) | `LowQual` (3), `q10` (4) | `LowQual` (3), `q10` (4) |
| 9 | `AAA`, `LowQual`, `q10`, `PASS` | default | `AAA` (2), `LowQual` (3), `PASS` (4), `q10` (5) | `AAA` (3), `LowQual` (4), `q10` (5) — **PASS dropped** | `AAA` (3), `LowQual` (4), `PASS` (5), `q10` (6) |
| 10 | `PASS`, then `StrandBias` | default | `LowQual` (2), `PASS` (3), `StrandBias` (4) | `LowQual` (3), `StrandBias` (4) — **PASS dropped** | `LowQual` (3), `PASS` (4), `StrandBias` (5) |

The literal lines (fixture 2), both sides after the fix:

```
GATK  : ['##FILTER=<ID=LowQual,Description="Low quality">',
         '##FILTER=<ID=PASS,Description="All filters passed">']        # indices 2,3
NATIVE: ['##FILTER=<ID=LowQual,Description="Low quality">',
         '##FILTER=<ID=PASS,Description="All filters passed">']        # indices 3,4
```

Fixture 4 confirms the verbatim-propagation rule and refutes "canonical PASS is
re-emitted":

```
GATK  : [..., '##FILTER=<ID=PASS,Description="Some other PASS text">']
NATIVE (after): [..., '##FILTER=<ID=PASS,Description="Some other PASS text">']
```

Fixture 9 shows the ordering rule is a plain sort of the whole line text:
`AAA` < `LowQual` < `PASS` < `q10`, independent of the input's order
(fixture 5 vs 6 differ only in the input's order and produce the identical GATK
group).

### 1.2 The absolute index of the group is still off by one (unchanged, not gated)

GATK puts the group at index 2 (right after `##fileformat`/`##ALT`); native puts
it at index 3, displaced by its retained `##contig`, which htsjdk sorts to
position 28.  This is the pre-existing whole-header ordering divergence recorded
in `.diag/round-header-filter.md` §7.2 — untouched here, and the gate asserts the
group's **content and internal order**, not its absolute index.

### 1.3 Does the same gap exist for other tools?

Not in scope and not measured this round.  The only other native tool that
mentions the canonical PASS text is none — `grep -rn 'All filters passed'
fastgatk-native/src/ fastgatk-core/` returns exactly one hit, the strip that this
round changed.  `##FILTER=<ID=PASS,...>` fixtures used by
`verify_combine_gvcfs_*.py` / `verify_reblock_*.py` are inputs of *different*
tools, which were not touched.

## 2. STEP 2 — the gate, extended and run BEFORE the fix

**Extended** the already-registered
`fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` rather than
adding a new script, for the same reason as the previous round: it already runs
pinned GATK and native on the same fixture with the same arguments, and it is
already registered in `fastgatk-native/CMakeLists.txt:1556-1560` under the
strict-gate filter of step 5(c).  `CMakeLists.txt` must not be edited this round,
so a brand-new file could not have been registered as a ctest gate at all.

Added (oracle diff `+120/-0`):

* a docstring section *"`##FILTER=<ID=PASS,...>` is propagated, not generated"*
  quoting `GenotypeGVCFsEngine.java:395`, `:398-399`, `:416`,
  `GATKVCFHeaderLines.java:90`, the measured fixture matrix, and the negative
  control (no input `##FILTER` line ⇒ no PASS line in GATK's output);
* `HEADER_WITH_PASS_FILTER`, `HEADER_WITH_NONSTANDARD_PASS_FILTER`,
  `HEADER_WITH_PASS_AND_Q10_FILTER` (the same fixture header with the respective
  `##FILTER` declarations inserted where a filter line would sit);
* three new **GATED** cases, each asserting the whole `##FILTER` group
  byte-for-byte **in order** through the assertion the previous round added
  (`if gatk_filters != native_filters: ... violations.append(...)`), so each
  fails both when a line is missing and when it is ordered differently:
  * `header-filter-line-with-input-pass-declaration`,
  * `header-filter-line-with-nonstandard-pass-description`,
  * `header-filter-line-with-pass-and-other-filter`.

### 2.1 Exit status and the literal lines from both sides

**BEFORE the fix** — **exit 1**, `"status": "divergence"`, 28 cases,
**3 violations**, every one of them the `##FILTER` group of a new case and
nothing else (`.diag/filter-pass-gate-before.log`); the pre-existing 25 cases
stayed green:

```
[header-filter-line-with-input-pass-declaration]
  ##FILTER header lines are not byte-identical:
  GATK=['##FILTER=<ID=LowQual,Description="Low quality">', '##FILTER=<ID=PASS,Description="All filters passed">']
  NATIVE=['##FILTER=<ID=LowQual,Description="Low quality">']

[header-filter-line-with-nonstandard-pass-description]
  ##FILTER header lines are not byte-identical:
  GATK=['##FILTER=<ID=LowQual,Description="Low quality">', '##FILTER=<ID=PASS,Description="Some other PASS text">']
  NATIVE=['##FILTER=<ID=LowQual,Description="Low quality">']

[header-filter-line-with-pass-and-other-filter]
  ##FILTER header lines are not byte-identical:
  GATK=['##FILTER=<ID=LowQual,Description="Low quality">', '##FILTER=<ID=PASS,Description="All filters passed">', '##FILTER=<ID=q10,Description="Quality below 10">']
  NATIVE=['##FILTER=<ID=LowQual,Description="Low quality">', '##FILTER=<ID=q10,Description="Quality below 10">']
```

The data rows of all three new cases were **already byte-identical** before the
fix (`chr1 2 . AA A 92.60 . AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63
GT:AD:DP:GQ:PL 0/1:0,20:20:99:100,0,100`), so the only divergence was the header.

**AFTER the fix: exit 0**, `"status": "pass"`, 0 violations, 28 cases
(26 gated + 2 reported-only), every case's `##FILTER` group byte-identical to
GATK's, including the three-element ordered group of
`header-filter-line-with-pass-and-other-filter`
(`.diag/filter-pass-gate-after.log`).

## 3. STEP 3 — the fix

`git diff --stat` = **2 files changed, 202 insertions(+), 4 deletions(-)**
(oracle `+120/-0`, tool `+82/-4`).

### 3.1 The first attempt, and why it was discarded (measured, NOT in the tree)

The obvious minimal edit — delete the `continue` at `:4528` — was built and
measured first.  It is **wrong**: because HTSlib injects the synthetic PASS line
into every header it parses, native then emitted

```
NATIVE (naive attempt) ['##FILTER=<ID=LowQual,Description="Low quality">',
                        '##FILTER=<ID=PASS,Description="All filters passed">']
```

for the `no-filter-line` fixture, i.e. it *introduced* a new divergence on every
input that declares no PASS line, and it still could not reproduce fixture 4
(it emitted HTSlib's canonical text instead of the input's
`Description="Some other PASS text"`).  It also moved the group to header index 1
(HTSlib's synthetic record sits right after `##fileformat`).  That build was
replaced by the one below; the archived evidence of the underlying cause is
`.diag/filter-pass-htslib-injection.log`, produced with the **final** binary.
(The probe logs of that intermediate build were overwritten by the re-run, so
the naive build's numbers above are quoted from the session, not from a log.)

### 3.2 The change that is in the tree

`fastgatk-native/src/genotype_gvcf_tool.cpp`, all inside the anonymous-namespace
header construction used by **both** writers (streaming `:5679`, aggregate
`:6502`):

```diff
+constexpr const char* kHtslibSyntheticPassFilterLine =
+    "##FILTER=<ID=PASS,Description=\"All filters passed\">";
+
+// (comment quoting GenotypeGVCFsEngine.java:395/:398-399/:416,
+//  GATKVCFHeaderLines.java:90 and htslib vcf.c:1271-1276 / :1092-1111)
+std::string input_pass_filter_line(const std::string& path) {
+    if (path.empty() || path == "-") return {};
+    htsFile* input = hts_open(path.c_str(), "r");
+    if (input == nullptr) return {};
+    std::string found;
+    if (input->format.format == vcf) {              // BCF has no text header
+        kstring_t line{0, 0, nullptr};
+        while (hts_getline(input, '\n', &line) >= 0) {
+            if (line.l == 0) continue;
+            const std::string text(line.s, line.l);
+            if (text.rfind("#CHROM", 0) == 0) break;
+            if (text.rfind("##FILTER=<ID=PASS,", 0) == 0 ||
+                text.rfind("##FILTER=<ID=PASS>", 0) == 0) { found = text; break; }
+        }
+        free(line.s);
+    }
+    hts_close(input);
+    return found;
+}
+
-std::string gatk_compatible_header_text(const std::string& formatted) {
+std::string gatk_compatible_header_text(const std::string& formatted,
+                                        const std::string& input_pass_filter) {
 ...
-        if (line == "##FILTER=<ID=PASS,Description=\"All filters passed\">") continue;
+        // HTSlib's own synthetic PASS declaration, not the input's ...
+        if (line == kHtslibSyntheticPassFilterLine) continue;
 ...
     std::vector<std::string> filter_lines;
     for (const auto& line : retained)
         if (line.rfind("##FILTER=", 0) == 0) filter_lines.push_back(line);
+    // The input's own PASS declaration, when the input file really has one.
+    if (!input_pass_filter.empty() &&
+        std::find(filter_lines.begin(), filter_lines.end(), input_pass_filter) ==
+            filter_lines.end())
+        filter_lines.emplace_back(input_pass_filter);
     if (std::find(filter_lines.begin(), filter_lines.end(), kGatkLowQualFilterLine) ==
         filter_lines.end())
         filter_lines.emplace_back(kGatkLowQualFilterLine);
```

and at both call sites the new argument:

```diff
                 write_vcf_text_line(output, gatk_compatible_header_text(
                     std::string(formatted_header.s == nullptr ? "" : formatted_header.s,
-                                formatted_header.l))) != 0) {
+                                formatted_header.l),
+                    input_paths.empty()
+                        ? std::string{}
+                        : input_pass_filter_line(input_paths.front()))) != 0) {
```

Why this is correct and minimal:

* **The strip stays, retargeted at the real owner.**  The line being skipped is
  HTSlib's synthetic record, which is present in the formatted header for every
  input; the input's genuine declaration is destroyed by HTSlib before native
  ever sees it, so it has to be read from the file.
* **Text and order are GATK's**, byte for byte: the input's own line goes into
  the `##FILTER` group the previous round built, which sorts the group by line
  text and pins LowQual's presence — exactly the mechanism that already
  reproduced `AAA, LowQual` and `LowQual, q10`.  No new ordering logic.
* **The LowQual interaction is handled by construction**: the input's PASS line
  is added *before* LowQual is ensured and *before* the group is sorted, so
  `LowQual` precedes `PASS` (measured GATK order) without a special case.
* **No over-declaration**: the raw read is what distinguishes "the input declared
  PASS" from "HTSlib injected it", so fixture 1 (and every existing gated case,
  whose inputs declare no PASS) is unchanged — verified by the gate staying green
  for the 25 pre-existing cases (§2.1) and by the full suite (§5).
* **Best-effort by design**: an empty, stdin (`-`), BCF, unreadable or
  GenomicsDB-expanded path yields `""`, which falls back to exactly the previous
  behaviour (drop the synthetic line).  No new failure mode is introduced — the
  helper cannot throw.
* Verified to work over **gzip/BGZF** inputs as well as plain VCF (raw read via
  `hts_open`/`hts_getline`, `input->format.format == vcf` guard):
  `.gz` fixtures with no PASS / PASS / non-standard PASS / `q10` produced
  `[LowQual]`, `[LowQual,PASS]`, `[LowQual,PASS(odd)]`, `[LowQual,q10]`
  respectively.

Not restructured: the genotyping engine, the FILTER-column handling, the AF/PL
kernels and their ABI, the `##INFO`/`##FORMAT` ordering logic, the GVCFBlock
stripping, the `##source` insertion, the non-compatibility diagnostic profile,
and every other tool.  `fastgatk-native/CMakeLists.txt` and every root `*.md` are
unmodified.

## 4. STEP 4 — stale assertions

**No registered test pinned the stripped behaviour; no test line needed
correcting.**  The lines checked, and why each was or was not stale:

| File / check | Why it could have been stale | Outcome |
| --- | --- | --- |
| `verify_genotype_gvcf_spandel_gatk_oracle.py:1234-1237` (the assertion added last round) | compares the `##FILTER` group byte-for-byte in order | **not stale — extended**: it is the one that caught this divergence (§2) and now pins it |
| `verify_gatk_genotype_gvcf.py:234-236` | set equality of *semantic* headers including `##FILTER=` | **not stale**: input is a HaplotypeCaller gVCF, which declares no PASS (`HaplotypeCallerEngine.java:557-559` adds only LowQual), so nothing changes; the comparison is order-insensitive |
| `verify_gatk_genotype_gvcf_legacy_qual.py:105`, `verify_gatk_genotype_gvcf_multisample.py:116` | `assert native_header == gatk_header`, order-sensitive | **not stale**: the lists are built from `sorted({...})`; same HaplotypeCaller-derived input |
| `verify_hc_dense_gvcf_genotype_gatk_oracle.py:128` | `normalized_joint_text(gatk) == normalized_joint_text(native)` compares every non-ignored header line **in order** | **not stale**: the input gVCF is produced by HaplotypeCaller at run time and declares no PASS; ran green in step 5(b2)/(c) |
| `verify_genotype_gvcf_assignment_gatk_oracle.py:27-40,118,134` | parses `##FILTER=<` header lines | **not stale**: `header_ids()` returns a *set* of `(kind, id)` pairs and the assertions are subset tests (`required_headers <= ...`), which an extra line cannot break |
| `verify_gatk_genotype_gvcf_inbreeding.py:137-149`, `verify_gatk_genotype_gvcf_multiallelic.py:123-134` | header **ID sets** (`required <= ids`) | **not stale**: subset assertions |
| `grep -rln 'All filters passed' fastgatk-native/scripts/` | could pin a fixture-header expectation | 4 hits, all `verify_combine_gvcfs_*.py` / `verify_reblock_*.py` — fixtures of *other* tools, unmodified |
| `fastgatk-native/src/hc_call.cpp:3074`, `:3795` | other tool's header | untouched; HaplotypeCaller declares only LowQual, matching `HaplotypeCallerEngine.java:559` |

The only test-script change in this round is the oracle itself (§2).

## 5. STEP 5 — gate results (mandatory order)

| Step | Command | Result |
| --- | --- | --- |
| a | `python3 fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py` | **exit 0** (28 cases: 26 gated + 2 reported-only, 0 violations); **exit 1 before** the fix (3 header violations, §2.1) |
| b1 | `python3 fastgatk-native/scripts/verify_genotype_gvcf.py` | `B1_EXIT=0` — `{"status": "pass", "output_records": 1}` (`.diag/filter-pass-verify-gvcf.log`) |
| b2 | `ctest --test-dir fastgatk-native/build -R 'genotype-gvcf' -V` | `B2_EXIT=0` — **15/15 passed**, incl. `fastgatk-genotype-gvcf-spandel-gatk-oracle` (`.diag/filter-pass-ctest-genotype-gvcf.log`) |
| c | the 17-name strict-gate filter from the task | `C_EXIT=0` — **19/19 matched tests passed**, incl. `fastgatk-genotype-gvcf-spandel-gatk-oracle` (`.diag/filter-pass-strict-gates.log`) |
| d | `fastgatk-native/scripts/run_regression.sh --label filter-pass-line` | see below |

`fastgatk-native/scripts/run_regression.sh --label filter-pass-line`
(omp + serial, `FASTGATK_REQUIRE_GATK_ORACLE=1` by default, ctest parallelism 8,
both backends in parallel):

```
REG_EXIT=0
| 后端   | 构建目录                              | 结果 | 通过/总数 | 耗时        |
| omp    | OpenMP (fastgatk-native/build)        | 通过 | 300/300   | 1294.02 sec |
| serial | Serial (fastgatk-native/build-serial) | 通过 | 300/300   | 1286.24 sec |
```

Evidence block `.diag/regression/20260912-084612/{omp.log,serial.log,evidence.md}`
(+ `.diag/filter-pass-regression.out` for `REG_EXIT`), run window
08:46:12 → 09:07:46, label `filter-pass-line`, git `e5acb66` + the two uncommitted
files of §7.  The runner's staleness check reported no warning
(`grep -c 警告` = 0), i.e. no source file is newer than the binaries it tested.

**Nothing was edited after the suite run.**  Provenance:

```
07:51:47  fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py
07:59:54  fastgatk-native/src/genotype_gvcf_tool.cpp
08:00:07  fastgatk-native/build/fastgatk-genotype-gvcf          (omp)
08:06:45  fastgatk-native/build-serial/fastgatk-genotype-gvcf   (serial)
08:46:12 → 09:07:46   regression window
```

`md5sum -c .diag/filter-pass-binaries.md5` re-confirms

```
f1f2ede3ce31e3d64d8afc6867618133  fastgatk-native/build/fastgatk-genotype-gvcf        (omp)
9ce79361729c804870e1eb1cd260a9d6  fastgatk-native/build-serial/fastgatk-genotype-gvcf (serial)
```

and both were built after the last source edit and before every gate above
(`.diag/filter-pass-build-{omp,serial}.log`).  The two changed test inputs are
frozen at

```
14b6f8f2ab8897353626590eed2da21e  fastgatk-native/src/genotype_gvcf_tool.cpp
81cbbb450736ee80ec82e03c22e94d36  fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py
```

The only file touched at/after the suite window is **this report** (written
08:57:18, i.e. during the run) plus the read-only evidence probe
`.diag/filter_pass_htslib_injection_probe.py` / `.diag/filter-pass-htslib-injection.log`
(run 08:5x against the same frozen binary).  Neither is a test input — no source,
test script, header or build input changed after the binaries were built, so no
rebuild and no re-run were required.

## 6. What remains unproven

* The gate pins the `##FILTER` group for the fixtures it runs (28 cases, 3 of
  them new).  §1.1 covers 10 fixture/mode combinations and 5 distinct input
  filter spellings, including the two PASS spellings.  It is a point
  measurement, not a proof for every input shape.
* **Multi-input runs are not gated.**  The raw read uses
  `input_paths.front()`, the same source that already provides every other
  retained header line (`output_header = bcf_hdr_dup(header)` of the first
  input, `:5483`), whereas GATK builds its output header from the **merged**
  header of all `-V` sources (`GenotypeGVCFsEngine.java:395`).  So if shard 1
  declares no PASS but shard 2 does, GATK emits PASS and native still will not.
  *Speculation:* the same applies to every other filter line and is therefore a
  pre-existing property of native's first-input header, not something this
  change introduces; it is not measured.
* **BCF input is unhandled by the new helper** (`input->format.format == vcf`
  guard), so a BCF input that declares PASS keeps today's behaviour (line
  dropped).  Not measured.
* **GenomicsDB / `gendb://` inputs**: the helper is best-effort and falls back
  silently, so if the expanded path is not a text VCF the previous behaviour is
  kept.  Not measured.
* **The naive fix was measured but its logs were overwritten** (§3.1); the
  numbers quoted there come from the session, not from an archived log.  The
  underlying cause is archived independently in
  `.diag/filter-pass-htslib-injection.log` using the final binary.
* The claim "no registered test pinned the stripped behaviour" (§4) is a grep
  plus the step-5 results, not a proof: a test could construct such an
  expectation dynamically.  Step 5(d) is the empirical check (300/300 on both
  backends with `FASTGATK_REQUIRE_GATK_ORACLE=1`), which is evidence, not proof.

## 7. Tree state

* Change **kept in the tree** (no commit, no branch, no branch switch).
* Modified — exactly two files:
  `fastgatk-native/src/genotype_gvcf_tool.cpp` (`+82/-4`) and
  `fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py`
  (`+120/-0`).
* NOT touched: `fastgatk-native/CMakeLists.txt`, every root `*.md`, Mutect2,
  every other tool, and every registered test script other than the oracle.
* Artifacts: probe `.diag/filter_pass_probe.py`
  (+ `.diag/filter-pass-probe-{before,after}.{log,json}`), injection probe
  `.diag/filter_pass_htslib_injection_probe.py`
  (+ `.diag/filter-pass-htslib-injection.log`), gate logs
  `.diag/filter-pass-gate-{before,after}.log`, build logs
  `.diag/filter-pass-build-{omp,serial}.log`, step-5 logs
  `.diag/filter-pass-{verify-gvcf,ctest-genotype-gvcf,strict-gates,regression}.log|out`,
  regression evidence `.diag/regression/20260912-084612/`, hashes
  `.diag/filter-pass-binaries.md5` and `.diag/filter-pass-testfiles.md5`.
