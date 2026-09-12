# Round: the output header's ORDER (htsjdk sorts the whole header; native did not)

Scope: the last header class left open by the two previous rounds — GATK's output
header is written by htsjdk **in sorted order**, native's compatibility header kept
the input's line order plus its own INFO rank list.  Measured consequences on the
gate fixture: `##contig` at position 27 (GATK) vs 1 (native), `INFO/AD` right after
`AC` (GATK) vs last (native), and the two `##INFO=<ID=DP,...>` lines swapped; plus
the `##FILTER` group's absolute index (GATK 2 vs native 3).  Gate first, fix second.

## 0. Verdict

* **htsjdk's comparator, read out of the pinned jar (not inferred):** the writer
  writes `##fileformat` first and then every other metadata line in
  `VCFHeaderLine.compareTo` = **`toString().compareTo(...)`, i.e. the full line text
  without the leading `##`** — with one override: `VCFContigHeaderLine.compareTo`
  orders two contig lines by their **`contigIndex`** (the order the *input* header
  declared them) and falls back to text against anything else.  `#CHROM` is written
  last.  §1 has the four `javap -c` excerpts and the measurement.
* **Confirmed on four measured fixtures with four different input orders**
  (contig-first, already-text-sorted, reversed, and dense mode): GATK's own output
  satisfies that comparator in every one, and native diverged at 14 / 5 / 14 / 14
  positions.
* **Blast radius: one order-sensitive registered test**, and it is order-sensitive
  only for the lines it does not ignore.  The change is confined to output
  serialization of one tool's GATK-compatibility header: no data-row path, no
  kernel ABI, no other tool, no non-compat profile.  §2 names the test and the
  18 scripts checked.
* **Fix kept in the tree:** one local change in
  `fastgatk-native/src/genotype_gvcf_tool.cpp` (`+87/-0`): a new
  `gatk_htsjdk_sorted_header()` helper that re-sorts the built header exactly the
  way htsjdk's writer sorts it, called once at the end of
  `gatk_compatible_header_text()`.
* **Gate:** new
  `fastgatk-native/scripts/verify_genotype_gvcf_header_order_gatk_oracle.py`
  (4 fixtures, 1 violation each).  **exit 1 before** the fix (§3.2, literal
  position-by-position diff), **exit 0 after** (0 mismatched positions in all 4).
* Step 5: (a) exit 0; (b1) exit 0; (b2) 16/16 exit 0; (c) 21/21 exit 0;
  (d) **302/302 on BOTH backends**, `REG_EXIT=0`.  Change **kept in the tree**.

## 1. STEP 1 — htsjdk's order, from the jar and from measurement

### 1.1 The four `javap -c` excerpts (pinned jar, bytecode verbatim)

`third_party/jdk17/bin/javap -p -c -classpath
third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar <class>`

```
htsjdk.variant.vcf.VCFHeaderLine  (base class of every line type except contigs)
  public int compareTo(java.lang.Object);
     0: aload_0
     1: invokevirtual #115   // Method toString:()Ljava/lang/String;
     4: aload_1
     5: invokevirtual #116   // Method java/lang/Object.toString:()Ljava/lang/String;
     8: invokevirtual #119   // Method java/lang/String.compareTo:(Ljava/lang/String;)I
    11: ireturn

htsjdk.variant.vcf.VCFContigHeaderLine  (the ONE override)
  public int compareTo(java.lang.Object);
     0: aload_1
     1: instanceof    #2     // class htsjdk/variant/vcf/VCFContigHeaderLine
     4: ifeq          22
     7: aload_0
     8: getfield      #39    // Field contigIndex:Ljava/lang/Integer;
    11: aload_1
    12: checkcast     #2
    15: getfield      #39    // Field contigIndex:Ljava/lang/Integer;
    18: invokevirtual #130   // Method java/lang/Integer.compareTo:(Ljava/lang/Integer;)I
    21: ireturn
    22: aload_0
    23: aload_1
    24: invokespecial #132   // Method htsjdk/variant/vcf/VCFSimpleHeaderLine.compareTo:(Ljava/lang/Object;)I
    27: ireturn

htsjdk.variant.vcf.VCFHeader
  public java.util.Set<VCFHeaderLine> getMetaDataInSortedOrder();
     0: aload_0
     1: new           #461   // class java/util/TreeSet
     4: dup
     5: aload_0
     6: getfield      #103   // Field mMetaData:Ljava/util/Set;
     9: invokespecial #462   // Method java/util/TreeSet."<init>":(Ljava/util/Collection;)V
    12: invokevirtual #458   // Method makeGetMetaDataSet:(Ljava/util/Set;)Ljava/util/Set;
    15: areturn

htsjdk.variant.variantcontext.writer.VCFWriter
  public static VCFHeader writeHeader(VCFHeader, Writer, String versionString, String streamName);
     // ... writer.write(<"##fileformat=" + versionString + "\n">)
    14: invokevirtual #181   // Method VCFHeader.getMetaDataInSortedOrder:()Ljava/util/Set;
    18: invokeinterface #187, 1 // InterfaceMethod java/util/Set.iterator:()Ljava/util/Iterator;
    ...
    47: aload 5
    49: invokevirtual #202   // Method VCFHeaderLine.getKey:()Ljava/lang/String;
    52: invokestatic  #208   // Method VCFHeaderVersion.isFormatString:(Ljava/lang/String;)Z
    55: ifeq 61
    58: goto 25                       // skip the ##fileformat line in the loop
    61: ... writer.write("##" + line + "\n")
    85: ... writer.write("#" + HEADER_FIELDS + samples + "\n")
```

### 1.2 The comparator, stated precisely

1. `##fileformat` is **not** part of the sort: the writer emits it first from the
   header version, and skips any `isFormatString(key)` line inside the loop
   (`makeGetMetaDataSet` injects exactly such a line, which is why the skip exists).
2. Every other line is ordered by `key + "=" + value`, i.e. **by key first, then by
   the rest of the line**, compared as Java `String.compareTo` (UTF-16 code units;
   identical to `std::string` byte order for every ASCII header line measured).
   There is no ID/type sub-ordering: it is the literal line text.
   *This is why `Description="Approximate read depth; some reads may have been
   filtered"` precedes `Description="Read depth"` for the two `INFO/DP` lines, and
   why the `INFO/AD` line sits between `AC` and `AF`.*
3. `VCFContigHeaderLine` is the single exception: contig lines sort among themselves
   by **`contigIndex`**, the index assigned while the *input* header was parsed,
   so they come out in the input's own contig order; against a non-contig line the
   base text comparison applies, so the contig block lands where `##contig=...`
   sorts: **after every `##INFO=`, before `##source=`**.
4. The container is a `TreeSet`, so two lines with equal `toString()` collapse —
   harmless here because the header *set* already de-duplicates by
   `equals()` (`getClass() && mKey.equals && mValue.equals`, disassembled).
5. `#CHROM` is written last, outside the sort.
6. Where GATK's own lines land: `##GATKCommandLine` is an ordinary line and sorts
   between `##FORMAT=` (F < G) and `##INFO=` (G < I) — measured at index 10 of 31
   on every fixture; `##source=GenotypeGVCFs` sorts last, after the contigs (its
   key `source` > `contig`); `##FILTER=<ID=LowQual,...>` at index 2, between `ALT`
   and `FORMAT`; the block's absolute index is therefore a consequence of the sort,
   not of the group.

### 1.3 Measurement (`.diag/header_order_probe.py`)

4 fixtures, byte-identical declarations, only the input order differs; pinned GATK
4.6.2.0 GenotypeGVCFs vs native `--gatk-compatible-annotations`.

| fixture | input header order | prediction `fileformat + sorted(rest) + #CHROM` reproduces GATK | native mismatches BEFORE | AFTER |
| --- | --- | --- | --- | --- |
| A `contig-first` (the registered fixture shape) | contig, ALT, INFO DP/AD, FORMAT… | **True** | **14** | **0** |
| B `already-sorted` | ALT, FORMAT…, INFO AD/DP, contig | **True** | **5** | **0** |
| C `reversed` | FORMAT…, INFO, ALT, contig | **True** | **14** | **0** |
| D (gate) + `--include-non-variant-sites` | contig, ALT, … | **True** | **14** | **0** |

Log `.diag/header-order-probe.log` (BEFORE, 16:32) / `.diag/header-order-probe-after.log`
(AFTER, 16:38).  The probe writes to a fixed JSON path, so
`.diag/header-order-probe.json` now holds the AFTER state; the BEFORE numbers are
in that run's `.log`.

### 1.4 The contig control (`.diag/header_order_contig_probe.py`)

The `VCFContigHeaderLine` override was measured, not assumed, because "sort by
full text" alone would reorder contigs of a real genome (`chr1 < chr10 < chr2`):

* reference dictionary `chr2, chr10, chr1` + input header declaring the same order
  → GATK emits `chr2, chr10, chr1` (`gatk_sorts_contigs_by_text: false`), native
  the same;
* reference dictionary `chr2, chr10, chr1` + input header declaring
  `chr1, chr10, chr2` → GATK **accepts** and emits `chr1, chr10, chr2`, i.e. the
  **input's** order, not the reference dictionary's.

So `contigIndex` = input declaration order, and native's "keep the input's contig
order" is already right; only the block's *position* had to change.  Log
`.diag/header-order-contig-probe.log`.

## 2. STEP 2 — blast radius

### 2.1 Which outputs change

`gatk_compatible_header_text()` has exactly two call sites, `genotype_gvcf_tool.cpp`
`:5941` (streaming writer) and `:6783` (aggregate writer), **both guarded by
`options.gatk_annotation_compatibility`**.  Therefore:

| output | changes? |
| --- | --- |
| native GenotypeGVCFs, `--gatk-compatible-annotations`, VCF/VCF.gz, both writers | **line ORDER of the `##` block only** (content already matched) |
| native GenotypeGVCFs default (non-compat diagnostic) profile | **no** — the function is not reached |
| native GenotypeGVCFs **data rows** | **no** — records are written by `bcf_write(output, output_header, record)` from the untouched in-memory `bcf_hdr_t`; only the text header written by `write_vcf_text_line` is re-sorted |
| every other tool, the kernels, the Kokkos ABI | **no** — nothing else calls the helper |

The change is therefore confined to output serialization.  (One residual: with
`-O out.bcf --gatk-compatible-annotations` the compat path already writes raw text
into the stream, before and after this round; no registered test exercises it — §7.)

### 2.2 The named order-sensitive tests, and the count

Method: every script that feeds native GenotypeGVCFs output through
`--gatk-compatible-annotations`, plus a repo-wide grep for the signature of an
order-sensitive header comparison (a `##GATKCommandLine=`/`##source=` ignore list,
`read_text()` whole-file equality, or a position-wise list comparison).

* **18 scripts pass `--gatk-compatible-annotations`.**  Classified:
  * **ORDER-SENSITIVE — 1 registered test:**
    `fastgatk-hc-dense-gvcf-genotype-gatk-oracle`
    (`verify_hc_dense_gvcf_genotype_gatk_oracle.py:40-49`, `:127-128`):
    `normalized_joint_text()` drops only `##GATKCommandLine=`, `##source=`,
    `##fileDate=`, `##contig=` and `##fastgatk_genotype_gvcfs_status=`, then
    compares **the remaining lines in order** plus every data row.
    *Why it is expected to survive:* its input gVCF is a GATK-HaplotypeCaller
      gVCF, i.e. already htsjdk-sorted, so native's output already coincides with
      the sorted order for every line that test keeps — and the lines whose
      position this round moves (`##contig=`) are ignored there.
      **Run explicitly: exit 0** (`DENSE_EXIT=0`,
      `.diag/header-order-dense-hc-oracle.log`), plus green in step 5(d).
  * **ORDER-INSENSITIVE — 17 scripts**, all of which compare header *sets*,
    header *ID sets*, or data rows only:
    `verify_gatk_genotype_gvcf.py:26-33` (`sorted({...})`),
    `verify_gatk_genotype_gvcf_legacy_qual.py:33-37` (`sorted({...})`),
    `verify_gatk_genotype_gvcf_multisample.py:33-37` (`sorted({...})`),
    `verify_gatk_genotype_gvcf_inbreeding.py:30-40` (ID set),
    `verify_gatk_genotype_gvcf_multiallelic.py:31-38` (ID set),
    `verify_genotype_gvcf_assignment_gatk_oracle.py:27-40`, `:114-138`
    (ID subset), `verify_genotype_gvcf_exclude_intervals_gatk_oracle.py`,
    `verify_genotype_gvcf_gp_input_gatk_oracle.py`,
    `verify_genotype_gvcf_include_non_variant_gatk_oracle.py` (docstring:
    "compares the complete non-header …"),
    `verify_genotype_gvcf_malformed_gatk_oracle.py:302-306` (rows),
    `verify_genotype_gvcf_max_alternate_alleles_gatk_oracle.py:17` (rows),
    `verify_genotype_gvcf_multisample_reference_confidence_oracle.py`,
    `verify_genotype_gvcf.py`, `verify_genotype_gvcf_spandel_gatk_oracle.py`
    (FILTER group in order + content as multiset, order-only reported),
    `verify_genotype_gvcf_spanning_deletion_gatk_oracle.py:22-26` (rows),
    `verify_genotype_gvcf_starts_in_intervals_gatk_oracle.py:28` (rows),
    `verify_hc_chr20_max_mnp_polyploid_gvcf_gatk_oracle.py:23,29` (rows).
* **5 further scripts use the genotype binary without the compat flag**
  (`verify_genomicsdb_bridge.py`, `verify_genomicsdb_import.py`,
  `verify_genomicsdb_import_sample_map_gatk_oracle.py`,
  `verify_genomicsdb_import_update_workspace_gatk_oracle.py`,
  `verify_optional_boolean_contract.py`) — not reached by the change; the two that
  read a header read only the `#CHROM` line.
* The repo-wide `##GATKCommandLine=` grep returns 25 hits, all in Mutect2/HC/
  filter-mutect scripts that never consume native GenotypeGVCFs output.

**Count: 1 order-sensitive registered test, named above; it is explainable and
was run explicitly.  The change is confined to output serialization.  Therefore
the fix was implemented (step 4) rather than deferred.**

## 3. STEP 3 — the gate, written and run BEFORE the fix

`fastgatk-native/scripts/verify_genotype_gvcf_header_order_gatk_oracle.py` (new,
**not registered** — `fastgatk-native/CMakeLists.txt` must not be edited this
round; the orchestrator registers it).  It runs pinned GATK and native on four
fixtures with byte-identical declarations in four different input orders (dense
mode included), normalises away the producer-identity lines, and then:

1. re-derives htsjdk's comparator **from GATK's own bytes** and asserts GATK still
   satisfies it (so the gate cannot pass because the model drifted);
2. asserts the identity-stripped headers of GATK and native are byte-identical
   **position by position**, `##fileformat` first, `#CHROM` last;
3. asserts both sides carry the contig lines in the input's declaration order.

Identity normalisation (`is_identity()`): every line whose key is
`GATKCommandLine`, `source` or `fileDate`, plus any line carrying a
`CommandLine="` or `Date="` field.  Those record the program, its command line and
its wall-clock date; native is a different program, and comparing them would be a
test of the producer, not of the file.

### 3.1 Exit status and the literal diff, BEFORE the fix

**exit 1**, `"status": "divergence"`, 4 cases, **4 violations** — one per case, all
of them the position assertion (`.diag/header-order-gate-before.log`):

```
[header-order-star-contig-first]   mismatched positions: 14   comparator reproduces GATK: True
[header-order-star-input-already-sorted] mismatched positions: 5    comparator reproduces GATK: True
[header-order-star-input-reversed] mismatched positions: 14   comparator reproduces GATK: True
[header-order-star-dense]          mismatched positions: 14   comparator reproduces GATK: True
```

The literal position-by-position diff for `header-order-star-contig-first`
(identity lines already removed on both sides; 29 compared lines each):

```
  [ 1] GATK   ##ALT=<ID=NON_REF,Description="Represents any possible alternate allele">
       NATIVE ##contig=<ID=chr1,length=100>
  [ 2] GATK   ##FILTER=<ID=LowQual,Description="Low quality">
       NATIVE ##ALT=<ID=NON_REF,Description="Represents any possible alternate allele">
  [ 3] GATK   ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
       NATIVE ##FILTER=<ID=LowQual,Description="Low quality">
  [ 4] GATK   ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
       NATIVE ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
  [ 5] GATK   ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
       NATIVE ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
  [ 6] GATK   ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
       NATIVE ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
  [ 7] GATK   ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
       NATIVE ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
  [ 8] GATK   ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
       NATIVE ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
  [ 9] GATK   ##FORMAT=<ID=RGQ,...>
       NATIVE ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
  [10] GATK   ##INFO=<ID=AC,Number=A,...>
       NATIVE ##FORMAT=<ID=RGQ,...>
  [11] GATK   ##INFO=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
       NATIVE ##INFO=<ID=AC,Number=A,...>
  [15] GATK   ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
       NATIVE ##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
  [16] GATK   ##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
       NATIVE ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
  [27] GATK   ##contig=<ID=chr1,length=100>
       NATIVE ##INFO=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
```

`header-order-star-input-already-sorted` differs at only 5 positions (10, 11, 15,
16, 27) — the `INFO/AD` position, the `INFO/DP` pair and `##contig` — which is
exactly the "the input is already sorted, but native still moves those lines"
control.

### 3.2 AFTER the fix

**exit 0**, `"status": "pass"`, 0 violations, **0 mismatched positions in all four
cases**, comparator still reproducing GATK's own order in all four
(`.diag/header-order-gate-after.log`).

The independent probe of the previous round
(`.diag/header_lines_probe.py`, 9 fixture/mode combinations) now reports
`(a) missing in native 0`, `(b) extra in native 0`, `(c) same key different text 0`
and **`(d) order-only 0`** for all eight real fixtures (was 12-18) —
`.diag/header-order-lines-probe-after.log|json`.  The ninth case
(`reordered-attribute-keys`) is the pre-existing control whose input GATK rejects
(exit 3), so GATK's output file is empty and every native line counts as "extra";
unchanged by this round.

## 4. STEP 4 — the fix

`git diff --stat` = **1 file changed, 87 insertions(+)**
(`fastgatk-native/src/genotype_gvcf_tool.cpp`; the new gate script is a new,
untracked file).  No restructuring, no kernel ABI, no `CMakeLists.txt`, no
Mutect2, no other tool, no root `*.md`.

```diff
+std::vector<std::string> gatk_htsjdk_sorted_header(std::vector<std::string> lines) {
+    const std::size_t total = lines.size();
+    std::vector<std::string> fileformat, chrom, contigs, others;
+    for (auto& line : lines) {
+        if (line.rfind("##fileformat=", 0) == 0)      fileformat.push_back(std::move(line));
+        else if (line.rfind("#CHROM", 0) == 0)        chrom.push_back(std::move(line));
+        else if (line.rfind("##contig=", 0) == 0)     contigs.push_back(std::move(line));
+        else                                          others.push_back(std::move(line));
+    }
+    std::sort(others.begin(), others.end());
+    const std::string contig_prefix = "##contig=";
+    const auto at = std::lower_bound(others.begin(), others.end(), contig_prefix);
+    std::vector<std::string> sorted;
+    sorted.reserve(total);
+    sorted.insert(sorted.end(), fileformat.begin(), fileformat.end());
+    sorted.insert(sorted.end(), others.begin(), at);
+    sorted.insert(sorted.end(), contigs.begin(), contigs.end());
+    sorted.insert(sorted.end(), at, others.end());
+    sorted.insert(sorted.end(), chrom.begin(), chrom.end());
+    return sorted;
+}
 ...
 std::string gatk_compatible_header_text(...) {
 ...
     if (!inserted_source) { ... output.insert(chrom, "##source=GenotypeGVCFs"); }
+    // GATK's output header is htsjdk's SORTED header, not the input's order with
+    // the tool's own lines appended ...
+    output = gatk_htsjdk_sorted_header(std::move(output));
     std::ostringstream result;
     for (const auto& line : output) result << line << '\n';
```

Why this is the correct and minimal shape:

* it is a **literal transcription** of the disassembled comparator: `##fileformat`
  pinned first (item 1); the rest `std::sort`ed by full line text (item 2, and
  `std::string::operator<` is byte order = Java `String.compareTo` for the ASCII
  lines measured); contigs taken out and re-inserted as a block, in their existing
  (input) order, at the `##contig=` text position (item 3 — `std::lower_bound` on
  the prefix is exactly "the first non-contig line that sorts after any contig
  line", because every contig line starts with that prefix and no other line does);
  `#CHROM` last (item 4/5);
* it is applied **after** the existing group insertions, which therefore still
  decide *which* lines are present and are now positionally inert — the smallest
  possible edit (one new function, one new call, no removal);
* it lives wholly inside `gatk_compatible_header_text()`, the
  `--gatk-compatible-annotations` boundary used by both writers, so the diagnostic
  profile, the genotyping engine, the FILTER-column handling from the previous
  round, the AF/PL kernels and their ABI, and every other tool are untouched;
* the duplicate `##INFO=<ID=DP,...>` and `##FORMAT=<ID=AD,...>` pairs that the
  previous round reproduced come out in GATK's order automatically, because the
  sort is by full text (`"Approximate read depth; …"` < `"Read depth"`), which was
  one of the two remaining visible symptoms.

## 5. STEP 5 — gate results (mandatory order)

| Step | Command | Result |
| --- | --- | --- |
| a | `python3 fastgatk-native/scripts/verify_genotype_gvcf_header_order_gatk_oracle.py` | **exit 0** (4 cases, 0 violations, 0 mismatched positions); **exit 1 before** the fix (4 violations, §3.1) |
| b1 | `python3 fastgatk-native/scripts/verify_genotype_gvcf.py` | `B1_EXIT=0` — `{"status": "pass", "output_records": 1}` (`.diag/header-order-verify-gvcf.log`) |
| b2 | `ctest --test-dir fastgatk-native/build -R 'genotype-gvcf' -V` | `B2_EXIT=0` — **16/16 passed** (`.diag/header-order-ctest-genotype-gvcf.log`) |
| c | the 19-name strict-gate filter from the task | `C_EXIT=0` — **21/21 passed** (`.diag/header-order-strict-gates.log`) |
| d | `fastgatk-native/scripts/run_regression.sh --label header-order` | `REG_EXIT=0` — **302/302 on BOTH backends** (§5.1) |
| extra | the one order-sensitive registered test, run directly | `DENSE_EXIT=0` — `{"status": "pass", "joint_normalized_full_text_exact": true}` (`.diag/header-order-dense-hc-oracle.log`) |

### 5.1 Dual-backend regression (d)

```
REG_EXIT=0
| 后端   | 构建目录                              | 结果 | 通过/总数 | 耗时        |
| omp    | OpenMP (fastgatk-native/build)        | 通过 | 302/302   | 1494.89 sec |
| serial | Serial (fastgatk-native/build-serial) | 通过 | 302/302   | 1381.02 sec |
```

Evidence block `.diag/regression/20260912-171711/{omp.log,serial.log,evidence.md}`
(+ `.diag/header-order-regression.out` for `REG_EXIT`), run window
17:17:11 → 17:42:06, label `header-order`, git `741a110` + the two uncommitted
files of §8.  Both backends ran the full 302-test suite with
`FASTGATK_REQUIRE_GATK_ORACLE=1` (the runner's default) and ctest parallelism 8;
the runner's staleness check reported no warning (`grep -c 警告` = 0 in both logs),
i.e. no source file is newer than the binaries it tested.

**Nothing was edited after the builds.**  Provenance (timestamps):

```
16:35:43  fastgatk-native/scripts/verify_genotype_gvcf_header_order_gatk_oracle.py (new)
16:36:39  fastgatk-native/src/genotype_gvcf_tool.cpp
16:37:11  fastgatk-native/build/fastgatk-genotype-gvcf          (omp, rebuilt)
16:37:24  fastgatk-native/build-serial/fastgatk-genotype-gvcf   (serial, rebuilt)
17:17:11 → 17:42:06   regression window (both backends)
```

`md5sum` (also `.diag/header-order-binaries.md5`):

```
8faefce59dd1e9e6bac7b7e5ba9f077a  fastgatk-native/src/genotype_gvcf_tool.cpp
e7fb159c8ce156c781bd958c84c0e8c8  fastgatk-native/scripts/verify_genotype_gvcf_header_order_gatk_oracle.py
a612bc85d2e99afd62d540345630220d  fastgatk-native/build/fastgatk-genotype-gvcf        (omp)
bbf534f59bc46227129e431b6ea4121c  fastgatk-native/build-serial/fastgatk-genotype-gvcf (serial)
```

The binaries were built (16:37) from the frozen source (16:36) and before every
step-5 run; the only files written afterwards are this report and the read-only
gate logs/JSON of step 5, none of which is a build or test input.  So no rebuild
and no re-run of the suite were required.

## 6. No stale assertions

No registered test pinned the old (input) order, so no test line needed
correcting; no existing test was weakened or edited.  The only new test file is
the gate itself, and the only source change is §4.  §2.2 lists every script
checked and why each is order-insensitive or unaffected; §5 is the empirical
check (302/302 both backends).

## 7. What remains unproven

* The gate pins the order on 4 fixtures/3 input orders/single-contig references.
  Multi-contig ordering is measured only by `.diag/header_order_contig_probe.py`
  (2 orders, 3 contigs) and is **not** part of the gate: the gated fixtures have
  one contig.  Untested with the comparator: an input that declares **no**
  `##contig` line at all, an input with non-ASCII bytes in a `Description` (Java
  compares UTF-16 code units, C++ compares UTF-8 bytes — they agree on ASCII and
  can disagree outside it), and a header carrying two lines with identical text.
* The `contigIndex` source is established by two measurements (input order, not
  reference-dictionary order).  I did not find the htsjdk call site that assigns
  it — `VCFHeaderLineTranslator.parseLine` in the shaded jar takes a
  `List<String>`, not an index, so the assignment is inside a class not present in
  the local jar.  If htsjdk ever used the dictionary index instead, native's
  "keep the input's order" would be wrong for an input whose contig order differs
  from the reference dictionary — measured once, in native's favour, but not
  proven from source.
* The non-compat diagnostic profile is claimed unchanged **by construction** (the
  function is unreachable without `--gatk-compatible-annotations`); it was not
  byte-compared before/after (the pre-fix binary was overwritten by the rebuild).
* `-O out.bcf --gatk-compatible-annotations` is not exercised by any registered
  test; that path already wrote raw text into the stream before this round, and
  this round neither fixed nor worsened it.
* §2.2's "order-insensitive" classification is a source reading plus the step-5
  results, not a proof: a test could build such an expectation dynamically.  Step
  5(d) is the empirical check, which is evidence, not proof.
* The `##FILTER` group's absolute index (GATK 2, native 3) was **not** separately
  measured this round; the sort moves the group to the `FILTER` key position
  (after `##ALT`, before `##FORMAT`), which is index 2 in GATK's output and is
  asserted positionally by the new gate on all four fixtures.

## 8. Tree state

* Change **kept in the tree** (no commit, no branch, no branch switch).
* Modified: `fastgatk-native/src/genotype_gvcf_tool.cpp` (`+87/-0`).
  New (untracked): `fastgatk-native/scripts/verify_genotype_gvcf_header_order_gatk_oracle.py`.
* NOT touched: `fastgatk-native/CMakeLists.txt`, every root `*.md`, Mutect2,
  every other tool, `third_party/`, and every registered test script.
* Artifacts: probes `.diag/header_order_probe.py`, `.diag/header_order_contig_probe.py`
  (+ logs/JSON `header-order-probe*.log|json`, `header-order-contig-probe.log`),
  gate logs `.diag/header-order-gate-{before,after}.log`, build logs
  `.diag/header-order-build-{omp,serial}.log`, step-5 logs
  `.diag/header-order-{verify-gvcf,ctest-genotype-gvcf,strict-gates,dense-hc-oracle}.log`,
  the previous round's probe re-run `.diag/header-order-lines-probe-after.{log,json}`,
  regression evidence `.diag/regression/20260912-171711/`, hashes
  `.diag/header-order-binaries.md5`.
