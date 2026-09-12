# Round 49 — Multi-input header divergence of GenotypeGVCFs (measurement only)

Date: 2026-09-13
Scope: **measurement and root-cause localisation only. No production file was modified.**
Artifacts fixed by this round:

| artifact | md5 |
|---|---|
| `fastgatk-native/build/fastgatk-genotype-gvcf` (OpenMP) | `4a65c0abb7db3424ecce8a7d0c2afdc2` |
| `fastgatk-native/build-serial/fastgatk-genotype-gvcf` (serial) | `a79162808e64f6d39ceca93eed81ad54` |
| `third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar` | `aef501a2e7c1f70622061799025dbd69` |
| `.diag/round49/ref/tiny.fa` (synthetic 2-contig reference) | `98d050b625f9d24bdd6203483f11560a` |

Scratch: `.diag/round49/` (inputs, outputs, scripts, logs).
Raw comparator log: `fastgatk-native/evidence/2026-09-13-round49/compare_all.txt`.

---

## 0. Headline: the suspicion as stated cannot be tested against the pinned oracle, because the pinned oracle refuses more than one `-V`

**Measured, not inferred.** GATK 4.6.2.0 `GenotypeGVCFs` is a `VariantLocusWalker`, and that
walker declares the driving variant input as a **single, non-repeatable** argument:

`gatk-source/src/main/java/org/broadinstitute/hellbender/engine/VariantLocusWalker.java:33-35`

```java
// NOTE: using String rather than FeatureInput<VariantContext> here so that we can keep this driving source
//       of variants separate from any other potential sources of Features
@Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME, shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME, doc = "A VCF file containing variants", common = false, optional = false)
public String drivingVariantFile;
```

and `gatk-source/src/main/java/org/broadinstitute/hellbender/tools/walkers/GenotypeGVCFs.java:112`
declares `public final class GenotypeGVCFs extends VariantLocusWalker`, with the class javadoc
(lines 65 and 99) stating verbatim:

> The GATK4 GenotypeGVCFs tool can take only one input track.
> ...
> * Cannot take multiple GVCF files in one command.

Runtime confirmation (exact command, run for every one of the 18 case pairs/triples that went through the comparison runner):

```
$ third_party/jdk17/bin/java -Xmx1g -jar third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
    GenotypeGVCFs -R .diag/round49/ref/tiny.fa \
    -V .diag/round49/in/t1_baseline/in1.g.vcf -V .diag/round49/in/t1_baseline/in2.g.vcf \
    -O /dev/null --create-output-variant-index false
A USER ERROR has occurred: Illegal argument value: Argument 'V/variant' cannot be specified more than once.
exit status 1
```

This was reproduced for all 18 of them (t1, t1b, t2..t17): `gatk-multiV <case> exit=1` in every case,
including pairs whose headers are byte-identical apart from the sample name. **So the premise
"GATK merges the headers of all inputs" is only reachable through a different entry point.**

The reachable GATK merge oracle is the multi-input data source that `CombineGVCFs`,
`GenomicsDBImport` and every other `MultiVariantWalker` use —
`gatk-source/src/main/java/org/broadinstitute/hellbender/engine/MultiVariantDataSource.java`
(lines 237-242, 249-260):

```java
/** Merge and sort the samples from each header requiring unique samples */
private SortedSet<String> getSortedSamples() {
    final Map<String, VCFHeader> headers = featureDataSources.stream()
            .collect(Collectors.toMap(ds -> ds.getName(), ds -> (VCFHeader) ds.getHeader()));
    return VcfUtils.getSortedSampleSet(headers, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);
}

private VCFHeader getMergedHeader() {
    ... 
    // Now merge the headers using htsjdk, which is pretty permissive, and which only works properly
    // because of the cross-dictionary validation done in validateAllSequenceDictionaries.
    return headers.size() > 1 ?
            new VCFHeader(VCFUtils.smartMergeHeaders(headers, true)) : headers.get(0);
}
```

I therefore used, as the reference oracle for "what a merge does":

1. `CombineGVCFs -V in1 -V in2 [-V in3]` — the end-to-end GATK multi-input path whose output
   header **is** `MultiVariantDataSource.getMergedHeader()`; and
2. `GenotypeGVCFs -V <combined.g.vcf>` — the canonical GATK cohort workflow
   (`CombineGVCFs` → `GenotypeGVCFs`), which yields the final user-visible VCF; and
3. a direct in-process probe of `VCFUtils.smartMergeHeaders(headers, true)` and
   `VcfUtils.getSortedSampleSet(map, REQUIRE_UNIQUE)` out of the **pinned jar itself**
   (`.diag/round49/Probe.java`), so the merge contract is measured from GATK's own bytecode,
   not quoted from memory.

The GATK startup banner in every run confirms the bundled library version:
`INFO CombineGVCFs - HTSJDK Version: 4.2.0`. I disassembled the bundled
`htsjdk.variant.vcf.VCFUtils.smartMergeHeaders` and diffed it against the htsjdk 4.2.0 jar found on
this machine; the only differences are constant-pool index widths / `ldc` vs `ldc_w` and branch
offsets — the control flow is identical. All htsjdk statements below are read off the **bundled**
class (`javap -p -c`), with 4.2.0 as the corroborating upstream version.

---

## 1. Method

### 1.1 Reference and inputs

Small synthetic reference (not the 63 Mbp chr20), 2 contigs, built so that header edits are
visible and runs are seconds long:

```
>chr1   length 1000
>chr2   length 800
```
`ref/tiny.fa` + `.fai` + `tiny.dict`. Inputs are 1-sample gVCFs, 3 records each on `chr1`
(ref block 1-80, het SNP at 81 with `ALT=<ALT>,<NON_REF>`, ref block 82-200), or on `chr2` for the
disjoint cases. Generator: `.diag/round49/make_inputs.py`, `.diag/round49/make_extra_inputs.py`.

### 1.2 Case matrix (19 fixture directories under `.diag/round49/in/<case>/`; 18 were run through the
full native + GATK comparison runner, `t18` was probe-only)

| case | axis deliberately injected |
|---|---|
| `t1_baseline` | identical headers except sample name (`SAMPLE_A` vs `SAMPLE_B`) |
| `t1b_three` | three inputs, `SAMPLE_A/B/C` |
| `t2_extra_info` | input 2 declares + uses `##INFO=<ID=EXINFO,...>`; input 1 declares a different extra `ONLY1INFO` |
| `t3_extra_format` | input 2 declares + uses `##FORMAT=<ID=EXFMT,...>`; input 1 declares `ONLY1FMT` |
| `t4_source` | different `##source=` lines (`fastgatk-round49-inputA` vs `SomeOtherCaller_v9.9`) |
| `t5_collision` | same sample name `SAME_NAME` in both inputs |
| `t6_no_nonref_decl` | input 2 missing `##ALT=<ID=NON_REF,...>` |
| `t15_nonref_in_second` | `##ALT=<ID=NON_REF,...>` **only** in input 2 |
| `t7_contig_length` | input 2 declares `chr1` length 900 vs input 1's 1000 |
| `t8_contig_order` | input 2 declares `(chr2, chr1)`, input 1 `(chr1, chr2)` |
| `t16_contig_order_first_reversed` | input 1 declares `(chr2, chr1)` — non-sorted first input |
| `t9_extra_contig` | `chr2` declared only in input 2 |
| `t10_three_mixed` | 3 inputs: extra FILTER in in1, extra INFO in in2, extra FORMAT + duplicate `SAMPLE_A` in in3 |
| `t11_extra_filter` | input 2 declares extra `##FILTER=<ID=Only2Filter,...>` |
| `t14_sample_order` | samples in non-lexicographic input order: `ZEBRA_SAMPLE` then `ALPHA_SAMPLE` |
| `t17_collision_diffdata` | same sample name, deliberately different PL payloads |
| `t12_disjoint_contigs` | `SAMPLE_A` only on `chr1`, `SAMPLE_B` only on `chr2` |
| `t13_disjoint_intervals` | `SAMPLE_A` on `chr1:1-200`, `SAMPLE_B` on `chr1:201-400` |
| `t18_version_mix` | `##fileformat=VCFv4.2` vs `VCFv4.3` (probe only) |

### 1.3 Commands (per case, `-V` repeated in the order the shell glob produced, i.e. `in1, in2[, in3]`)

```
NATIVE  fastgatk-native/build/fastgatk-genotype-gvcf -R <ref> -V in1.g.vcf -V in2.g.vcf \
            --gatk-compatible-annotations -O native.vcf
SERIAL  fastgatk-native/build-serial/fastgatk-genotype-gvcf  (same arguments)
GATK-M  java -Xmx1g -jar <jar> GenotypeGVCFs -R <ref> -V in1.g.vcf -V in2.g.vcf \
            -O gatk_multi_v.vcf --create-output-variant-index false
GATK-C  java -Xmx1g -jar <jar> CombineGVCFs  -R <ref> -V in1.g.vcf -V in2.g.vcf \
            -O gatk_combine.g.vcf
GATK-F  java -Xmx1g -jar <jar> GenotypeGVCFs -R <ref> -V gatk_combine.g.vcf \
            -O gatk_combined_then_gg.vcf --create-output-variant-index false
```

Runners: `.diag/round49/run_all.sh`, `.diag/round49/run_extra.sh`.
Comparator: `.diag/round49/compare.py`. The OpenMP and serial native binaries produced
**byte-identical** VCFs for all six extra cases (`cmp` reported IDENTICAL), so everything below is
not a scheduling artifact.

Header comparison is done on the htsjdk key rule read off `smartMergeHeaders`: the merge key of a
line is its key alone, or `key + "=" + ID` when the line implements `VCFHeaderIDLine`
(`contig=chr1`, `INFO=EXINFO`, `FORMAT=EXFMT`, `FILTER=Only2Filter`). Order is compared as the
relative order of the *shared* keys, and separately by quoting the raw files.

---

## 2. Ranked divergences (confirmed)

### D1 — Sample columns are emitted in INPUT order instead of GATK's sorted order  🔴 highest impact

This one silently permutes every data column of the VCF and is invisible in a single-sample run.

*Minimal reproducing pair*: `in1` header sample `ZEBRA_SAMPLE`, `in2` header sample `ALPHA_SAMPLE`
(case `t14_sample_order`), files otherwise byte-identical.

```
$ native -R ref -V in1.g.vcf -V in2.g.vcf --gatk-compatible-annotations -O native.vcf
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ZEBRA_SAMPLE	ALPHA_SAMPLE
chr1	81	.	C	T	500.50	.	AC=2;AF=0.500;AN=4;DP=30;ExcessHet=1.7609;MLEAC=2;MLEAF=0.500;QD=8.07	GT:DP:GQ:PL:MIN_DP	0/1:31:99:255,0,255:30	0/1:31:99:255,0,255:30

$ java ... CombineGVCFs -V in1.g.vcf -V in2.g.vcf -O gatk_combine.g.vcf     (the merge path)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ALPHA_SAMPLE	ZEBRA_SAMPLE

$ java ... GenotypeGVCFs -V gatk_combine.g.vcf -O gatk_combined_then_gg.vcf  (canonical GATK cohort path)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ALPHA_SAMPLE	ZEBRA_SAMPLE
chr1	81	.	C	T	500.50	.	AC=2;AF=0.500;AN=4;DP=60;ExcessHet=1.7609;MLEAC=2;MLEAF=0.500;QD=8.07	GT:AD:DP:GQ:PL	0/1:30,0:30:99:255,0,255	0/1:30,0:30:99:255,0,255
```

Required behaviour is explicit in the GATK source: `MultiVariantDataSource.getSortedSamples()`
returns a **`SortedSet`** (line 237-243, javadoc line 234-236 "Merge **and sort** the samples"), and
`VcfUtils.getSortedSampleSet` builds it with `new TreeSet<>()` and `addAll`. The probe confirms the
sorted result directly:

```
--- VcfUtils.getSortedSampleSet(map, REQUIRE_UNIQUE)  [MultiVariantDataSource.getSortedSamples]
  RESULT [ALPHA_SAMPLE, ZEBRA_SAMPLE]
```

Native instead assigns column indices by first-seen order
(`genotype_gvcf_tool.cpp:6533-6541`: `output_sample_indices.emplace(name, next)` with
`next = output_sample_indices.size()`).

All other 18 cases agree on the sample column set only because their names happen to be
lexicographically ordered by input order (`SAMPLE_A < SAMPLE_B < SAMPLE_C`). This is a near-miss,
not a pass.

### D2 — `##INFO` declarations that appear only in a later input are dropped  🔴 highest impact

*Minimal reproducing pair*: `t2_extra_info` — `in1` declares `##INFO=<ID=ONLY1INFO,...>`, `in2`
declares `##INFO=<ID=EXINFO,Number=1,Type=Integer,Description="Extra INFO only declared in input 2">`
and carries `EXINFO=5` on the chr1:81 record. Both are otherwise identical apart from the sample name.

```
native header: (no EXINFO anywhere)
gatk final header: ##INFO=<ID=EXINFO,Number=1,Type=Integer,Description="Extra INFO only declared in input 2">

native row:  chr1	81	.	C	T	500.50	.	AC=2;AF=0.500;AN=4;DP=30;ExcessHet=1.7609;MLEAC=2;MLEAF=0.500;QD=8.07	...
gatk   row:  chr1	81	.	C	T	500.50	.	AC=2;AF=0.500;AN=4;DP=60;EXINFO=5;ExcessHet=1.7609;MLEAC=2;MLEAF=0.500;QD=8.07	...
```

Comparator output (`compare_all.txt`, case `t2_extra_info`, "[native vs FINAL gatk VCF header]"):

```
      native-only keys (0):
      gatk-final-only keys (1):
         ##INFO=<ID=EXINFO,Number=1,Type=Integer,Description="Extra INFO only declared in input 2">
      same key, different text (0):
      shared-key relative ORDER identical (30 shared keys)
```

The declaration loss propagates into the data: because htslib will not emit a tag the output header
does not declare, the `EXINFO=5` annotation is silently lost. Proof that this is caused purely by
which input came first, not by htslib refusing the tag — the same two files with the `-V` order
swapped:

```
$ native -R ref -V in2.g.vcf -V in1.g.vcf --gatk-compatible-annotations -O native_swapped.vcf
$ grep EXINFO native_swapped.vcf
##INFO=<ID=EXINFO,Number=1,Type=Integer,Description="Extra INFO only declared in input 2">
chr1	81	.	C	T	500.50	.	AC=2;AF=0.500;AN=4;DP=30;ExcessHet=1.7609;MLEAC=2;MLEAF=0.500;QD=8.07;EXINFO=5	...
$ grep ONLY1INFO native_swapped.vcf
(absent)
$ grep '^#CHROM' native_swapped.vcf
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE_B	SAMPLE_A
```

Symmetric and decisive: the output header is the **first** input's header, and later inputs
contribute only sample names and contigs.

### D3 — `##FORMAT` declarations that appear only in a later input are dropped

*Minimal reproducing pair*: `t3_extra_format`. `in2` declares
`##FORMAT=<ID=EXFMT,Number=1,Type=Integer,...>` and carries `...:7` in the FORMAT column; `in1`
declares `ONLY1FMT` and carries nothing.

```
      native-only keys (0):
      gatk-final-only keys (1):
         ##FORMAT=<ID=EXFMT,Number=1,Type=Integer,Description="Extra FORMAT only declared in input 2">

native row FORMAT/fields: GT:DP:GQ:PL:MIN_DP   0/1:31:99:255,0,255:30   0/1:31:99:255,0,255:30
gatk   row FORMAT/fields: GT:AD:DP:EXFMT:GQ:PL 0/1:30,0:30:.:99:255,0,255   0/1:30,0:30:7:99:255,0,255
```

Note the GATK data consequence: the sample that does not carry `EXFMT` gets a literal `.`
(`0/1:30,0:30:.:99:...`), the one that does gets `7` — the union declaration makes mixed-coverage
FORMAT columns expressible at all. Native cannot express either (it drops the column entirely, for
the same reason as D2, and additionally drops the `7` value). Reproduced at 3 inputs in
`t10_three_mixed` (gatk-final-only: `EXFMT` **and** `EXINFO`).

### D4 — `##FILTER` declarations that appear only in a later input are dropped

*Minimal reproducing pair*: `t11_extra_filter`, `in2` declares
`##FILTER=<ID=Only2Filter,Description="Filter only declared in input 2">`.

```
      gatk-final-only keys (1):
         ##FILTER=<ID=Only2Filter,Description="Filter only declared in input 2">
```

Native emits only `##FILTER=<ID=LowQual,...>`, i.e. `in1`'s filter group. A record in a later input
that is filtered with `Only2Filter` would therefore be unrepresentable / would lose its FILTER value.

### D5 — The `##ALT=<ID=NON_REF,...>` declaration is not unioned from later inputs

The suspicion's own listed axis, tested in both directions:

* `t6_no_nonref_decl` (`in1` has `##ALT`, `in2` lacks it): **both agree** — the merged header keeps
  `in1`'s line, which is what GATK's first-wins rule does too (`gatk-final-only keys (0)`).
* `t15_nonref_in_second` (`##ALT` **only** in `in2`): **diverge**.

```
$ grep '^##ALT' native.vcf                                  (t15, -V in1 -V in2)
(none — native emits no ##ALT line at all)
$ grep '^##ALT' gatk_combined_then_gg.vcf
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
```

Comparator: `gatk-final-only keys (1): ##ALT=<ID=NON_REF,...>`. Native still writes `<NON_REF>` into
the ALT column of its records while declaring the symbolic allele nowhere — the output VCF is
internally inconsistent, not merely incomplete.

### D6 — Root cause confirmed: the output header is literally the first input's header

`--gatk-compatible-annotations` takes the non-streaming path (`Options::stream_by_locus` defaults to
`false`, `genotype_gvcf_tool.cpp:115`; `--stream-by-locus` is the only thing that sets it, line 662).
In that path the header is seeded once and never extended with declarations:

`fastgatk-native/src/genotype_gvcf_tool.cpp` — the per-input loop that begins at line 6512
(`for (const auto& input_path : input_paths) {`), with the header assembly at 6544-6614 (the loop
itself continues past 6619 into record processing):

```cpp
6544:            if (!output_header) {
6545:                output_header = bcf_hdr_dup(header);      // <-- FIRST input only, verbatim
...
6599:            }
6600:            // Merge contig dictionaries before materializing records. ...
6604:            merge_contig_dictionary(output_header, header);
6605:            for (int sample = 0; sample < header->n[BCF_DT_SAMPLE]; ++sample) {
...
6607:                if (bcf_hdr_id2int(output_header, BCF_DT_SAMPLE, name) < 0 &&
6608:                    bcf_hdr_add_sample(output_header, name) != 0) {
...
6613:            }
6614:            if (bcf_hdr_sync(output_header) != 0) {
```

Only two things are carried from inputs 2..N: contig lines (`merge_contig_dictionary`, defined at
line 995-1032) and sample names (`bcf_hdr_add_sample`, line 6605-6613). Nothing merges
`INFO` / `FORMAT` / `FILTER` / `ALT` / non-ID lines. The streaming path has the same shape at
lines 5871-5879 (`output_header = bcf_hdr_dup(header)`, then
`add_genotype_output_header_fields`), 5880, and 5901-5909.

The `-V`-swap experiment in D2 is the direct empirical confirmation of "first input only".

### D7 — Per-sample value for a sample with no input data at a locus: GATK `./.` vs native's full row

Explicitly requested as point (d), and it does interact with the earlier genotype-assignment work.

*Minimal reproducing pair*: `t12_disjoint_contigs` (`SAMPLE_A` has records only on `chr1`,
`SAMPLE_B` only on `chr2`); `t13_disjoint_intervals` is the same-contig variant.

```
        native rows:
chr1	81	.	C	T	247.64	.	AC=1;AF=0.500;AN=2;DP=30;...	GT:DP:GQ:PL:MIN_DP	0/1:31:99:255,0,255:30	./.:.:.:.,.,.:.
chr2	50	.	T	A	247.64	.	AC=1;AF=0.500;AN=2;DP=30;...	GT:DP:GQ:PL:MIN_DP	./.:.:.:.,.,.:.	0/1:31:99:255,0,255:30

        gatk GenotypeGVCFs(Combine(...)) rows:
chr1	81	.	C	T	247.64	.	AC=1;AF=0.500;AN=2;DP=30;...	GT:AD:DP:GQ:PL	0/1:30,0:30:99:255,0,255	./.
chr2	50	.	T	A	247.64	.	AC=1;AF=0.500;AN=2;DP=30;...	GT:AD:DP:GQ:PL	./.	0/1:30,0:30:99:255,0,255
```

**Answer: GATK writes a bare `./.`** — one token, no colons, no missing FORMAT placeholders. Native
writes the fully expanded `./.:.:.:.,.,.:.`, i.e. GT=`./.` plus five empty slots, and the PL slot is
rendered as `.,.,.` (three empty sub-values). Same result in `t13`:
`./.:.:.:.,.,.:.` (native) vs `./.` (GATK) at `chr1:81` and `chr1:300`.

A bare `./.` means "this sample has no data at this locus" and is *not* the same concept as a
no-call genotype with missing annotations, which is what native emits. Downstream parsers that
count FORMAT sub-fields, or that distinguish "no data" from "no-call with missing values", will
read these two files differently.

### D8 — Interface divergence: native accepts 2+ `-V` at all

Independent of any header content. Native accepts `-V` repeated
(`genotype_gvcf_tool.cpp:5832` / `6495`, `expand_genomicsdb_inputs(options.inputs, options.regions)`),
and its `--help` advertises `-V, --variant FILE  input VCF/GVCF (repeatable)`. The pinned GATK errors
out with `Argument 'V/variant' cannot be specified more than once.` So on this axis native is a
superset of the oracle, and any multi-`-V` invocation has **no GATK counterpart to be compatible
with** — the only defensible compatibility target is the merged header semantics reachable through
`CombineGVCFs` / `MultiVariantDataSource`.

(What native's multi-input support actually mirrors is the GenomicsDB shard case: for
`gendb://` inputs `expand_genomicsdb_inputs` fans a workspace out into per-shard VCFs whose headers
are identical by construction, because `GenomicsDBImport` writes one merged header. The divergences
above only become observable when a caller passes genuinely different VCFs — which is exactly what
the native `-V` repeatable interface invites.)

---

## 3. Axes where the suspicion was REFUTED (native already agrees with GATK)

Negative results, stated explicitly. All of these were measured, not reasoned.

### N1 — `##source`: GATK does **not** union differing `source` lines either  ✅ agree

*Pair*: `t4_source` (`##source=fastgatk-round49-inputA` vs `##source=SomeOtherCaller_v9.9`).

```
native:                ##source=GenotypeGVCFs
                       ##source=fastgatk-round49-inputA
gatk CombineGVCFs:     ##source=CombineGVCFs
                       ##source=fastgatk-round49-inputA
gatk final:            ##source=CombineGVCFs
                       ##source=GenotypeGVCFs
                       ##source=fastgatk-round49-inputA
```

`SomeOtherCaller_v9.9` is absent from GATK's output too. Comparator: `gatk-combine-only keys (0)`,
`same key, different text (0)`. The direct probe of `smartMergeHeaders` returns exactly 15 lines for
these two headers, and 15 is precisely `in1`'s line count — the second header contributes nothing on
this axis.

Mechanism, read off the bundled bytecode (`htsjdk.variant.vcf.VCFUtils.smartMergeHeaders`): a plain
(non-ID) `VCFHeaderLine` whose key is already present but whose value differs falls into the final
`else` branch, which only calls `HeaderConflictWarner.warn(...)` and then **does not** `map.put` —
first value wins, later values are dropped with a warning. (`CombineGVCFs`'s stderr in `t4` shows no
warning line, so the warning itself is a `java.util.logging` record I did not capture; the line set
is the measurement.)

### N2 — Contig lines from later inputs **are** merged  ✅ agree

*Pair*: `t9_extra_contig` (`chr2` declared only in `in2`).

```
native:  ##contig=<ID=chr1,length=1000>  ##contig=<ID=chr2,length=800>
gatk:    ##contig=<ID=chr1,length=1000>  ##contig=<ID=chr2,length=800>
```

`merge_contig_dictionary` (line 995-1032) does union later inputs' contigs, keeping first-seen
lengths. This is the one declaration group native already merges — so the fix for D2-D5 is a
generalisation of an existing, working pattern, not new machinery.

### N3 — Contig ORDER: both preserve the FIRST input's declaration order  ✅ agree

Tested in both directions, and this is where I expected a divergence and found none.

* `t8_contig_order`: `in1=(chr1,chr2)`, `in2=(chr2,chr1)` → native `chr1, chr2`; gatk final `chr1, chr2`.
* `t16_contig_order_first_reversed`: `in1=(chr2,chr1)`, `in2=(chr1,chr2)` → native `chr2, chr1`;
  gatk final `chr2, chr1`.

Neither side sorts contigs; both inherit the first input's order. `t16` exists specifically to
falsify "native keeps input order, GATK sorts" — it did not falsify into a divergence.

### N4 — Contig-length conflict: both **reject** the run  ✅ agree (different error text)

*Pair*: `t7_contig_length` (`chr1` = 1000 in `in1`, 900 in `in2`).

```
$ native ... -V in1.g.vcf -V in2.g.vcf ...
fastgatk-genotype-gvcf: BAD_INPUT: contig length mismatch across GenotypeGVCFs shards: chr1
exit status 2

$ java ... CombineGVCFs -V in1.g.vcf -V in2.g.vcf ...
A USER ERROR has occurred: Input files .../in2.g.vcf and .../in1.g.vcf have incompatible contigs:
Incompatible sequences found (chr1: 900) and (chr1: 1000).
exit status 2
```

GATK enforces this before merging
(`MultiVariantDataSource.validateAllSequenceDictionaries()` / `validateSequenceDictionaryRecords`,
lines 283-339, throwing `UserException.IncompatibleSequenceDictionaries`); native throws its own
`BAD_INPUT` from `merge_contig_dictionary` (lines 1016-1021). Both refuse to produce output, so on
the "does a length conflict get reconciled or refused?" axis the answer is the same: refused.
Only the message/exception taxonomy differs.

### N5 — Sample-name collision: both collapse to ONE column  ✅ agree

*Pairs*: `t5_collision` and `t17_collision_diffdata` (both inputs name their sample `SAME_NAME`).

```
native  #CHROM ... SAME_NAME
gatk    #CHROM ... SAME_NAME        (both CombineGVCFs and the final VCF)
```

I expected GATK to error here, because `MultiVariantDataSource` asks for
`GenotypeMergeType.REQUIRE_UNIQUE` and its javadoc says
"Merge and sort the samples from each header **requiring unique** samples". It does not error.
The probe measures why:

```
--- VcfUtils.getSortedSampleSet(map, REQUIRE_UNIQUE)
  RESULT [SAME_NAME]
```

`VcfUtils.getSortedSampleSet` builds a `TreeSet` and adds
`GATKVariantContextUtils.mergedSampleName(fileName, sampleName, uniquify)` per sample; the
disassembly shows `mergedSampleName(a, b, uniquify)` is `uniquify ? b + "." + a : b`, and with
`REQUIRE_UNIQUE` (`uniquify == false`) it returns the sample name unchanged, so the `TreeSet`
silently de-duplicates. No exception, no suffix, no warning. Native's own de-duplication
(`output_sample_indices`, lines 6533-6541) produces the same single column, so (a) agrees on the
collision axis.

Secondary, **data-level** (not header) observation on the same pair, reported for completeness and
*not* counted as a header divergence: GATK folds both payloads into the one surviving sample
(`gatk_combine.g.vcf` record at chr1:81 has `DP=60`, `PL=255,0,255`), while native keeps only the
first input's payload (its row has `DP=30`).

### N6 — Declaration ORDER of shared keys is identical  ✅ agree

For every case with a complete comparison (17 of the 18), the comparator reports
`shared-key relative ORDER identical (N shared keys)`: 17 primary comparisons against the final GATK
VCF plus 17 against the intermediate GVCF, i.e. 34 order comparisons, **0** reporting `DIFFERS`. GATK writes the header through htsjdk's sorted
header (`VCFHeader.getMetaDataInSortedOrder`, key/ID order) and native reproduces that order via
`gatk_htsjdk_sorted_header` (`genotype_gvcf_tool.cpp:4989-5021`, `std::sort` at 5006, applied at line 5209) plus the
per-group sorts at lines 5083-5100 (INFO), 5125 (FILTER), 5149-5152 (FORMAT). Concretely, the raw
filetext of native's `t2` output and GATK's final `t2` output agree on the position of every shared
line; only the missing `EXINFO` line differs. This matters for the fix proposal in §5: the ordering
machinery is already right, so a fix only has to fix the line **set**.

### N7 — Baseline: with identical headers, the final headers match exactly  ✅ agree

`t1_baseline` and `t12_disjoint_contigs`:

```
      [native vs FINAL gatk VCF header]
      native-only keys (0):
      gatk-final-only keys (0):
      same key, different text (0):
      shared-key relative ORDER identical (29 shared keys)
```

29 keys, no set difference, no text difference, no order difference — including the `MLEAC`/`MLEAF`
lines that are genuine GATK GenotypeGVCFs output annotations. (The comparator also prints a second
block against the intermediate `CombineGVCFs` GVCF; its "native-only `MLEAC`/`MLEAF`" lines are an
artifact of comparing a final VCF against an intermediate GVCF and are **not** a divergence. Only
the "[native vs FINAL gatk VCF header]" block is the verdict.)

### N8 — htslib's canonical header grouping happens to match htsjdk's sorted order  ✅ agree

Native's header line order is not the input file's line order (input: `fileformat, source, contig,
contig, FILTER, FORMAT…, INFO…, ALT`; native output: `fileformat, ALT, FILTER, FORMAT…, INFO…,
contig…, source`). That reordering is htslib's, and it coincides with htsjdk's key-sorted order, so
it is a non-divergence here. I state the mechanism as *not established* (I did not trace htslib
1.22.1's header writer); only the equality of the resulting order is measured.

---

## 4. Row-level (c) differences observed — and which of them are NOT header artifacts

Stated separately because it would be wrong to attribute them to the merge.

Present in **every** case, including the byte-identical-header baseline `t1_baseline`:

| axis | native | GATK |
|---|---|---|
| `INFO/DP` on chr1:81 | `DP=30` | `DP=60` (2 inputs), `DP=90` (3 inputs) |
| per-record FORMAT keys | `GT:DP:GQ:PL:MIN_DP` | `GT:AD:DP:GQ:PL` |
| missing sample slot | `./.:.:.:.,.,.:.` | `./.` (see D7) |

Since these appear when the two inputs' headers are identical, they are **independent** divergences
(summed-depth accounting, FORMAT key set/`AD` emission, missing-sample rendering), not consequences
of first-input-only header assembly. D7 is the only one of these three that I was asked about
directly; the other two are recorded as observations and out of this round's scope.

Differences that **are** direct consequences of the header divergence (absent in the baseline,
present exactly where the header sets differ):

* `t2`, `t10`: GATK's INFO contains `EXINFO=5`, native's does not (native baseline order).
* `t3`, `t10`: GATK's FORMAT column contains `EXFMT`, native's does not.
* `t15`: native emits `<NON_REF>` with no `##ALT` declaration.

---

## 5. Minimal fix site and what a merge would have to do (proposal — NOT implemented)

### 5.1 Exact sites

**Default (non-streaming) path — `run_tool`, `fastgatk-native/src/genotype_gvcf_tool.cpp`:**

* **6544-6546** — `if (!output_header) { output_header = bcf_hdr_dup(header); … }`
  This single `bcf_hdr_dup` of the first input's header is the whole defect: it is the only place the
  output header's declaration content is decided.
* **6547-6598** — the inlined `add_genotype_output_header_fields` equivalent that appends GATK's own
  annotations (guard: `bcf_hdr_id2int(output_header, BCF_DT_ID, id) < 0`).
* **6604** — `merge_contig_dictionary(output_header, header);` — the correct, already-existing
  precedent: called for **every** input, unions by contig name, throws on length conflict.
* **6605-6613** — sample union via `bcf_hdr_add_sample`.
* **6614** — `bcf_hdr_sync`.

**Streaming path — `run_streaming_genotype_gvcf`, same file (not on the measured code path, but the
same bug):**

* **5871-5879** — `output_header = bcf_hdr_dup(header)` / `add_genotype_output_header_fields()`
  (that helper is defined at **5655-5718**).
* **5880** — `merge_contig_dictionary(...)`.
* **5901-5909** — sample union; **5914** — `bcf_hdr_sync`.

**Where the header is written (so the fix must be upstream of it):**

* non-streaming: `7017-7027` `bcf_hdr_format(output_header, …)` →
  `gatk_compatible_header_text(...)`; `7029` `bcf_hdr_write(output, output_header)`. The writer is
  created at **7014**.
* streaming: `6192-6204` / `6205`, writer created at **6187**.
* `gatk_compatible_header_text` is defined at **5023** and `gatk_htsjdk_sorted_header` at **4989**.

### 5.2 What a merge has to do

Mirror `htsjdk.variant.vcf.VCFUtils.smartMergeHeaders(Collection<VCFHeader>, true)` exactly, per
`MultiVariantDataSource.getMergedHeader()` (`MultiVariantDataSource.java:249-260`). The algorithm,
read off the bundled bytecode of the pinned jar:

```
map = LinkedHashMap<String, VCFHeaderLine>       // insertion-ordered
for each header h, in input order:
  for each line l in h.getMetaDataInSortedOrder():
     key = l.getKey()
     if l instanceof VCFHeaderIDLine: key = key + "=" + l.getID()   // "INFO=EXINFO", "contig=chr1"
     if (!map.containsKey(key)):  map.put(key, l);  continue        // union: new IDs are ADDED
     prev = map.get(key)
     if l.equals(prev):           continue                          // duplicate: drop
     if l.getClass() != prev.getClass():  throw IllegalStateException
     if l instanceof VCFFilterHeaderLine:                            // FILTER
         if l.getID() != prev.getID(): throw;  else continue
     if l instanceof VCFCompoundHeaderLine:                          // INFO/FORMAT
         if !l.equalsExcludingDescription(prev):
             same type       -> warn(prev), prev.setNumberToUnbounded()
             Integer vs Float-> warn, keep prev
             Float vs Integer-> warn, keep prev
             else            -> throw IllegalStateException
         if !l.description.equals(prev.description): warn
         continue                                                    // first description wins
     // plain line (source, fileformat, fileDate, ...): FIRST VALUE WINS, later dropped + warn
```

The required consequences for native, in the order they matter:

1. **Sample columns must be `SortedSet` order.** `VcfUtils.getSortedSampleSet(map, REQUIRE_UNIQUE)`
   returns `new TreeSet<>()` of the raw sample names
   (`GATKVariantContextUtils.mergedSampleName(name, sample, uniquify) == uniquify ? sample+"."+name : sample`,
   and `REQUIRE_UNIQUE` passes `uniquify == false`). So native needs to sort `sample_names`
   (`genotype_gvcf_tool.cpp:5839`, `6492`) instead of numbering them by first appearance at
   6533-6541 / 6526-6542, and to drop the file-name suffix behaviour it does not have. Duplicate
   names must be silently collapsed (measured GATK behaviour), which native already does.
2. **Declaration union by `(group, ID)`, keyed exactly as above.** A new helper —
   e.g. `void merge_declaration_dictionary(bcf_hdr_t* output_header, const bcf_hdr_t* input_header)`
   — called next to `merge_contig_dictionary` at **6604** and **5880** (and folded into
   `merge_contig_dictionary`'s role), must append every `##INFO` / `##FORMAT` / `##FILTER` / `##ALT`
   line from `input_header` whose `key=ID` is not already declared in `output_header`. htslib
   provides the primitives in the bundled htslib 1.22.1
   (`third_party/htslib-build/htslib-src/htslib/vcf.h:693` `bcf_hrec_format`,
   `:596` `bcf_hdr_append`, `:712` `bcf_hdr_add_hrec`), so this is an iterate-`hrec` +
   `bcf_hrec_format` + `bcf_hdr_append` loop, plus `bcf_hdr_id2int(output_header, BCF_DT_ID, id) < 0`
   as the "already present" test — the same idiom already used at 6549-6598. `merge_contig_dictionary`
   (995-1032) is the existing template for the iteration and error handling.
3. **Non-ID lines keep the first input's value.** `##source`, `##fileformat`, `##fileDate`,
   `##reference`-style lines must **not** be unioned; native's current behaviour (keep input 1) is
   already the GATK behaviour (N1), so the fix must not "improve" this. The one caveat is the
   `##fileformat` version: `smartMergeHeaders` calls `enforceHeaderVersionMergePolicy`, which throws
   `IllegalArgumentException` when the accumulated version set has `size > 1` and contains `VCF4_3`.
   In my `t18` probe this did not fire, because the VCFv4.2 header reported
   `getVCFHeaderVersion() == null` (measured; the reason is a parser field I did not trace — see §6).
   Reproducing the throw would require reproducing that version bookkeeping too; I did not test a
   case that makes the policy fire.
4. **Ordering needs no change**, because the emitted header is re-sorted:
   `gatk_compatible_header_text` (5023) ends with
   `output = gatk_htsjdk_sorted_header(std::move(output))` (5209), and
   `gatk_htsjdk_sorted_header` (4989-5021, `std::sort` at 5006) does `std::sort` over all non-contig lines by their full
   text, which is htsjdk's `VCFHeaderLine.compareTo`. Verified against measurement: GATK's final `t2`
   header places `EXINFO` between `END` and `ExcessHet`, and native's sort would place the same line
   in the same position (`"##INFO=<ID=EXINFO,"` < `"##INFO=<ID=ExcessHet,"` since `'X' < 'x'`).
   The candidate-derived `info_order` ranking at 5070-5100 is *not* consulted for lines it does not
   know (they all get `rank == info_order.size()`), but it does not matter for the final output
   because of the 5209 re-sort. If the 5209 re-sort is ever bypassed, this becomes a second fix site.
5. **Contig handling needs no change**: N2/N3/N4 show union + first-order + length-conflict rejection
   already match GATK.
6. **The missing-sample row (D7)** is *not* a header fix: the header is what makes the FORMAT keys
   exist, but the choice of `./.` versus `./.:.:.:.,.,.:.` is made in the record encoder
   (`apply_genotype_assignment`, 1774+, and the writer's per-sample emission). Native must emit a
   single token for a sample with no input data at the locus. Flagging it here because it was asked
   and because it is the one header-adjacent divergence that a header merge alone will not fix.

### 5.3 Source references for the required behaviour

| behaviour | reference |
|---|---|
| GenotypeGVCFs accepts exactly one `-V` | `gatk-source/.../engine/VariantLocusWalker.java:33-35`; `gatk-source/.../tools/walkers/GenotypeGVCFs.java:65,99,112` |
| merged header construction | `gatk-source/.../engine/MultiVariantDataSource.java:249-260` |
| sorted/unique sample set | `MultiVariantDataSource.java:234-243`; `VcfUtils.getSortedSampleSet` (bytecode of pinned jar) |
| cross-dictionary validation (contig length) | `MultiVariantDataSource.java:277-339` |
| header-merge algorithm | `htsjdk.variant.vcf.VCFUtils.smartMergeHeaders(Collection<VCFHeader>, boolean)` — not vendored in `gatk-source`; read from `javap -p -c` of `third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar`, corroborated against htsjdk 4.2.0 (`HTSJDK Version: 4.2.0` per the GATK banner; disassembly identical modulo constant-pool/`ldc_w` artifacts) |
| htslib primitives for the native merge | `third_party/htslib-build/htslib-src/htslib/vcf.h:596,693,712` |

---

## 6. What I could not falsify, and what remains unverified

* **The `-V`-swap direction is symmetric**, but I verified it only for the INFO axis (`t2`). I did not
  separately swap `-V` for the FORMAT, FILTER and ALT cases; the mechanism (a single `bcf_hdr_dup`
  of `input_paths.front()`) is identical for all of them and the comparator already shows the
  symmetric loss, but the swap experiment itself is one case.
* **`smartMergeHeaders` conflict branches were never exercised end-to-end.** My `Number`/`Type`
  conflicts (e.g. `##INFO=<ID=X,Number=1,...>` vs `Number=.` or `Integer` vs `Float`) were not
  built as fixtures, so the `setNumberToUnbounded` / Integer-vs-Float warn-and-keep behaviour is
  read from bytecode, not measured from a GATK output. The `IllegalStateException` branches likewise.
* **The `##fileformat` version policy is unresolved.** In `t18` (`VCFv4.2` + `VCFv4.3`) the expected
  `IllegalArgumentException` did not fire; the probe showed `getVCFHeaderVersion() == null` for the
  v4.2 header, which makes `enforceHeaderVersionMergePolicy` return early (`ifnull 100`). I did not
  establish why the v4.2 header's version field is null. Reported as a measured observation with
  low confidence, not as a claim about GATK.
* **The htslib header reordering in native is not traced.** I measured that native's output order
  equals htsjdk's sorted order; I did not read htslib 1.22.1's `bcf_hdr_format`/writer to establish
  *why*, so I cannot say whether the agreement is guaranteed or coincidental. Treat N8 as
  "measured to agree on these 17 cases", not "structurally identical".
* **`--gatk-compatible-annotations` off.** All comparisons above used the flag, as specified.
  Without it native writes via `bcf_hdr_write` (7029) and *does* emit
  `##fastgatk_genotype_gvcfs_status=reference-block-materialization`, which
  `gatk_compatible_header_text` strips at line 5046. I did not run the flagless profile across the
  matrix, so I make no claim about it beyond that one measurement.
* **GenomicsDB path not exercised.** `GenomicsDBImport -V a -V b` + `GenotypeGVCFs -V gendb://…` is
  the third GATK multi-input route and the one native's `gendb://` support actually mirrors. I used
  `CombineGVCFs`, which shares the same `MultiVariantDataSource` header code, and did not run the
  GenomicsDB route. If a GenomicsDB workspace's merged header differs from `CombineGVCFs`'s, D2-D5
  would need re-measuring against it.
* **No GATK counterpart for D8.** Because `GenotypeGVCFs -V a -V b` is a hard error, "native differs
  from GenotypeGVCFs" on the multi-input axis is not a statement about wrong output; the comparison
  above is native vs the `CombineGVCFs`-then-`GenotypeGVCFs` pipeline, which is the closest
  reachable equivalent and the workflow GATK's own documentation points users at.
* **Only `chr1:81` was genotyped as a concrete variant.** Each pair has exactly one variant locus,
  so the (c) row comparisons are one row per case. FORMAT/INFO emission issues that only appear at
  reference-block-only loci (e.g. `t5`/`t17`, where GATK's combined GVCF emits `./.` with data)
  are only lightly covered.

---

## 7. Reproduction

```
# 1. build fixtures (deterministic; nothing under fastgatk-native/src is touched)
python3 .diag/round49/make_inputs.py
python3 .diag/round49/make_extra_inputs.py

# 2. run the 12-case sweep and the 6 follow-up cases
bash .diag/round49/run_all.sh
bash .diag/round49/run_extra.sh          # CASES= in that script was edited per batch; see logs

# 3. compare (writes the raw log shipped with this report)
python3 .diag/round49/compare.py > fastgatk-native/evidence/2026-09-13-round49/compare_all.txt
```

Probe (GATK's own merge helpers, run in-process against the pinned jar):

```
third_party/jdk17/bin/javac -cp third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
    -d .diag/round49/probe_classes .diag/round49/Probe.java
third_party/jdk17/bin/java -cp third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar:.diag/round49/probe_classes \
    Probe .diag/round49/in/t5_collision/in1.g.vcf .diag/round49/in/t5_collision/in2.g.vcf
```

Every command above was executed in this round; every header/row string quoted in this report is
copied from a file under `.diag/round49/out/<case>/`.

---

## 8. One-paragraph answer to the stated suspicion

**Partly confirmed, but not on the terms stated.** The suspicion is right that native builds its
output header from the first input only — proven by reading
`genotype_gvcf_tool.cpp:6544-6545` and, empirically, by swapping the two `-V` arguments and watching
`EXINFO` appear and `ONLY1INFO` vanish. It is also right that later inputs' `INFO`, `FORMAT`,
`FILTER` and `ALT` declarations are lost, and that this contradicts a proper merge. But the stated
GATK behaviour is wrong for the pinned oracle: GATK 4.6.2.0 `GenotypeGVCFs` cannot take more than one
`-V` at all (`A USER ERROR … Argument 'V/variant' cannot be specified more than once.`), so the
comparison had to be made against the reachable merge path (`CombineGVCFs`, i.e.
`MultiVariantDataSource` + `htsjdk.variant.vcf.VCFUtils.smartMergeHeaders`). Against that oracle,
three additional results came out that the suspicion did not predict: native **agrees** with GATK on
contig union, contig order, contig-length rejection, `##source` first-wins behaviour and
sample-collision collapsing; native **disagrees** on sample column **order** (input order vs
lexicographic sorted, invisible in every case whose names happened to sort into input order); and
for a sample with no input data at a locus native writes `./.:.:.:.,.,.:.` where GATK writes a bare
`./.`.
