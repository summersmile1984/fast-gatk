# R1 — Exact GATK semantics for a null-allele `--alleles` record, and its exit code

Author: R1 audit subagent. Date: 2026-09-11.
Scope: **read-only audit**. No production source, CMake file or root `*.md` was modified; nothing was
rebuilt; no `git commit`. All experiment data lived inside `tempfile.TemporaryDirectory`; the driver
scripts were heredoc'd into `/tmp` (which does not persist across bash invocations, so every
experiment ran inside one command). This report file is the only artifact created in the repo.

Oracle / binaries used (both build trees assumed current, untouched):

```
third_party/jdk17/bin/java -Xmx1g -jar third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar HaplotypeCaller
fastgatk-native/build/fastgatk-hc-call            (Kokkos/OpenMP)
```

Reused prior evidence (not re-derived): `.diag/track-b-alleles-findings.md` §D3 (the
`chr1 903 . CAGCAG CAG . PASS .` reproducer).

---

## Bottom line

1. `Event.makeMinimalRepresentation` (`Event.java:44-62`) removes **only a shared SUFFIX**, from
   **both** REF and ALT, with a count bounded by `min(|REF|,|ALT|)`; it **never** moves `start`
   (only `stop` shrinks). It is *not* a left-aligning/prefix-trimming routine.
2. It yields a null allele iff **the shorter allele is a whole suffix of the longer one** and
   neither allele has length 1. `Allele.create(new byte[0], …)` then throws
   `IllegalArgumentException: Null alleles are not supported` from htsjdk's `SimpleAllele`
   constructor — from `Event.java:60` when the **ALT** empties, from `Event.java:59` when the
   **REF** empties (both observed).
3. The throw happens **per region, lazily**: `HaplotypeCallerEngine.callRegion` → `line 761`
   `splitVariantContextToEvents` → `GATKVariantContextUtils.java:1777` → `Event.ofWithoutAttributes`
   → `Event.<init>` → `makeMinimalRepresentation`. It is reachable only because the **biallelic**
   fast path (`GATKVariantContextUtils.java:1719-1721`) returns the VC **untrimmed**; multiallelic
   records go through `trimAlleles` (`:1749`), whose *restore-one-base-at-end* rule
   (`:1469-1477`) prevents the empty allele — a multiallelic record with the identical null-looking
   pair is accepted (observed, case `c04`).
4. GATK exits **3**, produced by the catch-all `catch (final Exception e)` in
   `Main.mainEntry` (`Main.java:231-233`) → `ANY_OTHER_EXCEPTION_EXIT_VALUE = 3`
   (`Main.java:81`) → `handleNonUserException` = `printStackTrace()` (`Main.java:274-276`). It is a
   *deliberate constant* but a **bucket-level** one ("any uncaught non-user exception"), not a
   null-allele-specific code. Native's single error channel returns **2** (`hc_call.cpp:7006-7009`),
   so exit-code parity is a separate decision from abort-vs-continue.
5. Native currently **silently accepts** the record (exit 0, normal calls, e.g. 8 data rows on the
   chr20 fixture where GATK aborts). Minimal fix = reject at the `ForcedAllele` loading boundary
   with a non-zero exit; hunk in the last section (`hc_call.cpp`, **not applied**). Exact `3` would
   require a new exit-code channel — recommended against (see §5).

---

## 1. The trimming algorithm of `Event.makeMinimalRepresentation`

`gatk-source/src/main/java/org/broadinstitute/hellbender/utils/haplotype/Event.java:28-62` (verbatim,
incl. the `differentLastBase` helper at `:130-132`):

```java
28    public Event(final String contig, final int start, final Allele ref, final Allele alt) {
29        Utils.validateArg(ref.isReference(), "ref is not ref");
30        this.contig = contig;
31        this.start = start;
32        final Pair<Allele, Allele> minimalAlleles = makeMinimalRepresentation(ref, alt);
33        refAllele = minimalAlleles.getLeft();
34        altAllele = minimalAlleles.getRight();
35        stop = start + refAllele.length() - 1;
36    }
...
44    private static Pair<Allele, Allele> makeMinimalRepresentation(Allele ref, Allele alt) {
45        //check for minimal representation
46        if(ref.length() == 1 || alt.length() == 1 || differentLastBase(ref.getBases(), alt.getBases())) {
47            return new Pair<>(ref, alt);
48        }
49        Utils.validateArg(!ref.basesMatch(alt), "ref and alt alleles are identical");
50        final byte[] refBases = ref.getBases();
51        final byte[] altBases = alt.getBases();
52        int overlapCount = 0;
53        final int minLen = Math.min(refBases.length, altBases.length);
54        while (overlapCount < minLen && refBases[refBases.length - 1 - overlapCount] == altBases[altBases.length - 1 - overlapCount]) {
55            overlapCount++;
56        }
57        final byte[] newRefBases = Arrays.copyOf(refBases, ref.getBases().length - overlapCount);
58        final byte[] newAltBases = Arrays.copyOf(altBases, alt.getBases().length - overlapCount);
59        final Allele newRefAllele = Allele.create(newRefBases, true);
60        final Allele newAltAllele = Allele.create(newAltBases, false);
61        return(new Pair<>(newRefAllele, newAltAllele));
62    }
...
130    private static boolean differentLastBase(final byte[] ref, final byte[] alt) {
131        return ref.length == 0 || alt.length == 0 || ref[ref.length-1] != alt[alt.length-1];
132    }
```

Answers to the sub-questions, each grounded in the lines above:

* **Order / direction.** There is a single loop (`:54`) walking **right to left** counting equal
  bases; it removes a **shared suffix only**. There is no shared-*prefix* removal in `Event` at all
  (the javadoc at `:39` says "removes identical suffixes"). Prefix removal exists only in
  `GATKVariantContextUtils.trimAlleles` → `AlignmentUtils.normalizeAlleles`, and that path is
  **not** taken for a biallelic `--alleles` record (§3).
* **Both alleles?** Yes — the same `overlapCount` is subtracted from REF (`:57`) and ALT (`:58`);
  both are then re-created (`:59`, `:60`). A cheap early exit at `:46` skips everything when either
  allele has length 1 or the last bases differ.
* **Position.** `start` is stored unchanged at `:31` and is **never modified** by the trimming; only
  `stop` is recomputed at `:35` as `start + <trimmed REF length> - 1`. So the Event's interval
  shrinks from the right and its 5′ coordinate is preserved; no left-shifting/normalization.
* **Guard at `:49`.** `Utils.validateArg(!ref.basesMatch(alt), "ref and alt alleles are identical")`
  yields a *different* message (`ref and alt alleles are identical`). For the `--alleles` path it is
  unreachable: `REF == ALT` is rejected earlier by htsjdk while decoding (§2, case `s02`).
* The comment in `AssemblyResultSet.java:769` (`isEventPresentInAssembly`) confirms the intent:
  *"note that Events are forced to have a minimal representation"*.
* Maximum trim = `minLen`, i.e. the loop can reach `overlapCount == minLen` (`:54` is
  `overlapCount < minLen`), which is exactly why the shorter allele can vanish.

---

## 2. Precise precondition for a null allele, and edge cases

Let `R = REF.getBases()`, `A = ALT.getBases()` (as decoded by the VCF codec), `n = min(|R|,|A|)`.

`makeMinimalRepresentation` throws iff **all** of the following hold:

1. the record is **exactly biallelic** (one ALT). Reason: `GATKVariantContextUtils.java:1719-1721`
   returns the VC untouched for a biallelic record, whereas a multiallelic record's per-ALT
   biallelics go through `trimAlleles` (`:1749`), which *restores one base* when the trim would
   empty an allele (`:1469-1475`) — no empty allele can survive that path;
2. `|R| >= 2` **and** `|A| >= 2` (otherwise `:46` returns unchanged — the classic
   `REF=ACG ALT=A` deletion form can therefore never abort);
3. `R` and `A` share their last base (otherwise `differentLastBase` at `:46` returns unchanged);
4. `R` and `A` share their entire length-`n` suffix, i.e. **the shorter allele is a whole suffix of
   the longer one**. Then `overlapCount == n`, so `Arrays.copyOf(…, 0)` at `:57`/`:58` produces a
   zero-length array and `Allele.create(new byte[0], …)` throws.

Which line throws: **`:60`** if `|A| <= |R|` (ALT empties, `isRef=false`), **`:59`** if `|R| < |A|`
(REF empties, `isRef=true`). Both observed (below).

htsjdk side (htsjdk **4.2.0**, per `META-INF/MANIFEST.MF: htsjdk-Version: 4.2.0`; no htsjdk sources are
vendored, so this is bytecode + runtime-stack evidence, jar at
`~/.gradle/caches/modules-2/files-2.1/com.github.samtools/htsjdk/4.2.0/7f2d3ae…/htsjdk-4.2.0.jar`):

```
htsjdk.variant.variantcontext.SimpleAllele(byte[], boolean)   // javap -c
  24: aload_1 ; 25: invokestatic Allele.wouldBeNullAllele([B)Z
  28: ifeq 41
  31: new java/lang/IllegalArgumentException
  35: ldc  // String "Null alleles are not supported"
  40: athrow

htsjdk.variant.variantcontext.Allele.wouldBeNullAllele(byte[])  // javap -c
  return (bases.length == 1 && bases[0] == '-') || bases.length == 0;
```

`Allele.create(byte[], boolean)` reaches `new SimpleAllele(...)` directly for `length != 1`
(disassembly at `Allele.create`, offset 251), so a zero-length array always throws. Note the
`'-'` half of `wouldBeNullAllele`: a literal `-` is also a "null allele" for htsjdk, but for
`--alleles` a `-` ALT never reaches `Event` (the VCF codec rejects it during decode — case `e04`).

### Edge-case matrix (all observed on the pinned 4.6.2.0 jar; `OMP build/` binary for native)

chr1 = the track-b synthetic 1500 bp fixture (same generator as
`fastgatk-native/scripts/verify_hc_alleles_deep_boundary.py`: `random.seed(11)`, 10×`A` at 600,
8×`CAG` at 900), `-L chr1:800-1080`; chr20 = `fixtures/chr20/ref20mnp.fasta` + `mnp.bam`,
`-L 20:10019900-10020500`, reference at 10020228 = `GTAAGA`.

| # | `--alleles` record | why | GATK exit | GATK frame | native exit |
|---|---|---|---|---|---|
| s01 / c01 | `chr1:903 CAGCAG>CAG` | ALT = whole suffix → ALT empties | **3** | `makeMinimalRepresentation(Event.java:60)` | 0 |
| e02 (also run as `d02`) | `20:10020228 GTAAGA>AGA` | same, on the real fixture | **3** | `Event.java:60` | 0 |
| f02 | `20:10020228 GTAAG>AG` | ALT len 2 is a whole suffix | **3** | `Event.java:60` | 0 |
| s03 / d03 | `20:10020228 AGA>GTAAGA` | REF = whole suffix → REF empties | **3** | `Event.java:59` | 0 |
| s04 | `chr1:907 AG>CAG` | same, REF len 2 | **3** | `Event.java:59` | 0 |
| f01 | `20:10020228 GTAAGA>A` | ALT len 1 → early exit at `:46` | 0 | – | 0 |
| c14 | `chr1:904 A>GA` | REF len 1 → early exit at `:46` | 0 | – | 0 |
| c08 | `chr1:900 CAGCAGCAG>CAGC` | shared *prefix*, last bases differ | 0 (1 row) | – | 0 (0 rows)* |
| c16 | `chr1:903 CAGCAG>CAGCAT` | last bases differ → `:46` | 0 | – | 0 |
| s02 / c05 | `chr1:903 CAGCAG>CAGCAG` (`REF == ALT`) | htsjdk VCF decode rejects | **3** | `TribbleException: Duplicate allele added to VariantContext: CAGCAG` (`AbstractVCFCodec.parseVCFLine(463)`), **not** `Event` | 0 |
| e04 | `20:10020228 GTAAGA>-` | VCF decode rejects | **3** | `TribbleException: … unparsable vcf record with allele -` (`AbstractVCFCodec.checkAllele(678)`) | 0 |
| c09 | `chr1:903 CAGCAG><DEL>` | symbolic (last base `>`) → `:46` | 0 | – | 0 |
| e05 | `20:10020228 GTAAGA><NON_REF>` | symbolic → `:46` | 0 (7 rows) | – | 0 (7 rows) |
| c10 | `chr1:903 CAGCAG>*` | `*` len 1 → `:46` (and `AssemblyResultSet.isSymbolic`, `:772-774`, skips SPAN_DEL) | 0 | – | 0 |
| c15 | `chr1:903 CAGCAG>.` | NO_CALL, not a variant | 0 | – | 0 |
| c04 | `chr1:903 CAGCAG>CAG,TT` | **multiallelic** → `trimAlleles` restore-one-base → `CAGC/C` → `:46` | 0 | – | 0 |
| c11 | `chr1:903 CAGCAG>CAG` with `FILTER=LowQual` | `HaplotypeCallerEngine.java:760` drops it before `Event` | 0 | – | 0 |
| s05 / c12 | same + `--force-call-filtered-alleles` | filter bypassed | **3** | `Event.java:60` | 0 |
| c02 | s01 record, `-L chr1:500-780` | record never fetched into that region's FeatureContext | 0 | – | 0 |
| c03 | s01 record, `-L chr1:1300-1400` | region inactive → return before line 759 (`:754-757`) | 0 | – | 0 |
| e03 | `20:10020228 TTTTTT>TTT` (REF not in the FASTA) | reference is never validated for feature alleles | **3** | `Event.java:60` | 0 |
| e01 | control `20:10020228 G>A` (a real 4.6.2.0 call) | positive control | 0, 7 rows | – | 0, 6 rows* |

`*` = pre-existing, non-R1 divergences (see "Side observations").

So the *practical* detector for a native port is exactly: `|REF| >= 2 && |ALT| >= 2 && REF[|REF|-n:] == ALT[|ALT|-n:]`
with `n = min(|REF|,|ALT|)`, **for a record with exactly one ALT**.

Full unfiltered GATK stack trace (canonical `20:10020228 GTAAGA>AGA`, `-L 20:10019900-10020500`,
observed exit 3, output file present but header-only, 2636 bytes, 0 data rows):

```
java.lang.IllegalArgumentException: Null alleles are not supported
	at htsjdk.variant.variantcontext.SimpleAllele.<init>(SimpleAllele.java:56)
	at htsjdk.variant.variantcontext.Allele.create(Allele.java:191)
	at org.broadinstitute.hellbender.utils.haplotype.Event.makeMinimalRepresentation(Event.java:60)
	at org.broadinstitute.hellbender.utils.haplotype.Event.<init>(Event.java:32)
	at org.broadinstitute.hellbender.utils.haplotype.Event.ofWithoutAttributes(Event.java:66)
	at java.base/java.util.stream.ReferencePipeline$3$1.accept(ReferencePipeline.java:197)
	... (java.util.stream frames) ...
	at org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils.splitVariantContextToEvents(GATKVariantContextUtils.java:1777)
	at org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine.lambda$callRegion$9(HaplotypeCallerEngine.java:761)
	at org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine.callRegion(HaplotypeCallerEngine.java:763)
	at org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller.apply(HaplotypeCaller.java:284)
	at org.broadinstitute.hellbender.engine.AssemblyRegionWalker.processReadShard(AssemblyRegionWalker.java:200)
	at org.broadinstitute.hellbender.engine.AssemblyRegionWalker.traverse(AssemblyRegionWalker.java:173)
	at org.broadinstitute.hellbender.engine.GATKTool.doWork(GATKTool.java:1119)
	at org.broadinstitute.hellbender.cmdline.CommandLineProgram.runTool(CommandLineProgram.java:150)
	at org.broadinstitute.hellbender.cmdline.CommandLineProgram.instanceMainPostParseArgs(CommandLineProgram.java:203)
	at org.broadinstitute.hellbender.cmdline.CommandLineProgram.instanceMain(CommandLineProgram.java:222)
	at org.broadinstitute.hellbender.Main.runCommandLineProgram(Main.java:166)
	at org.broadinstitute.hellbender.Main.mainEntry(Main.java:209)
	at org.broadinstitute.hellbender.Main.main(Main.java:306)
```

---

## 3. Where the call sits in the `--alleles` loading path, vs. the native loader

### GATK (HaplotypeCaller)

```
HaplotypeCaller.apply(HaplotypeCaller.java:284)
  └─ AssemblyRegionWalker.processReadShard(:200)  / traverse(:173)      <- one region at a time
      └─ HaplotypeCallerEngine.callRegion(:680…)
           :754  if( ! region.isActive() ) return referenceModelForNoVariation(region, true, VCpriors);
           :759  final List<Event> givenAlleles = features.getValues(hcArgs.alleles).stream()
           :760          .filter(vc -> hcArgs.forceCallFiltered || vc.isNotFiltered())
           :761          .flatMap(vc -> GATKVariantContextUtils
           :761                  .splitVariantContextToEvents(vc, false, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, false).stream())
           :762          .filter(event -> event.getStart() >= region.getSpan().getStart())
           :763          .collect(Collectors.toList());
```

`splitVariantContextToEvents` (`GATKVariantContextUtils.java:1773-1778`):

```java
1776        return splitVariantContextToBiallelics(vc, trimLeft, genotypeAssignmentMethod, keepOriginalChrCounts).stream()
1777                .map(Event::ofWithoutAttributes).collect(Collectors.toList());
```

and `splitVariantContextToBiallelics` (`:1714-1755`) — **the decisive detail**:

```java
1717        if (!vc.isVariant()) {
1718            return Collections.emptyList();
1719        } else if (vc.isBiallelic())
1720            // non variant or biallelics already satisfy the contract
1721            return Collections.singletonList(vc);          // <- NO trimAlleles, VC unchanged
1722        else {
...
1749                final VariantContext trimmed = trimAlleles(builder.make(), trimLeft, true);  // <- multiallelic path only
```

`trimAlleles` (public overload, `:1455-1478`) is what saves multiallelic records:

```java
1469        final boolean emptyAllele = ranges.stream().anyMatch(r -> r.size() == 0);
1470        final boolean restoreOneBaseAtEnd = emptyAllele && startTrim == 0;
1471        final boolean restoreOneBaseAtStart = emptyAllele && startTrim > 0;
1473        // if the end trimming consumed all the bases, leave one base
1474        final int endBasesToClip = restoreOneBaseAtEnd ? endTrim - 1 : endTrim;
```
with `startTrim == 0` always here because `trimLeft=false` (`HaplotypeCallerEngine.java:761`) and
`AlignmentUtils.normalizeAlleles(sequences, ranges, 0, true)` is called with `maxShift = 0`
(`GATKVariantContextUtils.java:1465`, `AlignmentUtils.java:818-849`). Hence, for any multiallelic
record, one base is restored and the resulting (short) allele hits the `length()==1` guard in
`Event` — no abort. This is why the abort is specific to **exactly-biallelic** records.

Equivalent code exists in `RampedHaplotypeCallerEngine.java:248-253` and, for Mutect2,
`Mutect2Engine.java:266-269` (same `splitVariantContextToEvents(vc, false, …)` call, with
`MTAC.forceCallFiltered` as the filter gate and the inactive-region early return at `:260`).
**Not verified empirically for Mutect2** (source-level only).

### Native (`fastgatk-native/src/hc_call.cpp`)

`load_forced_alleles`, `hc_call.cpp:2213-2293` (read, not edited), called once from
`hc_call.cpp:6815-6816`:

```cpp
2213 std::vector<fastgatk::calling::ForcedAllele> load_forced_alleles(
2214     const std::string& path, const fastgatk::io::HeaderSummary& read_header,
2215     const bool include_filtered) {
...
2237         while (true) {                                  // eager: the WHOLE --alleles file
2238             const auto status = bcf_read(input, header, record);
...
2257             if (!include_filtered && feature_record_is_filtered(header, record)) { ... continue; }   // == :760 filter
2261             if (!concrete_feature_allele(record->d.allele[0])) { ... continue; }                    // REF must be concrete
2265             const std::string reference(record->d.allele[0]);
2266             for (int alternate_index = 1; alternate_index < record->n_allele; ++alternate_index) {
2268                 if (!concrete_feature_allele(alternate)) continue;   // drops ".", "*", "<…>"
2270                 if (reference == alternate_text) continue;           // <-- REF==ALT: native SKIPS, GATK aborts
2274                 fastgatk::calling::ForcedAllele forced;              // raw strings, no minimalisation
2277                 forced.reference = reference;
2278                 forced.alternate = alternate_text;
```

with `concrete_feature_allele` at `:2207-2211` (`value != "." && value != "*" && front()!='<' && back()!='>'`).
Fields: `fastgatk-native/include/fastgatk/calling/pipeline.hpp:43-53` (`tid`, `position` 0-based
anchor, `reference`, `alternate`). Consumed per region by
`calling_pipeline.cpp:16079` → `inject_hc_given_alleles_into_haplotype_paths` (`:12797`), which
iterates `options.forced_alleles` and injects raw REF/ALT as an indel/SNP on base haplotypes
(skipping `reference == alternate` again at `:12875`).

### Structural comparison

| aspect | GATK | native |
|---|---|---|
| when records are read | lazily, per region (`features.getValues`, `:759`) | eagerly, once, whole file (`:2237`) |
| scope of a bad record | only regions that actually fetch it **and** are active (`:754`) | every record in the file, regardless of `-L` |
| normalisation before use | biallelic → **none**; multiallelic → `trimAlleles` (right-trim, restore-1-base) | none |
| minimalisation of REF/ALT | `Event.<init>` → `makeMinimalRepresentation` (`Event.java:32,44-62`) | **absent** |
| empty allele | `IllegalArgumentException` → exit 3 | impossible (raw strings are passed through) |
| `REF == ALT` | htsjdk VCF-decode `TribbleException` → exit 3 | skipped (`:2270`, and `:12875`) |
| `'-'` ALT | htsjdk decode error → exit 3 | kept as a literal 1-char ALT (`concrete_feature_allele` accepts it) |

The equivalent native boundary for an `Event`-level check is therefore **inside
`load_forced_alleles`, immediately after `reference`/`alternate_text` are materialised
(`hc_call.cpp:2265-2270`)** — with the caveat that native's loader is file-wide rather than
region-scoped (§5).

---

## 4. The GATK process exit status

`gatk-source/src/main/java/org/broadinstitute/hellbender/Main.java`:

```java
66    private static final int COMMANDLINE_EXCEPTION_EXIT_VALUE = 1;
71    public static final int USER_EXCEPTION_EXIT_VALUE = 2;
81    private static final int ANY_OTHER_EXCEPTION_EXIT_VALUE = 3;
86    private static final int OUT_OF_MEMORY_EXIT_VALUE = 137;
76    public static final int PICARD_TOOL_EXCEPTION = 4;
```

```java
212        } catch (final CommandLineException e){
217            System.exit(COMMANDLINE_EXCEPTION_EXIT_VALUE);
218        } catch (final PicardNonZeroExitException e) {
221            System.exit(PICARD_TOOL_EXCEPTION);
222        } catch (final UserException e){
224            System.exit(USER_EXCEPTION_EXIT_VALUE);
225        } catch (final StorageException e) {
227            System.exit(ANY_OTHER_EXCEPTION_EXIT_VALUE);
228        } catch (final OutOfMemoryError e) {
230            System.exit(OUT_OF_MEMORY_EXIT_VALUE);
231        } catch (final Exception e){
232            handleNonUserException(e);
233            System.exit(ANY_OTHER_EXCEPTION_EXIT_VALUE);
234        }
274    protected void handleNonUserException(final Exception exception) {
275        exception.printStackTrace();
276    }
```

`IllegalArgumentException` is a plain `Exception` (not a `UserException`, `CommandLineException`,
`StorageException` or `PicardNonZeroExitException`), so it lands in the `:231-233` bucket: **the stack
trace on stderr (observed above) and `System.exit(3)`**. `Main.main(String[])` (`:305-307`) calls
`mainEntry`; the jar's `META-INF/MANIFEST.MF` sets `Main-Class: org.broadinstitute.hellbender.Main`,
and the `gatk` wrapper execs `java`, so 3 propagates unchanged to the shell. Observed exit code in
**every** failing case above: `3` (c01, c05, c06, c07, c12, s01–s05, d03, d04, e02, e03, e04, f02).

**Stable or incidental?** Both, at different levels:

* *Contractual in form*: GATK deliberately owns exit-code policy in one place (`Main.java:66-86`;
  "this is the only method that is allowed to call System.exit", `:202`), and 3 is a named constant,
  so it does not change run to run and is what any harness will see.
* *Incidental in meaning*: 3 is not "null allele" — it is the catch-all for **any** uncaught
  non-user exception (a `NullPointerException` deep inside the assembly would also give 3). Nothing
  in GATK ties this error class to 3, so it is a bucket membership, not an error-specific mapping.
  Two other *different* failures with the same `--alleles` shape also exit 3 (htsjdk duplicate-allele
  and unparsable-allele `TribbleException`s, cases `s02`/`e04`).
* A partial output file is **not** part of the contract: GATK streams VCF rows as regions complete,
  so an abort can leave a header-only file (observed: 2636 bytes, 0 data rows) or a file with a few
  rows already written (observed: 2 data rows, same interval, different run). Do not assert on it.

---

## 5. What native would have to do, and whether matching is desirable

**Mechanics.** Native has exactly one error channel: any thrown `std::exception` reaches
`hc_call.cpp:7006-7009`, prints `error: <what>` on stderr and `return 2;` (the process's `main`
returns 0 on success only). Native never returns 3 (no `return 3` / `exit(3)` in `hc_call.cpp`).
So a faithful abort is a 1-line-ish change with exit code **2**, and exit code **3** would require
new plumbing (e.g. a dedicated exception type mapped in `main`, or a documented `--parity-exit-codes`
mode).

**Matching exit 3 literally — for.** A differential harness that only compares exit codes would go
green; any downstream pipeline that treats "3" as "GATK crashed on input" behaves identically.

**Matching exit 3 literally — against.** Native's exit codes are its own documented contract (`2` for
`BAD_INPUT`), and 3 in GATK is not null-allele-specific; encoding another tool's catch-all constant
adds a special case with no semantic content, and `error: …` on stderr already carries the diagnosis.
Recommendation: **abort with non-zero (2) + an explicit message; record the 3-vs-2 difference as a
documented divergence**; only wire 3 if a shipped oracle/regression test asserts the exact code.

**Abort vs. continue — for aborting.** (a) It is 1:1 by construction: GATK never accepts the record,
so "native accepts it" can never be validated against the oracle for such inputs; the only possible
agreements are coincidences of zero-row output. (b) The silent path is not harmless: on the chr20
fixture native exits 0 and writes **8 data rows** while GATK writes nothing — a caller that
implements a vetted `--alleles` list would silently get a different site set. (c) A `--alleles` list
with such a record is a user error that GATK deliberately surfaces at load time; reproducing it is
the point of the project.

**Abort vs. continue — against aborting (and the mitigations).** (a) The abort is *region-scoped* in
GATK: a null record outside the called intervals, or inside a region that is not active, is simply
never materialised as an `Event` and the run succeeds (cases `c02`, `c03`, both exit 0 on both
sides). Native's loader reads the whole file eagerly, so a naive loader-level check would newly
diverge on exactly those inputs (native 2 vs GATK 0) — a *new* divergence introduced by the "parity"
patch. (b) Native's raw-REF/ALT injection is not `Event`-minimalised at all, so a large class of
records (multiallelic in particular, case `c04`) is already handled non-1:1 for other reasons;
sprinkling an `Event`-only precondition into a loader that does not otherwise model `Event` buys a
partial guarantee. (c) Accepting the record and simply *skipping* it (drop the offending ALT with a
warning) would keep the rest of the sites, but that is a third behaviour matching neither GATK nor
today's native.
   The `c02`/`c03` problem can be removed later by making the check region-aware — the calling
   interval is already known before the load (`set_calling_intervals`, `hc_call.cpp:2183-2194`, runs
   at `:6814`, one line before `load_forced_alleles` at `:6815`), so passing the interval bounds (or
   restricting the scan to records overlapping them) is a small follow-up, not a redesign.
   Even then, "record overlaps `-L`" is only an approximation of "the region that fetched it was
   active" (`HaplotypeCallerEngine.java:754`), which cannot be known at load time at all.

---

## Side observations (NOT the null-allele issue; flagged, not investigated)

* `c08` (`chr1:900 CAGCAGCAG>CAGC`): GATK exit 0 with **1** data row, native exit 0 with **0** rows —
  a data-row divergence outside the abort class (same family as the audit's D1/D2).
* `e01` (control `20:10020228 G>A`, a call the pinned GATK makes): GATK 7 rows including
  `20 10020228 . G A`; native 6 rows and **missing** that site (and `10020229`), while native *does*
  emit `10020228 G>A` when the `--alleles` file contains the null record instead (`e02`, 8 rows).
  Native's forced-allele injection therefore perturbs its own assembled calls; again a pre-existing
  `--alleles`-boundary divergence, not an R1 finding.
* Native's annotation values differ slightly on the chr20 fixture (e.g. `FS=1.680` vs `FS=0.000`,
  `BaseQRankSum=3.630` vs `3.635` at 10020228) — the known `--drop-alleles`-persistent O1 class.

---

## Minimal native patch (unified diff — **NOT applied**)

Against the current content of `fastgatk-native/src/hc_call.cpp`. Line numbers are the real ones
(`@@ -2263,6 +2263,26 @@`); context is the actual file text at 2263-2268.

```diff
--- a/fastgatk-native/src/hc_call.cpp
+++ b/fastgatk-native/src/hc_call.cpp
@@ -2263,6 +2263,26 @@ std::vector<fastgatk::calling::ForcedAllele> load_forced_alleles(
                 continue;
             }
             const std::string reference(record->d.allele[0]);
+            // GATK parity -- Event.makeMinimalRepresentation (Event.java:44-62) strips the
+            // shared suffix of REF and ALT.  When the shorter allele is a full suffix of the
+            // longer one the strip consumes it entirely, htsjdk's SimpleAllele constructor
+            // throws IllegalArgumentException("Null alleles are not supported") and GATK
+            // exits 3 (Main.java:81,231-233).  A biallelic --alleles record is NOT pre-trimmed
+            // (GATKVariantContextUtils.java:1719-1721); a multiallelic one is (line 1749),
+            // where trimAlleles restores one base instead of emptying an allele, so only
+            // exactly-biallelic records can abort.
+            if (record->n_allele == 2 && concrete_feature_allele(record->d.allele[1])) {
+                const std::string alternate_text(record->d.allele[1]);
+                const std::size_t shorter = std::min(reference.size(), alternate_text.size());
+                if (reference != alternate_text && shorter > 1 &&
+                    reference.compare(reference.size() - shorter, shorter, alternate_text,
+                                      alternate_text.size() - shorter, shorter) == 0)
+                    throw std::runtime_error(
+                        "BAD_INPUT: Null alleles are not supported: --alleles record " +
+                        std::string(contig) + ":" + std::to_string(record->pos + 1) + " " +
+                        reference + ">" + alternate_text +
+                        " minimises to an empty allele (Event.makeMinimalRepresentation)");
+            }
             for (int alternate_index = 1; alternate_index < record->n_allele; ++alternate_index) {
                 const auto* alternate = record->d.allele[alternate_index];
                 if (!concrete_feature_allele(alternate)) continue;
```

Notes: `<algorithm>` and `<string>` are already included (`hc_call.cpp:11` and `:32`; `std::min` is
already used at `:6811`, `std::to_string` at `:1189`), and `contig` / `record` / `reference` are in
scope. Resulting behaviour: stderr `error: BAD_INPUT: Null alleles are not supported: …` and exit
**2** (native's single error channel, `hc_call.cpp:7006-7009`; the only `return`s from this `main`
are `0` at `:6853/:6872/:7005` and `2` at `:7009` — there is no `return 3`/`exit(3)`).

Both hunks were validated **without touching the repository**: copied to a scratch dir, then
`git apply --check` in strict mode (both hunks individually and cumulatively, on a `git init` scratch
copy) and `patch -p1 --dry-run` both pass; the report's context lines were additionally byte-compared
against the current file. `patch -F0` (GNU patch 2.8) reports fuzz on the optional hunk alone — a GNU
`patch` quirk with this hunk shape (`git apply` strict accepts it); `git apply` is the authoritative
check.

Standalone verification of the predicate (compiled with `g++` in a scratch dir, not against the
repo) — all ten cases agree with the observed GATK outcomes:

```
OK   CAGCAG>CAG    -> true  (s01/c01 exit 3)      OK   GTAAGA>A  -> false (f01 exit 0)
OK   CAG>CAGCAG    -> true  (s03  exit 3)         OK   CAGCAG>CAGC -> false (c08 exit 0)
OK   AG>CAG        -> true  (s04  exit 3)         OK   CAGCAGCAG>CAGC -> false (c08 exit 0)
OK   GTAAG>AG      -> true  (f02  exit 3)         OK   AC>GC     -> false (MNP w/ shared suffix)
OK   GTAAGA>AGA    -> true  (e02 exit 3)          OK   GTAAGA>GTAAGA -> true (s02: GATK aborts at decode)
```

Optional, *related-but-separate* hunk (same divergence class: GATK exit 3 at **VCF decode**, native
silent). Include only if the goal is full `--alleles` malformed-record parity:

```diff
--- a/fastgatk-native/src/hc_call.cpp
+++ b/fastgatk-native/src/hc_call.cpp
@@ -2266,8 +2266,15 @@
             for (int alternate_index = 1; alternate_index < record->n_allele; ++alternate_index) {
                 const auto* alternate = record->d.allele[alternate_index];
                 if (!concrete_feature_allele(alternate)) continue;
                 const std::string alternate_text(alternate);
-                if (reference == alternate_text) continue;
+                // GATK parity -- htsjdk rejects a duplicate allele while decoding the VCF
+                // ("Duplicate allele added to VariantContext: <bases>", observed for REF == ALT)
+                // and GATK exits 3; native currently skips the record silently.
+                if (reference == alternate_text)
+                    throw std::runtime_error(
+                        "BAD_INPUT: duplicate --alleles allele " + std::string(contig) + ":" +
+                        std::to_string(record->pos + 1) + " " + reference + ">" +
+                        alternate_text);
                 const auto key = std::make_tuple(target->second, record->pos,
                                                  reference, alternate_text);
                 if (!seen.insert(key).second) continue;
```

---

## What remains uncertain

1. **Region-scoped abort.** The patch is file-eager, GATK is per-region (cases `c02`, `c03`). After
   the minimal patch, native would exit 2 where GATK exits 0 for a null record outside `-L` or in an
   inactive region. A region-aware variant was not written (needs the interval or per-region
   metadata in the loader); the exact GATK predicate is "the active region's FeatureContext fetched
   this record" (`HaplotypeCallerEngine.java:754,759`), which is not computable at load time.
2. **Multiallelic duplicates.** I verified `REF==ALT` for a *biallelic* record only. The htsjdk
   duplicate-allele rejection is a decode-level rule, so `REF=X ALT=X,Y` presumably also exits 3,
   but I did not run it.
3. **Mutect2.** `Mutect2Engine.java:266-269` uses the same call, but the boundary differs (no
   `event.getStart() >= region start` filter, different filter flag, different inactive check at
   `:260`); I did not run Mutect2, and the native Mutect2 loader path was not inspected.
4. **Partial output.** Observed 0 data rows (2636-byte header-only file) in one run and 2 data rows
   in another for the same interval/record; the cause (writer buffering vs. per-region flush) was not
   investigated. Treated as non-contractual.
5. **htsjdk source.** No htsjdk sources are vendored; the `SimpleAllele`/`wouldBeNullAllele`
   behaviour is evidenced by `javap -c` on `htsjdk-4.2.0.jar` plus the runtime stack trace naming
   `SimpleAllele.java:56`. Line numbers `Allele.java:191` / `SimpleAllele.java:56` are htsjdk's, not
   arm's-length verified against htsjdk git.
6. **`-` ALT.** GATK rejects it at decode (`e04`); native's `concrete_feature_allele` accepts it and
   `ForcedAllele` would carry a literal `-`. Not patched; flagged as the same "silent accept" class.

## Reproduction recipe (each block inside one shell invocation)

```bash
# chr20 fixture (fixtures/chr20 already has .fasta/.fai/.dict/.bam/.bai)
W=$(mktemp -d); printf '##fileformat=VCFv4.2\n##contig=<ID=20,length=63025520>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n20\t10020228\t.\tGTAAGA\tAGA\t.\tPASS\t.\n' > "$W/n.vcf"
third_party/htslib-build/htslib-src/bgzip -f "$W/n.vcf"
third_party/htslib-build/htslib-src/tabix -f -p vcf "$W/n.vcf.gz"
COMMON="-R fixtures/chr20/ref20mnp.fasta -I fixtures/chr20/mnp.bam -L 20:10019900-10020500 --alleles $W/n.vcf.gz --create-output-variant-index false --add-output-vcf-command-line false"
timeout 300 third_party/jdk17/bin/java -Xmx1g -jar third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar \
    HaplotypeCaller $COMMON -O "$W/g.vcf"; echo "GATK exit=$?"        # 3, "Null alleles are not supported"
timeout 300 fastgatk-native/build/fastgatk-hc-call $COMMON --threads 2 -O "$W/n2.vcf"; echo "NATIVE exit=$?"   # 0
```

For the track-b synthetic fixture, re-use
`fastgatk-native/scripts/verify_hc_alleles_deep_boundary.py --case tandem-repeat-one-unit-literal`
(same reference generator, `random.seed(11)`), which is the D3 reproducer
(`chr1 903 . CAGCAG CAG . PASS .`, `-L chr1:800-1080`).
