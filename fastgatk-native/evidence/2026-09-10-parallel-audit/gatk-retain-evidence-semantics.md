# GATK retainEvidence() semantics for the gVCF spanning-deletion record at `20:10020680`

Scope: read-only research. **No repository file was modified, and nothing was rebuilt.**
All scratch output went to a per-call `/tmp` directory. The only file created is this report.

Pinned reference: `gatk-source/` (GATK 4.6.2.0 source).
Pinned GATK binary used for the experiments below:
`third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar` with `third_party/jdk17/bin/java`.

Repro command used for every GATK experiment in this report (only `-L` start varies):

```
-HC  -R fixtures/chr20/ref20mnp.fasta -I fixtures/chr20/mnp.bam
     -L 20:<START>-10020710 --emit-ref-confidence GVCF --max-mnp-distance 1
     --native-pair-hmm-threads 1 --add-output-vcf-command-line false
     --debug-genotyper-output <dbg.txt> -O <out.g.vcf>
```

`--debug-genotyper-output` is a first-class GATK argument (`HaplotypeCallerArgumentCollection.java:53,226`,
initialized at `HaplotypeCallerEngine.java:233-235`); it is the authoritative instrumentation for this question.

---

## TL;DR

| Question | Answer |
|---|---|
| Q1 | The gVCF `*` row carries the VariantContext **already annotated by `HaplotypeCallerGenotypingEngine`**; `RAW_MQandDP` comes from `RMSMappingQuality.annotateRawData` and `FORMAT/SB` from `StrandBiasBySample.annotate` (enabled by `HaplotypeCallerEngine.filterReferenceConfidenceAnnotations`), both reading the **same `AlleleLikelihoods<GATKRead,Allele>`** that was used for genotyping. |
| Q2 | `retainEvidence(target::overlaps)` where the evidence is the **realigned-to-best-haplotype** read collection and `target` is `new SimpleInterval(mergedVC).expandWithinContig(2, dict)`. For this record that is `[10020678, 10020683]`. Population = reads whose **post-realignment** alignment interval overlaps that window. |
| Q3 | Yes — window-invariant. Measured on pinned GATK: **8 reads, identical read names, in both windows**, even though the AssemblyRegion partition differs (`20:10020490-10020710` vs `20:10020659-10020710`). 27 is an artefact of using raw BAM CIGAR coordinates instead of realigned coordinates. |
| Q4 | **(a) the context-mapped population is correct; (b) the geometric fallback is the defect.** The fallback is not merely window-dependent — it is *semantically wrong in every window* (it would give 27 in window A too, and it gives a population GATK never uses). |
| Q5 | None of the four listed candidates is correct as stated. See `RECOMMENDED PATCH`: give the symbolic sibling candidate a `candidate_reads` row so it owns the region context (no new PairHMM work, no new marginalization rows). |

---

## Q1 — Which read collection produces `RAW_MQandDP` and `FORMAT/SB` on this row?

### Q1.1 The gVCF row is the genotyping engine's annotated VariantContext, passed through unchanged

`ReferenceConfidenceModel.calculateRefConfidence` emits the already-annotated call verbatim for any
position that starts a variant call:

```java
// gatk-source/.../haplotypecaller/ReferenceConfidenceModel.java:241-249
final VariantContext overlappingSite = GATKVariantContextUtils.getOverlappingVariantContext(curPos, variantCalls);
final List<VariantContext> currentPriors = VCpriors.isEmpty() ? Collections.emptyList() : getMatchingPriors(curPos, overlappingSite, VCpriors);
if (overlappingSite != null && overlappingSite.getStart() == curPos.getStart()) {
    if (applyPriors) {
        results.add(PosteriorProbabilitiesUtils.calculatePosteriorProbs(overlappingSite, currentPriors,
                numRefSamplesForPrior, options));
    } else {
        results.add(overlappingSite);
    }
```

The `variantCalls` list is `calledHaplotypes.getCalls()` from `HaplotypeCallerEngine.java:1019-1022`, i.e. the
VariantContexts built by `HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods`. **No re-annotation happens
inside the reference-confidence model.** (In plain VCF mode the analogous line is
`HaplotypeCallerEngine.java:1031-1034`, which calls `RMSMappingQuality.finalizeRawMQ`; that path is *not* taken
in ERC mode.)

### Q1.2 Call path inside the genotyping engine (this is the crux)

```java
// gatk-source/.../haplotypecaller/HaplotypeCallerGenotypingEngine.java:192-199
AlleleLikelihoods<GATKRead, Allele> readAlleleLikelihoods = readLikelihoods.marginalize(alleleMapper);
final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
final SimpleInterval variantCallingRelevantOverlap = new SimpleInterval(mergedVC).expandWithinContig(hcArgs.informativeReadOverlapMargin, sequenceDictionary);

// We want to retain evidence that overlaps within its softclipping edges.
readAlleleLikelihoods.retainEvidence(read -> readQualifiesForGenotypingPredicate.test(read, variantCallingRelevantOverlap));

readAlleleLikelihoods.setVariantCallingSubsetUsed(variantCallingRelevantOverlap);
```

```java
// HaplotypeCallerGenotypingEngine.java:261-269
final GenotypesContext genotypes = calculateGLsForThisEvent(readAlleleLikelihoods, mergedVC, noCallAlleles, ref, loc - refLoc.getStart(), dragstrs);
final GenotypePriorCalculator gpc = resolveGenotypePriorCalculator(...);
final VariantContext call = calculateGenotypes(new VariantContextBuilder(mergedVC).genotypes(genotypes).make(), gpc, givenAlleles);
if( call != null ) {
    readAlleleLikelihoods = prepareReadAlleleLikelihoodsForAnnotation(readLikelihoods, perSampleFilteredReadList,
            emitReferenceConfidence, alleleMapper, readAlleleLikelihoods, call, variantCallingRelevantOverlap);

    VariantContext annotatedCall = makeAnnotatedCall(ref, refLoc, tracker, header, mergedVC,
            mergedAllelesListSizeBeforePossibleTrimming, readAlleleLikelihoods, call, annotationEngine, preFilteringAlleleLikelihoods);
```

```java
// HaplotypeCallerGenotypingEngine.java:584-608  (the annotation collection)
private AlleleLikelihoods<GATKRead, Allele> prepareReadAlleleLikelihoodsForAnnotation(
        final AlleleLikelihoods<GATKRead, Haplotype> readHaplotypeLikelihoods,
        final Map<String, List<GATKRead>> perSampleFilteredReadList,
        final boolean emitReferenceConfidence,
        final Map<Allele, List<Haplotype>> alleleMapper,
        final AlleleLikelihoods<GATKRead, Allele> readAlleleLikelihoodsForGenotyping,
        final VariantContext call,
        final SimpleInterval relevantReadsOverlap) {

    final AlleleLikelihoods<GATKRead, Allele> readAlleleLikelihoodsForAnnotations;

    // We can reuse for annotation the likelihood for genotyping as long as there is no contamination filtering
    // or the user want to use the contamination filtered set for annotations.
    // Otherwise (else part) we need to do it again.
    if (hcArgs.useFilteredReadMapForAnnotations || !configuration.isSampleContaminationPresent()) {
        readAlleleLikelihoodsForAnnotations = readAlleleLikelihoodsForGenotyping;
        // the input likelihoods are supposed to have been filtered to only overlapping reads so no need to
        // do it again.
    } else {
        readAlleleLikelihoodsForAnnotations = readHaplotypeLikelihoods.marginalize(alleleMapper);
        readAlleleLikelihoodsForAnnotations.retainEvidence(relevantReadsOverlap::overlaps);
```

For this run the first branch is taken: `useFilteredReadMapForAnnotations = false` by default
(`HaplotypeCallerArgumentCollection.java:251`) and no contamination fraction is configured
(`StandardCallerArgumentCollection.java:68`, `isSampleContaminationPresent()`), so
**the very same retained collection serves genotyping and annotation.**

`makeAnnotatedCall` hands it to the annotation engine:

```java
// HaplotypeCallerGenotypingEngine.java:532-539
static protected VariantContext makeAnnotatedCall(byte[] ref, SimpleInterval refLoc, FeatureContext tracker, SAMFileHeader header, VariantContext mergedVC, int mergedAllelesListSizeBeforePossibleTrimming, AlleleLikelihoods<GATKRead, Allele> readAlleleLikelihoods, VariantContext call, VariantAnnotatorEngine annotationEngine, final AlleleLikelihoods<GATKRead, Haplotype> preFilteringAlleleLikelihoods) {
    final SimpleInterval locus = new SimpleInterval(mergedVC);
    ...
    final VariantContext untrimmedResult =  annotationEngine.annotateContext(call, tracker, referenceContext, readAlleleLikelihoods, Optional.empty(), Optional.empty(), Optional.ofNullable(preFilteringAlleleLikelihoods), a -> true);
```

### Q1.3 `RAW_MQandDP` — class and method

* `RMSMappingQuality implements InfoFieldAnnotation, StandardAnnotation, ReducibleAnnotation`
  (`annotator/RMSMappingQuality.java:45`) — `StandardAnnotation` is part of HC's default annotation group
  (`HaplotypeCallerEngine.java:284-286`).
* In ERC mode the engine is constructed with `useRaw = true`:

```java
// tools/walkers/haplotypecaller/HaplotypeCaller.java:264-265
final VariantAnnotatorEngine variantAnnotatorEngine = new VariantAnnotatorEngine(makeVariantAnnotations(),
        hcArgs.dbsnp.dbsnp, hcArgs.comps,  hcArgs.emitReferenceConfidence != ReferenceConfidenceMode.NONE, false);
```
  so `VariantAnnotatorEngine.addInfoAnnotations` takes the raw branch:

```java
// annotator/VariantAnnotatorEngine.java:391-392
if (useRawAnnotations && annotationType instanceof ReducibleAnnotation) {
    annotationsFromCurrentType = ((ReducibleAnnotation) annotationType).annotateRawData(ref, newGenotypeAnnotatedVC, likelihoods);
```

* The annotation itself:

```java
// annotator/RMSMappingQuality.java:105-119
public Map<String, Object> annotateRawData(final ReferenceContext ref,
                                           final VariantContext vc,
                                           final AlleleLikelihoods<GATKRead, Allele> likelihoods){
    Utils.nonNull(vc);
    if (likelihoods == null || likelihoods.evidenceCount() == 0) {
        return Collections.emptyMap();
    }
    ...
    calculateRawData(vc, likelihoods, myData);
```

```java
// annotator/RMSMappingQuality.java:213-229
private void calculateRawData(final VariantContext vc,
                             final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                             final ReducibleAnnotationData rawAnnotations){
    //GATK3.5 had a double, but change this to an long for the tuple representation (square sum, read count)
    long squareSum = 0;
    long numReadsUsed = 0;
    for (int i = 0; i < likelihoods.numberOfSamples(); i++) {
        for (final GATKRead read : likelihoods.sampleEvidence(i)) {
            long mq = read.getMappingQuality();
            if (mq != QualityUtils.MAPPING_QUALITY_UNAVAILABLE) {
                squareSum += mq * mq;
                numReadsUsed++;
            }
        }
    }
    rawAnnotations.putAttribute(Allele.NO_CALL, Arrays.asList(squareSum, numReadsUsed));
}
```

**`likelihoods.sampleEvidence(i)` is exactly the post-`retainEvidence` main evidence list.** For this record
that is 8 reads, all MAPQ 60 → `28800,8`. The `8` is `likelihoods.evidenceCount()`, i.e. the *retained*
population and nothing else.

### Q1.4 `FORMAT/SB` — class and method

In ERC mode GATK forces the per-sample strand annotation on and removes the INFO-level strand tests:

```java
// tools/walkers/haplotypecaller/HaplotypeCallerEngine.java:264-279
public static Collection<Annotation> filterReferenceConfidenceAnnotations(Collection<Annotation> annotations) {
    logger.info("Tool is in reference confidence mode and the annotation, the following changes will be made ... 'StrandBiasBySample' will be enabled. 'ChromosomeCounts', 'FisherStrand', 'StrandOddsRatio' and 'QualByDepth' annotations have been disabled");
    // Override user preferences and add StrandBiasBySample
    if (!annotations.contains(new StrandBiasBySample())) {
        annotations.add(new StrandBiasBySample());
    }
    // Override user preferences and remove ChromosomeCounts, FisherStrand, StrandOddsRatio, and QualByDepth Annotations
    return annotations.stream()
            .filter(c -> (
                    c.getClass() != (ChromosomeCounts.class) &&
                    c.getClass() != (FisherStrand.class) &&
                    c.getClass() != (StrandOddsRatio.class) &&
                    c.getClass() != (QualByDepth.class))
            ).collect(Collectors.toList());
}
```

```java
// annotator/StrandBiasBySample.java:61-85
public void annotate(final ReferenceContext ref, final VariantContext vc, final Genotype g,
                     final GenotypeBuilder gb, final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
    ...
    if ( likelihoods == null || !g.isCalled() ) { ... return; }

    final int[][] table = FisherStrand.getContingencyTable(likelihoods, vc, 0, Arrays.asList(g.getSampleName()));

    gb.attribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY, getContingencyArray(table));
}
```

`FisherStrand` inherits the static from `StrandBiasTest`, which consumes the likelihood object and nothing else:

```java
// annotator/StrandBiasTest.java:186-205
private static int[][] getContingencyTable( final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                                           final Allele ref, final List<Allele> allAlts,
                                           final int minCount, final Collection<String> samples) {
    final int[][] table = new int[ARRAY_DIM][ARRAY_DIM];
    for (final String sample : samples) {
        final int[] sampleTable = new int[ARRAY_SIZE];
        likelihoods.bestAllelesBreakingTies(sample).stream()
                .filter(ba -> ba.isInformative())
                .forEach(ba -> updateTable(sampleTable, ba.allele, ba.evidence, ref, allAlts));
        if (passesMinimumThreshold(sampleTable, minCount)) {
            copyToMainTable(sampleTable, table);
        }
    }
    return table;
}
```

Genotype annotations are applied by `VariantAnnotatorEngine.annotateContext` with the same `likelihoods`
object before INFO annotations are computed (`annotator/VariantAnnotatorEngine.java:350-370`).

**Conclusion for Q1:** the row's `RAW_MQandDP` (INFO) and `SB` (FORMAT) are both functions of the *single*
`AlleleLikelihoods<GATKRead,Allele>` instance created at `HaplotypeCallerGenotypingEngine.java:192` and reduced
by `retainEvidence` at line 197. There is no second, wider read collection anywhere in this path.

---

## Q2 — The `retainEvidence(...)` rule, quoted, and what it means for a `*` record

### Q2.1 The predicate

```java
// HaplotypeCallerGenotypingEngine.java:194
final SimpleInterval variantCallingRelevantOverlap = new SimpleInterval(mergedVC).expandWithinContig(hcArgs.informativeReadOverlapMargin, sequenceDictionary);

// HaplotypeCallerGenotypingEngine.java:197
readAlleleLikelihoods.retainEvidence(read -> readQualifiesForGenotypingPredicate.test(read, variantCallingRelevantOverlap));
```

```java
// HaplotypeCallerGenotypingEngine.java:342-350
private BiPredicate<GATKRead, SimpleInterval> composeReadQualifiesForGenotypingPredicate(final HaplotypeCallerArgumentCollection hcArgs) {
    if (hcArgs.applyBQD || hcArgs.applyFRD) {
        return (read, target) -> softUnclippedReadOverlapsInterval(read, target);
    } else {
        // NOTE: we must make this comparison in target -> read order because occasionally realignment/assembly produces
        // reads that consume no reference bases and this can cause them to overlap adjacent
        return (read, target) -> target.overlaps(read);
    }
}
```

Neither `--apply-bqd` nor `--apply-frd` is passed → **`target.overlaps(read)`**, i.e. the *read's current*
aligned interval `[getStart(), getEnd()]` vs the target interval.

* Margin: `informativeReadOverlapMargin = 2` (`AssemblyBasedCallerArgumentCollection.java:214`).
* Interval expansion: `expandWithinContig(padding, dict)` = `[start-padding, end+padding]` clamped to the contig
  (`SimpleInterval.java:347-367`).
* Removal mechanics — the predicate is applied to the main evidence list, and to the already-disqualified list:

```java
// utils/genotyper/AlleleLikelihoods.java:1179-1196
public void retainEvidence(final Predicate<? super EVIDENCE> predicate) {
    Utils.nonNull(predicate);
    final int sampleCount = samples.numberOfSamples();

    for (int s = 0; s < sampleCount; s++) {
        // Remove evidence from the primary data
        final List<EVIDENCE> sampleEvidence = this.evidenceBySampleIndex.get(s);
        final int[] removeIndices = IntStream.range(0, sampleEvidence.size())
                .filter(i -> !predicate.test(sampleEvidence.get(i)))
                .toArray();
        removeEvidenceByIndex(s, removeIndices);

        // If applicable also apply the predicate to the filters
        final List<EVIDENCE> sampleFiltered = filteredEvidenceBySampleIndex.get(s).stream()
                .filter(predicate).collect(Collectors.toList());
        filteredEvidenceBySampleIndex.set(s, sampleFiltered);
    }
}
```

`removeEvidenceByIndex` deletes rows outright (`AlleleLikelihoods.java:1300-1325`); the "disqualified" counter
printed by the debugger is `filteredSampleEvidence(0).size()`.

### Q2.2 Which interval, for *this* record?

`mergedVC` at line 194 is the merged EventMap VC **after** `replaceSpanDels` and **before** `<NON_REF>` is added:

```java
// HaplotypeCallerGenotypingEngine.java:174-177
final List<VariantContext> eventsAtThisLocWithSpanDelsReplaced = replaceSpanDels(eventsAtThisLoc,
        Allele.create(ref[loc - refLoc.getStart()], true), loc);
VariantContext mergedVC = AssemblyBasedCallerUtils.makeMergedVariantContext(eventsAtThisLocWithSpanDelsReplaced);
```
```java
// HaplotypeCallerGenotypingEngine.java:384-397
static VariantContext replaceWithSpanDelVC(final VariantContext variantContext, final Allele refAllele, final int loc) {
    if (variantContext.getStart() == loc) {
        return variantContext;
    } else {
        VariantContextBuilder builder = new VariantContextBuilder(variantContext)
                .start(loc).stop(loc)
                .alleles(Arrays.asList(refAllele, Allele.SPAN_DEL))
                .genotypes(GenotypesContext.NO_GENOTYPES);
        return builder.make();
    }
}
```

The debugger shows the actual merged VC for this event:

```
Event at: [VC HC0 @ 20:10020680-10020681 Q. of type=MIXED alleles=[CA*, *, AT] attr={} GT=[] filters= with 8 reads and 20 disqualified
```

So `new SimpleInterval(mergedVC)` = `20:10020680-10020681`, expanded by 2 → **`[10020678, 10020683]`**.
Crucially the `*` allele is *inside* this merged VC: the spanning deletion does **not** shrink, widen or
otherwise special-case the retained interval. The population is a single one shared by REF, the concrete ALT
and `*`.

### Q2.3 The reads are the **realigned** reads, not the raw BAM alignments

```java
// HaplotypeCallerEngine.java:959-966
//Realign reads to their best haplotype.
final SWParameters readToHaplotypeSWParameters = hcArgs.getReadToHaplotypeSWParameters();
if (!(hcArgs.pileupDetectionArgs.generatePDHaplotypes && !hcArgs.pileupDetectionArgs.useDeterminedHaplotypesDespitePdhmmMode)) {
    final Map<GATKRead, GATKRead> readRealignments = AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype(subsettedReadLikelihoodsFinal, assemblyResult.getReferenceHaplotype(), assemblyResult.getPaddedReferenceLoc(), aligner, readToHaplotypeSWParameters);
    subsettedReadLikelihoodsFinal.changeEvidence(readRealignments);
}
```

```java
// AssemblyBasedCallerUtils.java:107-119
public static Map<GATKRead, GATKRead> realignReadsToTheirBestHaplotype(final AlleleLikelihoods<GATKRead, Haplotype> originalReadLikelihoods, ...) {
    final Collection<AlleleLikelihoods<GATKRead, Haplotype>.BestAllele> bestAlleles = originalReadLikelihoods.bestAllelesBreakingTies(HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY);
    final Map<GATKRead, GATKRead> result = new HashMap<>(bestAlleles.size());
    for (final AlleleLikelihoods<GATKRead, Haplotype>.BestAllele bestAllele : bestAlleles) {
        final GATKRead originalRead = bestAllele.evidence;
        final Haplotype bestHaplotype = bestAllele.allele;
        final boolean isInformative = bestAllele.isInformative();
        final GATKRead realignedRead = AlignmentUtils.createReadAlignedToRef(originalRead, bestHaplotype, refHaplotype, paddedReferenceLoc.getStart(), isInformative, aligner, readToHaplotypeSWParameters);
        result.put(originalRead, realignedRead);
    }
    return result;
}
```

`createReadAlignedToRef` **overwrites the read's start** and gives it a
read→haplotype→reference CIGAR:

```java
// utils/read/AlignmentUtils.java:103-113
final Cigar readToRefCigar = applyCigarToCigar(swCigar, haplotypeToRef);
final CigarBuilder.Result leftAlignedReadToRefCigarResult = leftAlignIndels(readToRefCigar, refHaplotype.getBases(), readMinusSoftClips.getBases(), readStartOnReferenceHaplotype);
final Cigar leftAlignedReadToRefCigar = leftAlignedReadToRefCigarResult.getCigar();
// it's possible that left-alignment shifted a deletion to the beginning of a read and removed it, shifting the first aligned base to the right
copiedRead.setPosition(copiedRead.getContig(), readStartOnReference + leftAlignedReadToRefCigarResult.getLeadingDeletionBasesRemoved());
...
copiedRead.setCigar(newCigar);
```

The collection passed to `assignGenotypeLikelihoods` is the post-`changeEvidence` object
(`HaplotypeCallerEngine.java:980-994`; `filterAlleles` defaults to `false` —
`AssemblyBasedCallerArgumentCollection.java:377` — so `subsettedReadLikelihoodsFinal` *is* the object mutated
above), and per-event `marginalize` copies the evidence list before `retainEvidence` touches it:

```java
// utils/genotyper/AlleleLikelihoods.java:734-744
final List<List<EVIDENCE>> newEvidenceBySampleIndex = new ArrayList<>(sampleCount);
for (int s = 0; s < sampleCount; s++) {
    newEvidenceBySampleIndex.add(new ArrayList<>(evidenceBySampleIndex.get(s)));
}
final AlleleLikelihoods<EVIDENCE, B> result = new AlleleLikelihoods<>(new IndexedAlleleList(newAlleles), samples, newEvidenceBySampleIndex, ...);
```

**Answer to Q2:** the retained population is *the reads of the candidate's AlleleLikelihoods collection —
after realignment to their best haplotypes — whose current alignment interval overlaps
`mergedVC ± informativeReadOverlapMargin`. It is **not** "reads that overlap the event in the original BAM",
and it is **not** "reads carrying a PairHMM likelihood row for that particular allele" (there is no such
per-allele read subset in GATK: all alleles of a merged event share one evidence list).

### Q2.4 Empirical confirmation of the exact rule (pinned GATK, `--debug-genotyper-output`)

Region owning the event, window A (`-L 20:10019901-10020710`), the 11 main-evidence reads with their
**post-realignment** CIGARs and coordinates (`HaplotypeCallerEngine.java:968-978` prints
`readLikelihoods.sampleEvidence(0)`), and `retainEvidence` target `[10020678, 10020683]`:

```
read 0: 20GAVAAXX100126:6:3:6181:93202   cigar: 41M60H  mapQ: 60 loc: [10020660-10020700] ... length:41
read 1: 20FUKAAXX100202:1:8:20035:107624 cigar: 39M62H  mapQ: 60 loc: [10020662-10020700] ... length:39
read 2: 20GAVAAXX100126:1:62:21171:66429 cigar: 38M63H  mapQ: 60 loc: [10020663-10020700] ... length:38
read 3: 20FUKAAXX100202:4:43:2630:81606  cigar: 38M63H  mapQ: 60 loc: [10020663-10020700] ... length:38
read 4: 20GAVAAXX100126:8:44:10223:188523 cigar: 29M72H  mapQ: 60 loc: [10020672-10020700] ... length:29
read 5: 20FUKAAXX100202:3:45:1441:112187 cigar: 28M73H  mapQ: 60 loc: [10020673-10020700] ... length:28
read 6: 20FUKAAXX100202:8:24:3734:131116 cigar: 20M81H  mapQ: 60 loc: [10020681-10020700] ... length:20
read 7: 20FUKAAXX100202:7:3:12109:51312  cigar: 18M83H  mapQ: 60 loc: [10020683-10020700] ... length:18
read 8: 20FUKAAXX100202:8:5:18882:61156  cigar: 16M85H  mapQ: 60 loc: [10020685-10020700] ... length:16
read 9: 20FUKAAXX100202:1:8:2485:178522  cigar: 14M87H  mapQ: 60 loc: [10020687-10020700] ... length:14
read 10: 20FUKAAXX100202:2:24:13465:45294 cigar: 11M90H  mapQ: 60 loc: [10020690-10020700] ... length:11
```

Reads 0–7 overlap `[10020678,10020683]`; reads 8–10 start after `10020683` and do not.
Reads 0–7 are **exactly** the 8 read names that appear in the event's read-allele matrix — and exactly the 8
counted by `RAW_MQandDP=28800,8`.

Two independent controls:

1. `20FUKAAXX100202:7:3:12109:51312` has **raw BAM coordinates `10020684-10020784`** (pysam over
   `fixtures/chr20/mnp.bam`), i.e. it does **not** overlap `[10020678,10020683]` by original CIGAR, yet it *is*
   retained. After realignment it sits at `10020683-10020700`. This proves the predicate sees the realigned read.
2. Reads `...:8:5:18882:61156` (raw `10020686-10020786`), `...:1:8:2485:178522` (raw `10020688-10020788`) and
   `...:2:24:13465:45294` (raw `10020691-10020791`) are **not** retained; their realigned starts (`10020685`,
   `10020687`, `10020690`) are still outside the target.

---

## Q3 — Is GATK's population window-invariant? (Yes)

### Q3.1 Source argument

The predicate is `variantCallingRelevantOverlap.overlaps(read)` where the interval derives only from
`mergedVC` (a pure function of the EventMap at that locus) and the read objects are the region's realigned
reads. There is no `-L`/window/tile term anywhere in it. `informativeReadOverlapMargin` is a constant 2.

### Q3.2 Measured on pinned GATK

| `-L 20:<start>-10020710` | GATK AssemblyRegion owning the event | GATK retained evidence for `20:10020680` | `RAW_MQandDP` | `SB` |
|---:|---|---|---:|---|
| 10019901 (window A) | `20:10020490-10020710` | `with 8 reads and 20 disqualified` | `28800,8` | `0,0,3,3` |
| 10020381 (window B) | `20:10020659-10020710` | `with 8 reads and 20 disqualified` | `28800,8` | `0,0,3,3` |

The two debug traces are identical for that region: the same 33 reads handed to PairHMM and the same 11
main-evidence reads with the same names/cigars/coordinates
(checked programmatically: the sorted read-name sets of the two regions are equal).
The retained 8 names are identical in both windows.

**Why 8, and not 27:** 27 is the number of non-duplicate reads in `fixtures/chr20/mnp.bam` whose *original*
alignment overlaps `[10020678,10020683]` (33 overlapping reads minus 6 duplicates). GATK never looks at that
set: it looks at the region's realigned reads, of which 11 remained in the primary evidence list after
read disqualification, and only 8 of those overlap the expanded VC interval. The window does not enter the
computation; the read population entering it is the same in both windows, even though the AssemblyRegion
partition is not.

> Note on AssemblyRegion partitioning: it genuinely differs between the two windows
> (`20:10020490-10020710` vs `20:10020659-10020710`), and even the assembled haplotype spans differ
> (`[10020390-10020810]` vs `[10020559-10020810]`). That is why the *native* implementation's read/context
> bookkeeping drifts with the window. GATK absorbs the drift because the retained evidence is defined by the
> read/VC geometry alone. So GATK's invariance here is not "the region never changes" — it is
> "the region change does not alter this predicate".

---

## Q4 — Which native population models GATK?

**(a) the context-mapped population (`context_mapping_evidence`, 8 reads in window A) is the population
GATK uses. (b) the geometric fallback (`read_overlaps_annotation_interval && annotation_read_survives`,
27 reads in window B) is the defect — and it is defective in *both* windows, not only the wrong one.**

Reasoning:

1. GATK's rule is "read's **current (realigned)** alignment interval overlaps `mergedVC ± 2`".
   Native's context path (`calling_pipeline.cpp:4144-4150`) mirrors it: it prefers the candidate-local
   realigned interval (`context->start < event_end && context->end > event_start`) and only falls back to the
   source CIGAR when no realigned interval exists. Its interval arithmetic
   (`event_start = position - margin`, `event_end = position + refLen + margin`, 0-based half-open at
   `calling_pipeline.cpp:4116-4123`) is the exact 0-based equivalent of GATK's `[start-2, stop+2]`.
2. The fallback path (`calling_pipeline.cpp:4295-4297`) uses `read_overlaps_annotation_interval`, which is
   built on the **original source-CIGAR** interval (`reads.positions[record]` / `io::reference_end`,
   `calling_pipeline.cpp:6709-6715`). That is precisely the population GATK *does not* use; on this fixture it
   is 27 reads.
3. The fallback is not "wider but harmless": it changes the population that `RMSMappingQuality`-equivalent and
   `StrandBiasBySample`-equivalent annotations are computed over, which is why window B reports `97200,27`
   and `0,0,0,0`. In window A the fallback would produce the same 27 too — window A is correct *only* because
   the context path is taken there. Making native "window-invariant" by always taking the fallback would
   therefore break window A as well.

The deeper defect is the **sentinel** (`calling_pipeline.cpp:10548`): the `*` candidate never appears in
`candidate_reads`, hence never gets
`likelihood_candidate_read_context_ordinals[candidate]` (`calling_pipeline.cpp:10566`), hence
`have_context_mapping_evidence` is false for the whole annotation call
(`calling_pipeline.cpp:4112-4153`), hence *both* the MQ gate and the BestAllele-based strand counters lose
their population.

---

## Q5 — Candidate patches

### Q5.0 The code under consideration (read in full)

* `calling_pipeline.cpp:4108-4160` — context mapping; `if (context_ordinal == UINT32_MAX) continue;` at 4129;
  `have_context_mapping_evidence = true;` at 4152.
* `calling_pipeline.cpp:4280-4312` — the strand gate
  `retained_for_hc_allele_annotations = !has_likelihood_rows || read_overlaps_hc_genotyping_interval(...)`
  (4280-4282), `annotation_read_survives` (4288-4292), the MQ gate `mapping_evidence` (4293-4297),
  `if (!retained_for_hc_allele_annotations) continue;` (4310).
* `calling_pipeline.cpp:10540-10620` — `candidate_read_context_ordinals` allocation/assignment and
  `candidate_read_realignments` construction.
* `calling_pipeline.cpp:12141-12189` — `qualified_annotation_depth` (the other consumer of the context ordinal).
* `calling_pipeline.cpp:16950-16980` — the `Result` fields that carry candidate contexts and realignments.
* `hc_call.cpp:4505-4532` — render-time recomputation
  (`annotations = calls[best]->annotations;` then `if ((include_spanning_deletion || group.max_alt_subset) && annotation_reads != nullptr) { ... annotations = *recalculated; }`).

Additional facts established while reading (they decide the evaluation):

* `allele_read_likelihoods` / `reference_read_likelihoods` are allocated for **every** candidate, initialised
  to `-inf` (`calling_pipeline.cpp:7211-7217`), and are written **only from the `candidate_reads` loop**
  (`calling_pipeline.cpp:11867-11953`). A candidate with no `candidate_reads` row therefore has
  all-`-inf` rows.
* `has_finite_likelihood_row` requires **both** `alt_likelihood` and `ref_likelihood` to be finite
  (`calling_pipeline.cpp:4224-4225`); when it is false the BestAllele classification is skipped entirely
  (4240-4269) and every strand counter stays 0 (`calling_pipeline.cpp:4315-4321`). **This is the mechanism of
  `SB=0,0,0,0` in window B, and it is independent of the interval predicates.**
* `candidate_loci` groups candidates by `(tid, position, reference)` and reserves one extra marginalization
  column for the symbolic `*` **per locus** (`calling_pipeline.cpp:11092-11116`); `marginalization_row_ids` is
  keyed by `(group_ordinal, locus_id, source_record)` (`calling_pipeline.cpp:11519-11532`), i.e. **one row per
  (context, locus, read), not per candidate**. So the `*` column of the locus is already computed for the
  concrete sibling's rows; it is only *never copied* into `allele_read_likelihoods[*]`.

### Q5.1 Evaluation of the four proposed candidates

#### (i) Make the fallback use `read_overlaps_hc_genotyping_interval` (line 6725)

**Incorrect.**
`read_overlaps_hc_genotyping_interval` differs from `read_overlaps_annotation_interval` only by the 1-based /
0-based shift (`variant_start = candidate.position + 1`, `calling_pipeline.cpp:6735-6742`), giving
`[10020678, 10020684]` instead of `[10020677, 10020683]`. It is still a **source-CIGAR** predicate and still
admits ~27 reads. Verdict: it reproduces GATK's 8 in **neither** window. Window A is unaffected (A never
enters the fallback branch, so the edit is a no-op there); window B changes from 27 to ~27/28, still ≈3.4×
GATK's 8. It is not a fix, only a one-base shift of a predicate that should not be deciding this population
at all. (It also has a second consumer — the post-PairHMM `retained_for_candidate` gate at
`calling_pipeline.cpp:11513-11516` — so the predicate itself should not be repurposed casually.)

#### (ii) Sentinel ⇒ treat the candidate as having no MQ evidence (`mapping_count = 0`)

**Incorrect.** It would make the MQ field window-invariant but *wrong in window B* (GATK: `28800,8`; this: no
`RAW_MQandDP` at all, since `annotateRawData` returns empty when `evidenceCount() == 0`). It does nothing for
`SB`, which is driven by the missing likelihood rows, not by the MQ gate. It also silently suppresses real
evidence for every genuine "pileup-only candidate with no context" case that the code's own comment at
4283-4287 was written to preserve. Net effect: window A unchanged, window B still ≠ GATK.

#### (iii) Fix upstream so the spanning-deletion candidate is paired in `candidate_reads` for window B too

**Correct in direction.** This is the only family of changes that can restore **both** fields, because:
it removes the sentinel at its source (`calling_pipeline.cpp:10566`), which restores
`context_mapping_evidence` (→ 8 reads for `RAW_MQandDP`), **and** it is the only way the
`allele_read_likelihoods[candidate]` rows stop being all-`-inf` (→ finite BestAllele → `SB=0,0,3,3`).
As stated ("pair it in `candidate_reads`"), though, it is under-specified and can be implemented at very
different blast radii: appending *all* `contig_candidates` for assembly-region groups (essentially replacing
`ordered` with `contig_candidates` at `calling_pipeline.cpp:9152`) would attach every non-assembled pileup /
halo candidate in the region, which is much wider than needed. See RECOMMENDED PATCH for the narrow version
(only position-sharing siblings of already-assembled candidates).

#### (iv) At `hc_call.cpp:4530`, only overwrite the stored annotation when the recomputed one is non-empty

**Incorrect / a no-op for this bug.** The recomputed annotation for this record in window B is
**not** empty: it is `RAW_MQandDP=97200,27`, a fully populated annotation object. The literal reading of the
proposal therefore changes nothing. Even a "no evidence" reading does not help, because the in-run pass is
*also* zero-evidence for this record: running the native binary with
`FASTGATK_DEBUG_ANNOTATION_POSITION=10020679` shows, for window B,

```
[FASTGATK_ANNOTATION_SUMMARY] candidate=5 position=10020679 ... ref_forward=0 ref_reverse=0 alt_forward=0 alt_reverse=0 ...
[FASTGATK_ANNOTATION_SUMMARY] candidate=1 position=10020679 ... ref_forward=0 ref_reverse=0 alt_forward=3 alt_reverse=1 ...
[FASTGATK_ANNOTATION_SUMMARY] candidate=5 position=10020679 ... ref_forward=0 ref_reverse=0 alt_forward=0 alt_reverse=0 ...
```
versus window A, where the same candidate 5 produces `alt_forward=3 alt_reverse=1` and then
`alt_forward=3 alt_reverse=3` (the second pass feeds `SB=0,0,3,3`). Keeping the in-run value would keep a
zero-strand annotation, and the in-run `mapping_count` is produced by the same defective gate.
Blast radius is also non-trivial: it would freeze stale annotations in every legitimate zero-evidence case.

### RECOMMENDED PATCH

**Target:** `calling_pipeline.cpp:9152-9153` — make every symbolic sibling allele of an assembled EventMap
locus participate in that region's `candidate_reads` population.

Why this specific site: `pairhmm_event_candidates` is what `candidate_reads` is built from
(`calling_pipeline.cpp:9153-9326`), and `candidate_reads` is what sets
`candidate_read_context_ordinals[candidate]` (`calling_pipeline.cpp:10548-10566`) and what fills
`allele_read_likelihoods[candidate]` (`calling_pipeline.cpp:11867-11953`). The `*` allele is emitted by the
writer as a sibling of the concrete ALT (`hc_call.cpp:4534-4538`, `include_spanning_deletion`) and shares the
concrete sibling's locus (`calling_pipeline.cpp:11092-11104`), but it never seeds a local haplotype path, so it
is absent from `ordered`.

Blast radius (analysis, not measurement — see uncertainty):
* No new PairHMM requests: `request_indices` de-duplicates `(read_id, haplotype)`
  (`calling_pipeline.cpp:9315-9322`).
* No new marginalization rows: `marginalization_row_ids` is keyed by `(group_ordinal, locus_id, source_record)`
  (`calling_pipeline.cpp:11519-11532`), and the sibling shares all three with an already-paired candidate. The
  `*` column is already in `best_by_row_allele` (`calling_pipeline.cpp:11111-11115, 11882-11907`).
* Only these get populated for the previously-unpaired candidate: `candidate_read_context_ordinals`,
  `candidate_read_realignments` (`calling_pipeline.cpp:10562-10613`), and the copies at 11949-11952.
* One additional gate the new candidate rows must pass: the post-PairHMM allele boundary at
  `calling_pipeline.cpp:11513-11516`

  ```cpp
  const bool retained_for_candidate = somatic_mode ||
      read_overlaps_hc_genotyping_interval(
          reads, candidate_read.source_record,
          candidates[candidate_read.candidate], informative_read_overlap_margin);
  if (!retained_for_candidate) continue;
  ```

  This is a **source-CIGAR** filter, and it is evaluated with the `*` candidate's own reference length. It
  cannot *remove* a marginalization row (rows are per `(group_ordinal, locus_id, source_record)`, and the
  concrete sibling already owns that row), so the `*` candidate's context/realignment still exist; but it can
  prevent the row copy at 11949-11952 for individual reads. In window A the identical structure already
  produces GATK-identical output, which is why this is expected to be benign — it is nevertheless the first
  thing to inspect if `RAW_MQandDP` becomes `28800,8` while `SB` stays `0,0,0,0`.

* Window A must be byte-identical: candidate 5 is already in `ordered` there (it has a real context ordinal),
  so the new branch adds nothing.

```diff
--- a/fastgatk-native/src/calling_pipeline.cpp
+++ b/fastgatk-native/src/calling_pipeline.cpp
@@ -9150,7 +9150,38 @@
         // assembled alleles; attaching them here would widen the public
         // candidate set without adding a state to any haplotype.
-        const auto& pairhmm_event_candidates = group.assembly_region ? ordered : contig_candidates;
-        for (const auto candidate : pairhmm_event_candidates) {
+        // GATK marginalizes ONE read x haplotype AlleleLikelihoods per merged
+        // VariantContext and calls retainEvidence() on it for all of that VC's
+        // alleles at once (HaplotypeCallerGenotypingEngine.java:192-197), so a
+        // symbolic sibling such as the spanning-deletion '*' published next to
+        // a concrete ALT owns a likelihood row and a retained read population
+        // of its own.  Such a sibling never seeds a local haplotype path and is
+        // therefore absent from `ordered`; leaving it out keeps
+        // candidate_read_context_ordinals at the UINT32_MAX sentinel
+        // (calling_pipeline.cpp:10548), which makes the annotation gate fall
+        // back to a source-CIGAR interval GATK never uses
+        // (calling_pipeline.cpp:4293) and leaves allele_read_likelihoods rows
+        // at -inf (calling_pipeline.cpp:7211, 11867).  Pair the siblings of
+        // every ordered candidate so each published allele owns the region's
+        // context; this adds no PairHMM request and no marginalization row,
+        // because both are keyed by locus and read, not by candidate.
+        std::vector<std::size_t> assembly_event_candidates;
+        const std::vector<std::size_t>* pairhmm_event_candidates = &contig_candidates;
+        if (group.assembly_region) {
+            assembly_event_candidates = ordered;
+            for (const auto candidate : contig_candidates) {
+                if (std::find(assembly_event_candidates.begin(),
+                              assembly_event_candidates.end(), candidate) !=
+                    assembly_event_candidates.end())
+                    continue;
+                const bool shares_event_locus =
+                    std::any_of(ordered.begin(), ordered.end(), [&](const auto sibling) {
+                        return candidates[sibling].position == candidates[candidate].position;
+                    });
+                if (shares_event_locus) assembly_event_candidates.push_back(candidate);
+            }
+            pairhmm_event_candidates = &assembly_event_candidates;
+        }
+        for (const auto candidate : *pairhmm_event_candidates) {
             // PairHMMLikelihoodCalculationEngine scores the complete
             // trimmed AssemblyRegion read set before its fragment-level
             // event marginalization.  Selecting only reads that touch each
```

Required verification after applying (I did **not** build or run anything; the tree is unmodified):

1. Rebuild the OpenMP HC target.
2. For `-L 20:10019901-10020710` (window A) diff against pinned GATK: **all 48 rows must stay byte-identical**,
   especially `10020680 RAW_MQandDP=28800,8 SB=0,0,3,3`, `INFO/DP`, `AD`, `PL`, `PGT/PID/PS`.
3. For `-L 20:10020381-10020710` (window B): the `10020680` row must become `28800,8` / `0,0,3,3`, and no other
   row may change.
4. Regression windows from the D2 sweep (`10020421`, `10020431`, which currently show the third state
   `SB=1,2,6,16`) should also converge to GATK.

### Safety net if the recommended patch is judged too invasive for this session

The only *narrow* change that is guaranteed not to touch window A and that removes the raw-CIGAR population
from the MQ gate is to reuse the sibling's context ordinal when the candidate's own is the sentinel, at
`calling_pipeline.cpp:4124-4153` — i.e. when `likelihood_candidate_read_context_ordinals[index]` is
`UINT32_MAX`, look up the realignments of another candidate in the same group that does have an ordinal, and
build `context_mapping_evidence` from those realigned intervals. This is semantically defensible (in GATK all
alleles of a merged VC share one evidence collection) and it fixes `RAW_MQandDP`.
**It will not fix `SB`**, because `SB` additionally requires finite `allele_read_likelihoods[candidate]` rows,
which only the upstream pairing change can supply. It should therefore be treated as a diagnostic/interim
change, not as parity.

---

## What remains uncertain (explicit)

1. **I did not compile or run the recommended patch.** The evaluation of its blast radius is a source-level
   argument (`request_indices` de-duplication at 9315-9322; `marginalization_row_ids` keyed by
   `(group_ordinal, locus_id, source_record)` at 11519-11532; `candidate_loci` grouping at 11092-11104), not a
   measurement. It must be validated by rebuilding and re-diffing both windows (and the wider regression set).
2. **Why exactly the `*` candidate is missing from `ordered` in window B is inferred, not directly observed.**
   The round-4 instrumentation proves the sentinel (`likelihood_candidate_read_context_ordinals[*] ==
   UINT32_MAX`), and `candidate_reads` is built only from `pairhmm_event_candidates` (9152-9153), and
   `ordered` is narrowed to the "combinable" set at `calling_pipeline.cpp:8107` (`ordered =
   std::move(combinable);`). I did not instrument `ordered` itself. **Speculation:** in window B the
   spanning-deletion path is not among the combinable candidates because the region's assembled haplotype
   population differs (the pinned GATK debug shows the assembly spans differ between the windows). The
   experiment that settles it: print `ordered` (indices + `candidates[i].position` + alternate) and compare it
   with `group.candidates` for the two windows.
3. **Whether the recommended patch also makes window B's `SB` equal `0,0,3,3` is argued, not measured.**
   It depends on the marginalization kernel producing a finite value in the `*` column
   (`calling_pipeline.cpp:11111-11115, 11892-11907`) for the sibling's rows in window B. If the `*` column is
   `-inf` there, `has_finite_likelihood_row` stays false and `SB` stays `0,0,0,0` even though
   `RAW_MQandDP` becomes `28800,8`; in that case the remaining gap is a haplotype-population parity issue
   (assembly), not an annotation-scope issue. A second place to check is the source-CIGAR row gate at
   `calling_pipeline.cpp:11513-11516`, which is evaluated with the `*` candidate's own reference length.
   The experiment that settles it:
   `FASTGATK_DEBUG_ALLELE_MARGINALIZED...` / the `[FASTGATK_ANNOTATION_SUMMARY]` trace
   (`calling_pipeline.cpp:11921-11933`, `4508-4520`) after the patch, checking finiteness of
   `allele_read_likelihoods[*]`.
4. **The `disqualified` accounting is inferred from the debug trace.** GATK prints
   `readAlleleLikelihoods.evidenceCount()` = 8 and `filteredSampleEvidence(0).size()` = 20
   (`HaplotypeCallerGenotypingEngine.java:214`). `retainEvidence` only filters the *already*-disqualified list
   (`AlleleLikelihoods.java:1191-1194`), so the 20 are reads that had already been moved there by
   `filterPoorlyModeledEvidence` and that overlap the interval. The region's PairHMM input dump shows 33 reads
   and the main-evidence dump shows 11, which is consistent with (11 main + 22 disqualified = 33) → (8 + 20 =
   28 after `retainEvidence`); the 5 reads neither retained nor disqualified are disqualified reads outside the
   interval, dropped entirely from the filtered list. This arithmetic is self-consistent but was deduced, not
   printed directly.
5. **Population semantics of `FORMAT/DP` in native are out of scope here.** GATK's `FORMAT/DP` is
   `informative BestAllele` count (`GenotypeBuilder`/`DP` from the genotyping matrix) and native derives it via
   `derive_multiallelic_depths` (`hc_call.cpp:4492-4509`); both currently agree with GATK (8) and must be
   re-checked after the patch, since a newly populated `allele_read_likelihoods[*]` could shift an informative
   count.

---

## Appendix A — raw GATK oracle (pinned 4.6.2.0 jar)

```
window A  -L 20:10019901-10020710
  assemblyRegion: 20:10019901-10020199
  assemblyRegion: 20:10020200-10020489
  assemblyRegion: 20:10020490-10020710
    Event at: [VC HC2 @ 20:10020679-10020680 Q. of type=MNP alleles=[AC*, TA] attr={} GT=[] filters= with 7 reads and 20 disqualified
    Event at: [VC HC0 @ 20:10020680-10020681 Q. of type=MIXED alleles=[CA*, *, AT] attr={} GT=[] filters= with 8 reads and 20 disqualified
  20  10020680  .  CA  AT,*,<NON_REF>  154.25  .  DP=8;ExcessHet=0.0000;MLEAC=1,1,0;MLEAF=0.500,0.500,0.00;RAW_MQandDP=28800,8
      GT:AD:DP:GQ:PGT:PID:PL:PS:SB  1|2:0,4,2,0:6:84:1|0:10020679_AC_TA:249,84,99,165,0,159,252,99,168,264:10020679:0,0,3,3

window B  -L 20:10020381-10020710
  assemblyRegion: 20:10020381-10020658
  assemblyRegion: 20:10020659-10020710
    Event at: [VC HC2 @ 20:10020679-10020680 Q. of type=MNP alleles=[AC*, TA] attr={} GT=[] filters= with 7 reads and 20 disqualified
    Event at: [VC HC0 @ 20:10020680-10020681 Q. of type=MIXED alleles=[CA*, *, AT] attr={} GT=[] filters= with 8 reads and 20 disqualified
  20  10020680  .  CA  AT,*,<NON_REF>  154.25  .  DP=8;ExcessHet=0.0000;MLEAC=1,1,0;MLEAF=0.500,0.500,0.00;RAW_MQandDP=28800,8
      GT:AD:DP:GQ:PGT:PID:PL:PS:SB  1|2:0,4,2,0:6:84:1|0:10020679_AC_TA:249,84,99,165,0,159,252,99,168,264:10020679:0,0,3,3
```

The 8 retained read names for the `20:10020680` event, identical in both windows:
`20GAVAAXX100126:6:3:6181:93202`, `20FUKAAXX100202:1:8:20035:107624`,
`20GAVAAXX100126:1:62:21171:66429`, `20FUKAAXX100202:4:43:2630:81606`,
`20GAVAAXX100126:8:44:10223:188523`, `20FUKAAXX100202:3:45:1441:112187`,
`20FUKAAXX100202:8:24:3734:131116`, `20FUKAAXX100202:7:3:12109:51312`.

## Appendix B — native evidence collected for this report (binary run, no rebuild)

```
FASTGATK_DEBUG_ANNOTATION_POSITION=10020679 fastgatk-native/build/fastgatk-hc-call ... -L 20:<start>-10020710

window A (10019901):
  [FASTGATK_ANNOTATION_SUMMARY] candidate=5 position=10020679 ... alt_forward=3 alt_reverse=1 ...
  [FASTGATK_ANNOTATION_SUMMARY] candidate=5 position=10020679 ... alt_forward=3 alt_reverse=3 ...
  20  10020680  .  CA  AT,*,<NON_REF>  ...  RAW_MQandDP=28800,8  ... SB 1|2:0,4,2,0:6:84:...:0,0,3,3

window B (10020381):
  [FASTGATK_ANNOTATION_SUMMARY] candidate=5 position=10020679 ... alt_forward=0 alt_reverse=0 ...
  [FASTGATK_ANNOTATION_SUMMARY] candidate=1 position=10020679 ... alt_forward=3 alt_reverse=1 ...
  [FASTGATK_ANNOTATION_SUMMARY] candidate=5 position=10020679 ... alt_forward=0 alt_reverse=0 ...
  20  10020680  .  CA  AT,*,<NON_REF>  ...  RAW_MQandDP=97200,27 ... SB 1|2:0,4,2,0:6:84:...:0,0,0,0
```

`fixtures/chr20/mnp.bam`, reads overlapping `20:10020678-10020683` (pysam): 33 total, **6 duplicate-flagged**
→ 27 non-duplicate, all MAPQ 60. `27 * 60^2 = 97200`, which is exactly the native fallback value: the fallback
population is the raw-BAM overlap set, and the duplicate filter is applied (27, not 33).
