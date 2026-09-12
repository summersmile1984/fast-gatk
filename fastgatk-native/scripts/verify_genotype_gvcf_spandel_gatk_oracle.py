#!/usr/bin/env python3
"""Strict pinned-GATK oracle for the orphan spanning-deletion gVCF fixture used
by ``fastgatk-native/scripts/verify_genotype_gvcf.py``.

The gap this gate pins
----------------------
``verify_genotype_gvcf.py`` builds, at lines 373-381, a single-sample gVCF
record whose ALT list is ``*,G,<NON_REF>`` and asserts at line 391 that native
GenotypeGVCFs emits exactly one record whose ALT contains the symbolic
spanning deletion ``*``.

Pinned GATK 4.6.2.0 emits **zero** records for that fixture by default and, with
``--include-non-variant-sites`` (short name ``-all-sites``, both spellings bound
to the same ``@Argument`` at ``GenotypeGVCFs.java:116-117,126-127``), exactly one
record that carries ``ALT='.'`` and ``GT='./.'``.  The old assertion therefore
pinned native-only output rather than the GATK contract.

GATK's rule (read from ``gatk-source/`` and then confirmed by measurement)
-------------------------------------------------------------------------
1. ``GenotypingEngine.calculateOutputAlleleSubset`` only publishes an ALT that
   is individually plausible at ``--standard-min-confidence-threshold-for-calling``
   (default 30), and it drops a symbolic ``*`` outright when no emitted deletion
   covers the locus (``GenotypingEngine.java:311-327``).  For this fixture both
   ``*`` and ``G`` are pruned, so the surviving ALT set is empty.
2. With an empty ALT set the site stays "monomorphic" and
   ``passesEmitThreshold()`` is false because the best guess is the reference
   (``GenotypingEngine.java:425-427``), so ``calculateGenotypes`` returns
   ``null`` and the default traversal writes nothing
   (``GenotypingEngine.java:167-169``).
3. ``--include-non-variant-sites`` sets ``OutputMode.EMIT_ALL_ACTIVE_SITES``
   (``GenotypeGVCFsEngine.java:381-383``), which makes ``emitAllActiveSites()``
   true (``GenotypingEngine.java:421-423``) and skips that ``null`` return.  The
   rebuilt VariantContext is REF-only: ``subsetToRefOnly()`` supplies the
   genotypes (``GenotypingEngine.java:190``) and
   ``cleanupGenotypeAnnotations(result, true, false)`` then replaces each of
   them with ``./.`` and drops every FORMAT field but GT
   (``GenotypeGVCFsEngine.java:191-194`` and ``:479-491``).  ``INFO/DP`` survives
   from the merged record and the Number=A ``MLEAC``/``MLEAF`` vectors are
   written as missing values, giving the literal row

       chr1 2 . A . 127.78 . DP=20;MLEAC=.;MLEAF=. GT ./.

   ``GenotypeGVCFs.apply()`` writes it because the record is not
   "spanning-deletion only" -- it has no ALT at all
   (``GenotypeGVCFs.java:324-330``, ``GATKVariantContextUtils.java:2089-2091``).

The orphan-``*``-with-surviving-ALT cases
-----------------------------------------
``surviving-concrete-alt`` was carried as REPORTED ONLY until the round
documented in ``.diag/round-prefer-pls.md`` promoted it.  It pins the other half
of ``GATKVariantContextUtils.makeGenotypeCall()``'s PREFER_PLS split: when the
likelihood row that ``AlleleSubsettingUtils.subsetAlleles()`` projected onto the
KEPT alleles is not informative (``isInformative`` is ``sum(log10GL) < -0.1``,
``GATKVariantContextUtils.java:54-58``), the row is ignored and the call is
projected from the source genotype with ``bestMatchToOriginalGT()``
(``:333-338``, definition ``:397-403``) -- a surviving source allele keeps its
identity, every pruned allele becomes the reference, the source copy order and
phase bit are preserved, and no GQ is assigned.  The companion cases pin the
position-independent projection (``...-star-not-first``), the copy order
(``...-reversed-source-gt``) and the phase bit (``...-phased-source-gt``), which
together distinguish "project the source call" from "derive the call from the
PLs" (whose argmax over this fixture is the ``(*,G)`` cell).

The deletion-ownership cases
----------------------------
Whether a symbolic ``*`` is "spurious" is an ownership question, not a
likelihood question.  ``GenotypingEngine.calculateOutputAlleleSubset()`` computes

    isSpuriousSpanningDeletion = isSpanningDeletion(allele) && !isVcCoveredByDeletion(vc)

(``GenotypingEngine.java:314``) and drops the allele when it is spurious
(``:316``).  ``isVcCoveredByDeletion()`` (``:365-371``) requires

    loc.getStart() < vc.getStart() && vc.getStart() <= loc.getEnd()

over the deletions collected by ``recordDeletions()`` (``:343-357``), and
``recordDeletions()`` is only reached at ``:179`` -- *after* the output-allele
subset was already decided at ``:155``.  Two consequences are pinned here:

1. a deletion allele carried by the very record being genotyped can never own
   that record's own ``*`` (it is not in the upstream set yet, and the strict
   ``<`` would exclude it anyway).  ``surviving-deletion-alt`` is exactly that
   shape (``AA`` / ``*,A,<NON_REF>``, deletion ``A`` owns the best genotype):
   GATK prunes ``*``, publishes ``ALT='A'`` and the source call projected to
   ``0/1``, with ``FORMAT=GT:AD:DP:PL`` and no ``GQ``.
2. an *upstream* deletion still owns a downstream ``*``, but only when it
   genuinely spans the locus: ``...-noncovering-upstream`` (deletion spans 2-3,
   ``*`` at 4) prunes, ``...-covering-upstream`` (deletion spans 2-4, ``*`` at 4)
   keeps ``*`` and publishes ``1/2``.

Native's ownership state is a pre-computed list of every input record carrying a
concrete deletion allele, tested with ``span.begin <= record.pos && record.pos <
span.end``; the inclusive start made a record's own deletion own its own ``*``.
The fix makes the start strictly exclusive, matching ``:369``.

The spanning-deletion-only record
---------------------------------
Ownership is not the last filter.  Two further rules apply once the surviving ALT
list is known, and neither depends on who owns the ``*``:

1. ``GenotypingEngine`` refuses to build the call at all when the output allele
   set is exactly ``[SPAN_DEL]`` and the traversal is not
   ``EMIT_ALL_ACTIVE_SITES``: ``if (!emitAllActiveSites() &&
   outputAlternativeAlleles.alleles.size() == 1 &&
   Allele.SPAN_DEL.equals(...)) return null;`` (``GenotypingEngine.java:173-175``;
   ``emitAllActiveSites()`` at ``:421-423``, configured by
   ``GenotypeGVCFsEngine.createMinimalArgs()`` at ``:381-383``).
2. ``GenotypeGVCFs.apply()`` writes the regenotyped record only when
   ``forceOutput || !GATKVariantContextUtils.isSpanningDeletionOnly(...)``
   (``GenotypeGVCFs.java:326-328``); ``isSpanningDeletionOnly()`` is
   ``getAlternateAlleles().size() == 1 && isSpanningDeletion(allele0)``
   (``GATKVariantContextUtils.java:2089-2091``).  ``forceOutput`` is
   ``includeNonVariants || inForceOutputIntervals`` (``:324-326``).

So a locus whose only surviving allele is ``*`` is never published by default,
while under ``--include-non-variant-sites`` it is published as the REF-only
no-call row: the allele subset is still computed by
``calculateOutputAlleleSubset()``, and ``passesThreshold()`` applies to ``*``
exactly as it does to a concrete ALT (``:312-316``).  ``covered-star-only-record``
therefore emits only the upstream deletion row, and
``covered-star-implausible-plus-concrete-alt`` shows that an implausible ``*`` is
dropped *allele-wise* while the record survives on its concrete ALT.

One divergence in this area is deliberately left REPORTED ONLY, because it is a
different root cause: ``recordDeletions()`` only ever records deletions that were
actually **emitted** (``:178-179``), while native's ownership state is a
pre-computed span list built from the **input** records.  When GATK drops an
upstream deletion record for failing the allele threshold, the downstream ``*``
becomes spurious and is pruned, but native still counts it as owned:
``unemitted-upstream-deletion-star-plus-concrete-alt``.  Repairing it needs the
ordered per-locus "emitted deletions" state, so it is out of scope here.

The dense-mode QUAL token
-------------------------
The REF-only row that dense mode materializes carries QUAL ``Infinity`` when the
locus is monomorphic *and* the model puts the whole posterior mass on "no
variant present".  ``GenotypingEngine`` selects the confidence formula by
monomorphy (``:158-163``):

    log10Confidence = !outputAlternativeAlleles.siteIsMonomorphic
                          || configuration.annotateAllSitesWithPLs
                      ? AFresult.log10ProbOnlyRefAlleleExists() + 0.0
                      : AFresult.log10ProbVariantPresent() + 0.0;

and it is that second arm which is assigned to the record
(``builder.log10PError(log10Confidence)`` at ``:183``, with the phred value only
used for the LowQual test at ``:184``).  ``log10ProbOnlyRefAlleleExists()`` is
``-0.0`` here (the posterior is 1.0), and
``MathUtils.log10OneMinusPow10`` returns ``Double.NEGATIVE_INFINITY`` for a zero
argument, so the QUAL is a genuine ``Double.POSITIVE_INFINITY``.  htsjdk renders
it through ordinary decimal formatting -- ``VCFEncoder`` writes
``formatQualValue(vc.getPhredScaledQual())`` whenever ``vc.hasLog10PError()``,
and ``formatQualValue`` is ``String.format(Locale.US, "%.2f", qual)`` with a
trailing ``.00`` stripped -- so Java's ``%f`` prints the token ``Infinity``
(``-Infinity``/``NaN`` would print their own tokens).  The token is therefore
value-driven, not a special case in the writer: the finite sibling
(``include-non-variant-sites`` over ``STAR_RECORD``) is the *same* branch with a
finite ``log10ProbVariantPresent()`` and reads ``127.78``.

The FILTER column
-----------------
GATK's genotyper does not inherit the source record's filters.  It rebuilds the
call from scratch:

    final VariantContextBuilder builder = new VariantContextBuilder(
            callSourceString(), vc.getContig(), vc.getStart(), vc.getEnd(),
            outputAlleles);            // GenotypingEngine.java:181
    builder.log10PError(log10Confidence);
    if ( ! passesCallThreshold(phredScaledConfidence) ) {
        builder.filter(GATKVCFConstants.LOW_QUAL_FILTER_NAME);   // :184-186
    }

because that constructor is given only the source string, the coordinates and
the output alleles, no filter state is copied from ``vc``: the new context is
unfiltered and ``filtersWereApplied`` is false, so htsjdk's ``VCFEncoder``
writes ``.`` (its ``addFilterString`` emits ``.`` exactly when
``!vc.filtersWereApplied()``).  The only filter a GenotypeGVCFs output record
can carry is the engine's own ``LowQual`` from ``passesCallThreshold``, a test
on the *recomputed* phred confidence against
``--standard-min-confidence-threshold-for-calling`` (``:430-431``).

The dense-mode non-PASS case therefore reads

    chr1 3 . A . Infinity . DP=7;MLEAC=.;MLEAF=. GT ./.

even though the source leaf is ``FILTER=LowQual`` -- the QUAL verifies that the
recomputed confidence passes the threshold, so no filter is applied.  The
measured divergence this gate was added for is native writing the inherited
filter instead, rendered as a *wrong name*: HTSlib keeps one shared ID
dictionary for ``##FILTER``/``##INFO``/``##FORMAT``, the source record was
decoded with a header in which HTSlib had just auto-registered the undeclared
``LowQual`` (``vcf_parse_filter`` appends a dummy ``##FILTER`` line while the
record is parsed), and that index is resolved against the *writing* header at
``vcf_format1``/``bcf_write`` time -- a header duplicated before any record was
parsed.  Measured proof: with ``##FORMAT=<ID=RGQ,...>`` present the token reads
``RGQ``, with that declaration removed the same row reads ``GQ``, and an input
FILTER with an entirely different name (``StrandBias``) still reads ``RGQ``.
Two layers are pinned separately below: the undeclared-header case (wrong name)
and the declared-header case (right name, still wrong because GATK writes ``.``).

The half of the contract native does not implement -- the engine's own
``LowQual`` -- is pinned as REPORTED ONLY by
``weak-locus-lowqual-filter-not-implemented``.

The output header's ``##FILTER`` lines
--------------------------------------
The FILTER *column* is not the whole contract: GATK's writer **always declares
the filter it may apply**.  ``GenotypeGVCFsEngine.setupVCFWriter()`` ends its
header construction with

    headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.LOW_QUAL_FILTER_NAME));
                                                    // GenotypeGVCFsEngine.java:416

as the last ``add()`` before ``new VCFHeader(headerLines, ...)``
(``:418-419``) and before ``vcfWriter.writeHeader(outputHeader)`` (``:420``).
Nothing guards that call: it is not conditional on the input header, on
``--include-non-variant-sites``, on ``--keep-combined``, on ``--dbsnp`` or on
the output mode, and ``setupVCFWriter()`` has exactly one caller
(``GenotypeGVCFs.java:305``), so the line is emitted for every run of the tool.
The text is the standard GenotypeGVCFs filter line
``##FILTER=<ID=LowQual,Description="Low quality">``
(``GATKVCFHeaderLines.java:89``:
``addFilterLine(new VCFFilterHeaderLine(LOW_QUAL_FILTER_NAME, "Low quality"))``
with ``LOW_QUAL_FILTER_NAME = "LowQual"`` at ``GATKVCFConstants.java:179``).

Position: htsjdk's ``VCFHeader`` writer emits its lines in sorted order (key
before value, ``##fileformat`` pinned first), so the whole ``##FILTER`` group
appears immediately after ``##fileformat``/``##ALT`` and the lines *inside* the
group are ordered by their full text -- ``LowQual`` therefore precedes an input
declaration such as ``q10`` (``'L'`` < ``'q'``) and follows one such as ``AAA``.
Measured on the fixture of this gate (``.diag/header-filter-probe.log``): with
no ``##FILTER`` line in the input, GATK's output header carries exactly one, at
index 2 of 31; with ``##FILTER=<ID=q10,...>`` declared it carries exactly
``['##FILTER=<ID=LowQual,Description="Low quality">',
'##FILTER=<ID=q10,Description="Quality below 10">']`` at indices 2-3.

This gate therefore compares the ``##FILTER`` header lines of both outputs
**byte for byte, in order**, for every gated case.  Every other header
difference is recorded as an observation and does not gate: native's
compatibility header deliberately preserves input order and its own INFO rank
list rather than reproducing htsjdk's fully sorted header, and it never writes
``##GATKCommandLine`` (see ``header_observations`` in each case's result).

``##FILTER=<ID=PASS,...>`` is propagated, not generated
-------------------------------------------------------
The reserved ``PASS`` filter is not a special case for GenotypeGVCFs.  The
writer's header set is seeded from the input header's own lines --

    final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());
                                                    // GenotypeGVCFsEngine.java:395

-- and the only line the tool removes is the GVCF block band
(``headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().startsWith(GVCF_BLOCK))``,
``:398-399``); every addition afterwards is ``##INFO``/``##FORMAT`` except the
LowQual filter at ``:416``.  So an input ``##FILTER=<ID=PASS,Description="All
filters passed">`` is written back **verbatim**, whatever its ``Description``
says, and its position inside the group follows htsjdk's sorted emission
(``LowQual`` < ``PASS`` < ``q10``).

This is *GATK propagation*, not an htsjdk rule and not htsjdk synthesis: the
byte string ``All filters passed`` does not occur in the pinned
``gatk-package-4.6.2.0-local.jar`` nor in ``htsjdk-4.2.0.jar``, htsjdk's
``VCFHeader`` never references ``VCFConstants.PASSES_FILTERS_v4``, and GATK's
own PASS line (``GATKVCFHeaderLines.java:90``, used only by VQSR and
``LabeledVariantAnnotationsWalker.java:317``) carries a different description
(``"Site contains at least one allele that passes filters"``).  Measured on this
gate's fixture matrix (``.diag/filter-pass-probe-before.log``): an input with no
``##FILTER`` line at all produces **no** PASS line in GATK's output, which
settles that nothing synthesizes it.

Pinned GATK 4.6.2.0, measured (``##FILTER`` group, in order, with the whole
header index of the first element):

    input: (none)                       -> [LowQual]                  @2
    input: PASS "All filters passed"    -> [LowQual, PASS]            @2
    input: PASS "Some other PASS text"  -> [LowQual, PASS(that text)] @2
    input: PASS + q10                   -> [LowQual, PASS, q10]       @2
    input: q10 + PASS (input order)     -> [LowQual, PASS, q10]       @2
    input: LowQual only                 -> [LowQual]                  @2
    input: AAA, LowQual, q10, PASS      -> [AAA, LowQual, PASS, q10]  @2

The group is therefore asserted as a whole **ordered list**, so a case fails
both when a line is missing and when it is ordered differently.

The output header's other content lines
---------------------------------------
The ``##FILTER`` group was only one face of the header divergence.  Beyond the
FILTER column the writer also

1. **re-serializes every header line it parsed.**  htsjdk parses a
   ``##INFO=<...>``/``##FORMAT=<...>``/``##ALT=<...>`` line into a
   ``VCFCompoundHeaderLine`` (``VCFHeaderLineTranslator.parseLine``, htsjdk
   4.2.0) and writes it back out of its own fields, so a ``Description`` the
   input left unquoted comes back quoted.  Measured: the gate fixture's
   ``##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>`` is written by
   GATK as ``##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">``.
   The attribute order inside ``<>`` is *not* re-ordered -- htsjdk rejects an
   input that does not already use its order (*"Tag Description in wrong order
   (was #2, expected #4)"*, ``VCF4Parser.parseLine``), so a valid input never
   needs re-ordering, only re-quoting.
2. **declares the annotations it always uses, unconditionally.**  The writer
   seeds its header from the input (``:395``) and then adds
   ``annotationEngine.getVCFAnnotationDescriptions(false)`` (``:401``),
   ``genotypingEngine.getAppropriateVCFInfoHeaders()`` (``:402``), the two MLE
   lines and RGQ (``:404-406``) and the standard INFO/DP line (``:407``, *"needed
   for gVCFs without DP tags"*).  GATK's ``VCFHeaderLine`` set de-duplicates by
   the **whole line** (ID + Number + Type + Description), not by ID: measured,
   GATK's header carries the input's ``##INFO=<ID=DP,...Description="Read
   depth">`` *and* ``##INFO=<ID=DP,...Description="Approximate read depth; some
   reads may have been filtered">``, and likewise two ``##FORMAT=<ID=AD,...>``
   lines.  So a declaration is added when its exact text is absent, and a
   same-ID-different-text input line is **kept alongside** it.

Measured on this gate's fixture, GATK declares these lines that native did not
(cited to the GATK source that produces each one):

    ##INFO=<ID=BaseQRankSum,...>       BaseQualityRankSumTest + GATKVCFHeaderLines.java:151
    ##INFO=<ID=MQRankSum,...>          MappingQualityRankSumTest + GATKVCFHeaderLines.java:168
    ##INFO=<ID=ReadPosRankSum,...>     ReadPosRankSumTest + GATKVCFHeaderLines.java:195
    ##INFO=<ID=DP,..."Approximate read depth; some reads may have been filtered">
                                       VCFStandardHeaderLines via GenotypeGVCFsEngine.java:407
    ##INFO=<ID=MLEAC,...>              GATKVCFHeaderLines.java:148 via GenotypeGVCFsEngine.java:404
    ##INFO=<ID=MLEAF,...>              GATKVCFHeaderLines.java:149 via GenotypeGVCFsEngine.java:405
    ##FORMAT=<ID=AD,..."Allelic depths for the ref and alt alleles in the order listed">
                                       DepthPerAlleleBySample via VariantAnnotation.java:22-33

The ``##INFO=<ID=MLEAC/MLEAF,...>`` pair was the *wording* half: native declared
the same keys with a short description of its own
(``Description=Maximum likelihood allele count``), which GATK never writes.

This gate therefore also asserts, for every gated case, that the **multiset of
header content lines is equal** -- every ``##`` line except the two
producer-identity lines (``##GATKCommandLine``, which records the command line of
the GATK process that ran, and ``##source``, which names the program).  Order is
deliberately *not* part of that assertion (native preserves the input's order and
htsjdk sorts the whole header, a pre-existing divergence recorded as an
observation); the ``##FILTER`` group keeps its additional ordered assertion
above.

The ``--annotate-with-num-discovered-alleles`` data-row contract
--------------------------------------------------------------
``--annotate-with-num-discovered-alleles`` (``@Argument`` at
``GenotypeCalculationArgumentCollection.java:73``, backing field
``ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED`` at ``:74``, default ``false``) makes
``GenotypingEngine.composeCallAttributes()`` add exactly one ``INFO`` key:

    if ( configuration.genotypeArgs.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED ) {
        attributes.put(GATKVCFConstants.NUMBER_OF_DISCOVERED_ALLELES_KEY,
                       vc.getAlternateAlleles().size());
    }
                                        // GenotypingEngine.java:464-465

``NUMBER_OF_DISCOVERED_ALLELES_KEY`` is the literal ``NDA``
(``GATKVCFConstants.java:73``); the declaration GATK writes for it is
``GATKVCFHeaderLines.java:205``, and it is added **only** under that flag
because it goes through ``GenotypingEngine.getAppropriateVCFInfoHeaders()``
behind the same boolean (``:90-94``), called from
``GenotypeGVCFsEngine.setupVCFWriter()`` (``:401-402``).

The count is read from the **``vc`` argument of ``calculateGenotypes()``**, not
from ``reducedVC`` (built at ``:137-144`` for ``--max-alternate-alleles``) and
not from the published allele subset (``outputAlternativeAlleles``, computed at
``:155``), so it is the ALT count of the *merged input* record -- before ``*``
ownership pruning, before the confidence test at ``:312-316`` and before the
max-ALT reduction.  Three consequences, all measured
(``.diag/nda-probe.log``, ``.diag/nda-probe2.log``,
``.diag/nda-gate-probe-before.log``):

1. ``<NON_REF>`` is **not** counted, because the merged ``vc`` handed to the
   engine has already had it removed: ``GenotypeGVCFsEngine`` calls
   ``merger.merge(variantsToProcess, loc, ref.getBase(), true, false)``
   (``GenotypeGVCFsEngine.java:136``), whose 4th parameter is
   ``removeNonRefSymbolicAllele``, and ``collectTargetAlleles()`` then adds
   ``Allele.NON_REF_ALLELE`` only when that flag is false
   (``ReferenceConfidenceVariantContextMerger.java:339-345``).  A one-ALT input
   therefore reports ``NDA=1``, not 2: the value counts *concrete* discovered
   alleles.
2. a symbolic ``*`` **is** counted, because ``collectTargetAlleles()`` re-adds
   ``Allele.SPAN_DEL`` for a spanning event
   (``ReferenceConfidenceVariantContextMerger.java:340-342``): ``A *,G,<NON_REF>``
   reports ``NDA=2`` even when GATK publishes only ``G``, and
   ``A *,<NON_REF>`` reports ``NDA=1`` even when the dense-mode row is REF-only.
   The rule is "input alleles", not "published alleles".
3. a locus whose merged ALT set is empty writes **no** ``NDA`` at all: that
   denormalized REF-only record comes from the ``regenotypeVC`` branch for a
   non-variant ``originalVC`` (``GenotypeGVCFsEngine.java:148``/``:154``) and never
   reaches ``composeCallAttributes()``.  Measured row ``chr1 2 . A . . . DP=20
   GT:AD ./.:20`` (``.diag/nda-probe2.log``, ``homref-nonref-only-dense``).

Native's ``annotate_num_discovered_alleles()`` (``genotype_gvcf_tool.cpp``,
called first in the compute stage, before ``apply_gatk_max_alternate_alleles()``
and ``apply_gatk_output_allele_subset()``, so the ordering and the input-vs-
published rule were already right) guarded the call with
``record.alleles.size() < 3`` and so silently required at least **two** ALT
alleles; for every ``NDA=1`` shape it wrote nothing.  The ``nda-annotation-*``
cases below pin the flag's whole contract: the row bytes including the ``NDA``
key, the input-allele counting rule (pruning must not change the value), the
REF-only dense row, two samples, and the control that neither the key nor its
declaration appears without the flag.

Two records at one locus: the input shape GATK does not support
--------------------------------------------------------------
``GenotypeGVCFs`` takes exactly ONE input track -- "1) a single single-sample
GVCF 2) a single multi-sample GVCF created by CombineGVCFs or 3) a GenomicsDB
workspace created by GenomicsDBImport" (``GenotypeGVCFs.java:65-66``); a second
``-V`` is refused outright ("Argument 'V/variant' cannot be specified more than
once", measured).  In that documented input there is exactly **one record per
(contig, position, REF)** per locus, because CombineGVCFs merges the per-sample
records of a locus before writing (``CombineGVCFs extends
MultiVariantWalkerGroupedOnStart``, ``CombineGVCFs.java:82``, one merged
``vcfWriter.add(mergedVC)`` per locus at ``:412-420``; measured: two
single-sample GVCFs with ``G`` and ``T`` at chr1:2 combine into the single
record ``chr1 2 . A T,G,<NON_REF>``).

What GATK does when that invariant is violated depends on the traversal mode,
and neither mode merges the records:

* **default** (no ``--include-non-variant-sites`` and no
  ``--force-output-intervals``): ``GenotypeGVCFs.onTraversalStart()`` calls
  ``changeTraversalModeToByVariant()`` (``GenotypeGVCFs.java:284-285``), so
  ``VariantLocusWalker.traverse()`` takes its by-variant branch and calls

      apply(variant, Collections.singletonList(variant), ...)
                                          // VariantLocusWalker.java:132-142

  i.e. **each record is its own locus**, the merger receives a one-element list
  (``GenotypeGVCFsEngine.java:128`` -> ``:136``) and one output row is written
  per input record (``GenotypeGVCFs.java:320-331``).  Two records at chr1:2
  therefore produce **two rows**, and each row carries the other sample as a
  bare ``./.`` because that record held no call for it.
* **group-by-locus** (``--include-non-variant-sites`` or
  ``--force-output-intervals``): ``apply()`` receives every variant overlapping
  the one-base locus (``VariantLocusWalker.java:154-172``) and then
  ``GenotypeGVCFsEngine.getVariantSubsetToProcess()`` **throws**:

      // since this tool only accepts a single input source, there should never be
      // more than one variant at a given starting locus
      throw new IllegalStateException(
              String.format(
                      "Variant input contains more than one variant starting at location: %s",
                      new SimpleInterval(matchingStart.get(0))));
                                          // GenotypeGVCFsEngine.java:349-368

  Measured: exit 3, no data row.  So the merge that
  ``ReferenceConfidenceVariantContextMerger.merge()`` performs over several
  records (a locus-level union: ``collectTargetAlleles()`` at ``:325-348``) is
  reachable only in this mode, and only for records that do **not** start at the
  locus (a spanning event, ``:150-151``) -- never for two records that start
  there.

Native instead coalesces every record sharing a ``record_key`` of
``rid:pos:REF`` (``genotype_gvcf_tool.cpp:1180-1187``; aggregate group loop
``:6538-6650`` and streaming group loop ``:5859-5880``/``:5969-5972``), keeps
the FIRST record of the group, drops every later record whose sample set
overlaps an already-accepted one (``:6639-6646``, streaming ``:5969-5972``) and
publishes one row per group.  That is the CombineGVCFs
model, not GenotypeGVCFs' by-variant model, so on an input with two records at
one locus the two tools cannot agree: GATK writes one row per record, native one
row per locus.  On the documented input -- one record per locus, ALT union
already computed -- the two models coincide and the rows are byte-identical,
which is what ``cross-sample-alt-union-single-record`` gated below pins.

The three ``same-position-two-records-*`` cases are therefore REPORTED ONLY:
they pin the measured GATK truth for the unsupported shape (including its
``NDA=1`` rows and the group-by-locus rejection) without asserting parity, so
the gate records the divergence and still passes.  See
``.diag/round-crosssample-merge.md`` for the full measurement.

Scope and comparison contract
-----------------------------
Pinned GATK and native run with identical arguments on the same plain
(unindexed .vcf.gz-refusing) VCF input and the same generated reference; the
data rows must be byte-identical.  Native additionally receives
``--gatk-compatible-annotations`` because the default native profile is an
explicitly non-GATK diagnostic output (it publishes RCQ/RCP and skips GATK's
AF-based output-allele pruning); every registered genotype-gvcf oracle uses the
same convention.

Exit status: 0 when every gated case passes, non-zero otherwise.
``--expect-divergence`` turns the run into a diagnostic that always exits 0.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import json
import os
import pathlib
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import oracle_guard  # noqa: E402

# Copied verbatim from verify_genotype_gvcf.py:328-337 so the gate pins the
# same fixture bytes as the registered contract test.
HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##INFO=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSTAR
"""

# The same header with two samples, for the multi-sample statement of the
# dense-mode rule.  GATK's group-by-locus traversal merges both samples into one
# record, and the materialized REF-only row carries one './.' per sample.
HEADER_TWO_SAMPLES = HEADER.replace("\tSTAR\n", "\tS1\tS2\n")

# The same header with an explicit FILTER declaration, so the source record's
# FILTER id is registered while the header is being READ rather than appended by
# HTSlib's vcf_parse_filter() when the first record is decoded.  The two
# spellings separate "the token names the wrong id" from "the source FILTER is
# carried at all": the undeclared header renders the inherited index against a
# different ID dictionary (see the FILTER section of this docstring), the
# declared one renders the right name, and GATK writes `.` for both.
HEADER_WITH_LOWQUAL_FILTER = HEADER.replace(
    "##FORMAT=<ID=GT", '##FILTER=<ID=LowQual,Description="Low quality">\n'
                       "##FORMAT=<ID=GT")

# The same header with an input filter that is NOT LowQual, so the gate also
# pins the *order* of the ##FILTER group: GATK's writer emits LowQual
# unconditionally (GenotypeGVCFsEngine.java:416) and htsjdk sorts the group by
# full line text, so GATK writes LowQual BEFORE q10 ('L' < 'q') even though q10
# is the only one the input declared.
HEADER_WITH_Q10_FILTER = HEADER.replace(
    "##FORMAT=<ID=GT", '##FILTER=<ID=q10,Description="Quality below 10">\n'
                       "##FORMAT=<ID=GT")

# The input declares the reserved ``PASS`` filter.  This is a header-propagation
# case, not a filter-generation case: GATK's output header is built as
# ``new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder())``
# (GenotypeGVCFsEngine.java:395) and nothing removes a ``FILTER`` line but
# ``GVCFBlock`` (:398-399), so whatever the input declared -- PASS included -- is
# written back out verbatim.  ``All filters passed`` is not a GATK or htsjdk
# constant (neither the pinned gatk-package jar nor htsjdk 4.2.0 contains that
# byte string), so the line and its text can only come from the input.
HEADER_WITH_PASS_FILTER = HEADER.replace(
    "##FORMAT=<ID=GT", '##FILTER=<ID=PASS,Description="All filters passed">\n'
                       "##FORMAT=<ID=GT")

# The same declaration with a DIFFERENT Description.  GATK keeps the input's text
# verbatim -- it neither normalizes nor replaces a ``PASS`` line, because neither
# it nor htsjdk ever synthesizes one (only GenotypeGVCFsEngine.java:416's LowQual
# line is added, and that is a different ID).
HEADER_WITH_NONSTANDARD_PASS_FILTER = HEADER.replace(
    "##FORMAT=<ID=GT", '##FILTER=<ID=PASS,Description="Some other PASS text">\n'
                       "##FORMAT=<ID=GT")

# PASS together with a non-reserved filter: the whole group must be reproduced in
# htsjdk's sorted order, i.e. LowQual (added by the tool) first, then PASS, then
# q10 -- measured, see the docstring section on the output header's FILTER lines.
HEADER_WITH_PASS_AND_Q10_FILTER = HEADER.replace(
    "##FORMAT=<ID=GT", '##FILTER=<ID=PASS,Description="All filters passed">\n'
                       '##FILTER=<ID=q10,Description="Quality below 10">\n'
                       "##FORMAT=<ID=GT")

# GenotypeGVCFsEngine.java:416 + GATKVCFHeaderLines.java:89 +
# GATKVCFConstants.java:179, byte for byte.
GATK_LOWQUAL_FILTER_LINE = '##FILTER=<ID=LowQual,Description="Low quality">'

# GATKVCFHeaderLines.java:205 + GATKVCFConstants.java:73, byte for byte.  Unlike
# the lines above this one is conditional: it is added through
# GenotypingEngine.getAppropriateVCFInfoHeaders() (GenotypingEngine.java:90-94)
# only when --annotate-with-num-discovered-alleles is set, so it is asserted
# per-case against that flag rather than against every case.
GATK_NDA_HEADER_LINE = (
    '##INFO=<ID=NDA,Number=1,Type=Integer,Description="Number of alternate '
    'alleles discovered (but not necessarily genotyped) at this site">')

# The same fixture header with every Description already in htsjdk's canonical
# quoted form, i.e. the shape a GATK-produced gVCF has.  It separates the two
# halves of the header divergence: here nothing needs re-quoting, so only the
# missing *declarations* can make the case fail.
HEADER_QUOTED_DESCRIPTIONS = HEADER
for _description in ("Represents any possible alternate allele", "Read depth",
                     "Allele depths", "Genotype", "Likelihoods",
                     "Genotype quality"):
    HEADER_QUOTED_DESCRIPTIONS = HEADER_QUOTED_DESCRIPTIONS.replace(
        f"Description={_description}>", f'Description="{_description}">')

# The fixture header with an MLEAC/MLEAF pair already declared using GATK's own
# (long) wording, which is what a HaplotypeCaller gVCF carries.  GATK's header
# set de-duplicates by the whole line, so this input's lines are kept and the
# tool's identical pair is NOT added a second time -- the case therefore also
# pins that native must not over-declare MLEAC/MLEAF.
HEADER_WITH_GATK_MLE_LINES = HEADER.replace(
    "##FORMAT=<ID=GT",
    '##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood '
    'expectation (MLE) for the allele counts (not necessarily the same as the '
    'AC), for each ALT allele, in the same order as listed">\n'
    '##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood '
    'expectation (MLE) for the allele frequency (not necessarily the same as '
    'the AF), for each ALT allele, in the same order as listed">\n'
    "##FORMAT=<ID=GT")

# Every header line GATK 4.6.2.0 GenotypeGVCFs declares for itself, i.e. the
# lines of its output header that do not come from the input, measured byte for
# byte on this gate's fixture and cited to their GATK source.
# GenotypeGVCFsEngine.setupVCFWriter(): `annotationEngine.
# getVCFAnnotationDescriptions(false)` (:401), `genotypingEngine.
# getAppropriateVCFInfoHeaders()` (:402), MLEAC/MLEAF/RGQ (:404-406),
# `VCFStandardHeaderLines.getInfoLine(DEPTH_KEY)` (:407) and the LowQual filter
# (:416).  The annotation lines resolve through GATKVCFHeaderLines (or, for
# AD/DP, htsjdk's VCFStandardHeaderLines -- GATKVCFHeaderLines.java:18-43).
GATK_DECLARED_HEADER_LINES = [
    GATK_LOWQUAL_FILTER_LINE,
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
    '##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">',
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">',
    '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">',
    '##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">',
    '##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">',
    '##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher\'s exact test to detect strand bias">',
    '##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">',
    '##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">',
    '##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">',
    '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">',
    '##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">',
    '##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">',
    '##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">',
    '##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">',
]

# verify_genotype_gvcf.py:378-381, byte for byte.
STAR_RECORD = (
    "chr1\t2\t.\tA\t*,G,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:12,8,0,0:0,0,100,100,100,100,100,100,100,100\n"
)

# The concrete ALT G owns the best genotype (PL index 4 = the (*,G) cell), so
# only the orphan '*' is pruned.  GATK keeps ALT='G' and publishes 0/1.
G_PLAUSIBLE_RECORD = (
    "chr1\t2\t.\tA\t*,G,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/2:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100\n"
)

# Same locus with the orphan '*' in the middle of the ALT list, so the pruned
# allele is neither first nor last.  bestMatchToOriginalGT() projects the source
# call by allele identity, not by index, so the published row is unchanged.
G_STAR_SECOND_RECORD = (
    "chr1\t2\t.\tA\tG,*,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0,0:100,100,100,100,0,100,100,100,100,100\n"
)

# A concrete DELETION owns the best genotype of the very same record.  GATK
# prunes the orphan '*' here too (see the deletion-ownership section of this
# docstring), so ALT='A' and GT 0/1 are published.
G_DELETION_ALT_RECORD = (
    "chr1\t2\t.\tAA\t*,A,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/2:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100\n"
)

# Same shape, but the deletion allele sits at a locus that an UPSTREAM deletion
# does NOT cover (the record at 2 spans 2-3 while this one starts at 4), so the
# record's own deletion still cannot own its own '*'.
DEL_ALT_NONCOVERING_UPSTREAM = (
    "chr1\t2\t.\tAA\tA,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\n"
    "chr1\t4\t.\tAA\t*,A,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/2:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100\n"
)

# Positive control: the same '*' plus own deletion, but now an upstream deletion
# emitted at 2 spans 2-4 (REF AAA -> ALT A), which DOES cover locus 4, so GATK
# keeps '*,A' and genotypes 1/2.  This separates "the record's own deletion does
# not own it" from "no deletion may own it".
DEL_ALT_COVERING_UPSTREAM = (
    "chr1\t2\t.\tAAA\tA,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\n"
    "chr1\t4\t.\tAA\t*,A,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/2:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100\n"
)

# '*', an own deletion and nothing else plausible: pruning '*' leaves a single
# ALT that fails the standard-confidence threshold, the site turns monomorphic
# and the default traversal writes no record at all.
SELF_DELETION_ONLY_RECORD = (
    "chr1\t2\t.\tAA\t*,A,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:10,0,10,0:0,100,100,100,100,100,100,100,100,100\n"
)

# A '*' that a strictly upstream deletion emits as covering, and which therefore
# SURVIVES the ownership test -- but the record has no other ALT, so GATK never
# publishes it (GenotypingEngine.java:173-175 in the default output mode, and
# GenotypeGVCFs.java:327-328 for the dense mode that reaches apply()).
STAR_ONLY_COVERED_RECORD = (
    "chr1\t2\t.\tAA\tA,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\n"
    "chr1\t3\t.\tA\t*,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,100,100,0,100,100\n"
)

# The same two records with two samples: the dense-mode materialization is a
# property of the locus, not of the sample count, and GATK emits one './.' per
# sample on the REF-only row.
STAR_ONLY_COVERED_TWO_SAMPLES_RECORD = (
    "chr1\t2\t.\tAA\tA,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\t"
    "0/1:20:0,20,0:100,0,100,100,100,100\n"
    "chr1\t3\t.\tA\t*,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,100,100,0,100,100\t"
    "0/1:20:0,20,0:100,100,100,0,100,100\n"
)

# The same shipping rule must not depend on the dropped record's own FILTER,
# INFO/DP or AD: the star-only record is now non-PASS and shallower, and GATK
# still writes only the upstream deletion row.
STAR_ONLY_COVERED_NONPASS_RECORD = (
    "chr1\t2\t.\tAA\tA,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\n"
    "chr1\t3\t.\tA\t*,<NON_REF>\t.\tLowQual\tDP=7\t"
    "GT:DP:AD:PL\t0/1:7:0,7,0:100,100,100,0,100,100\n"
)

# The same non-PASS star-only locus with two samples, for the dense-mode FILTER
# statement of the multi-sample path.
STAR_ONLY_COVERED_NONPASS_TWO_SAMPLES_RECORD = (
    "chr1\t2\t.\tAA\tA,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\t"
    "0/1:20:0,20,0:100,0,100,100,100,100\n"
    "chr1\t3\t.\tA\t*,<NON_REF>\t.\tLowQual\tDP=7\t"
    "GT:DP:AD:PL\t0/1:7:0,7,0:100,100,100,0,100,100\t"
    "0/1:7:0,7,0:100,100,100,0,100,100\n"
)

# The ordinary (non-dense) variant path with a non-PASS source FILTER: the same
# record as G_PLAUSIBLE_RECORD, only the FILTER column changed.  GATK regenotypes
# it and publishes `.`; a source FILTER is never carried into the output.
G_PLAUSIBLE_NONPASS_RECORD = (
    "chr1\t2\t.\tA\t*,G,<NON_REF>\t.\tLowQual\tDP=20\t"
    "GT:DP:AD:PL\t0/2:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100\n"
)

# A locus whose best guess is the reference and whose phred confidence stays
# below --standard-min-confidence-threshold-for-calling (30): GATK applies its
# own FILTER=LowQual.  Used by the REPORTED ONLY case that pins the half of the
# FILTER contract native does not implement.
WEAK_LOCUS_RECORD = (
    "chr1\t2\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:19,1:0,0,40\n"
)

# A covered '*' that is NOT the only surviving ALT: the record also carries the
# concrete ALT G, whose (A,G) cell is the best of the source row.  Both alleles
# pass the count threshold and the '*' has an upstream owner, so GATK KEEPS the
# record and genotypes it 1/2.  Without this case a fix that dropped every
# record carrying a '*' would look correct.
COVERED_STAR_PLUS_CONCRETE_RECORD = (
    "chr1\t2\t.\tAA\tA,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\n"
    "chr1\t3\t.\tA\t*,G,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/2:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100\n"
)

# The same locus with the source row's best cell moved to (A,G): the '*' is now
# implausible as well as covered, so GATK prunes it from the output allele set
# (AF threshold at GenotypingEngine.java:312-316) but still publishes the
# surviving concrete ALT.  This separates "pruned because implausible" from
# "pruned because unowned" and proves the record-level rule is not "drop the
# whole record when it carries a '*'".
COVERED_STAR_IMPLAUSIBLE_PLUS_CONCRETE_RECORD = (
    "chr1\t2\t.\tAA\tA,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\n"
    "chr1\t3\t.\tA\t*,G,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/2:20:0,0,20,0:100,100,100,0,100,100,100,100,100,100\n"
)

# The previous round's fixture D, sharpened with a concrete ALT so that the
# residual structural difference is visible in the DEFAULT output mode.  The
# upstream deletion record at 2 is implausible, so GATK drops it and never calls
# recordDeletions() for it; the '*' at 3 is plausible, and native still counts
# the INPUT record at 2 as owning the locus.
IMPLAUSIBLE_UPSTREAM_STAR_PLUS_CONCRETE_RECORD = (
    "chr1\t2\t.\tAA\tA,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/0:20:20,0,0:0,100,100,100,100,100\n"
    "chr1\t3\t.\tA\t*,G,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/2:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100\n"
)

# The same locus with the source call written in the opposite copy order.  The
# measured GATK row keeps that order, so the fallback must project the source
# genotype rather than any PL-argmax or sorted genotype.
G_REVERSED_GT_RECORD = (
    "chr1\t2\t.\tA\t*,G,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t2/0:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100\n"
)

# ... and the phased spelling of the same call, which GATK also preserves.
G_PHASED_GT_RECORD = (
    "chr1\t2\t.\tA\t*,G,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t2|0:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100\n"
)

GATK_DEFAULT_ROW = "chr1\t2\t.\tA\t.\t127.78\t.\tDP=20;MLEAC=.;MLEAF=.\tGT\t./."

# Measured with pinned GATK 4.6.2.0 (see .diag/round-prefer-pls.md).
GATK_G_ROW = ("chr1\t2\t.\tA\tG\t82.26\t.\t"
              "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11\t"
              "GT:AD:DP:PL\t0/1:0,20:20:0,0,0")
GATK_G_REVERSED_ROW = ("chr1\t2\t.\tA\tG\t82.26\t.\t"
                       "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11\t"
                       "GT:AD:DP:PL\t1/0:0,20:20:0,0,0")
GATK_G_PHASED_ROW = ("chr1\t2\t.\tA\tG\t82.26\t.\t"
                     "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11\t"
                     "GT:AD:DP:PL\t1|0:0,20:20:0,0,0")

# Measured with pinned GATK 4.6.2.0 (see .diag/round-star-ownership.md).  When
# the record's own deletion allele is pruned from owning the locus, the orphan
# '*' goes away and only the concrete deletion allele is published.
GATK_DEL_ROW = ("chr1\t2\t.\tAA\tA\t82.19\t.\t"
                "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11\t"
                "GT:AD:DP:PL\t0/1:0,20:20:0,0,0")
GATK_DEL_UPSTREAM_LOCUS_ROW = ("chr1\t2\t.\tAA\tA\t92.60\t.\t"
                               "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63\t"
                               "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")
# The dense-mode materialization of a locus whose only surviving allele was the
# symbolic spanning deletion (measured, see .diag/round-star-only-record.md).
#
# The QUAL token is the literal `Infinity`.  It is not a special case in the VCF
# writer: GenotypingEngine:158-163 takes
#
#     log10Confidence = AFresult.log10ProbVariantPresent() + 0.0
#
# for a monomorphic site, and MathUtils.log10OneMinusPow10(0.0) is
# Double.NEGATIVE_INFINITY, so log10PError really is -Infinity and the phred
# value really is Double.POSITIVE_INFINITY.  htsjdk's encoder then renders it
# through plain decimal formatting -- VCFEncoder.formatQualValue() is
# `String.format(Locale.US, "%.2f", qual)` with a trailing ".00" removed
# (htsjdk VCFEncoder, verified by javap on the pinned gatk-package jar: the
# QUAL column is written as formatQualValue(vc.getPhredScaledQual()) when
# vc.hasLog10PError()) -- and Java's `%f` on an infinite double yields
# "Infinity".  Hence the unusual token.
GATK_STAR_ONLY_DENSE_ROW = ("chr1\t3\t.\tA\t.\tInfinity\t.\t"
                            "DP=20;MLEAC=.;MLEAF=.\tGT\t./.")
# The two-sample statement of the same locus: the first row is genotyped from
# both samples and the materialized REF-only row repeats the token.
GATK_DEL_UPSTREAM_TWO_SAMPLES_ROW = (
    "chr1\t2\t.\tAA\tA\t190.46\t.\t"
    "AC=2;AF=0.500;AN=4;DP=20;ExcessHet=1.7609;MLEAC=2;MLEAF=0.500;QD=4.76\t"
    "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100\t0/1:0,20:20:99:100,0,100")
GATK_STAR_ONLY_DENSE_TWO_SAMPLES_ROW = ("chr1\t3\t.\tA\t.\tInfinity\t.\t"
                                       "DP=20;MLEAC=.;MLEAF=.\tGT\t./.\t./.")
# The FILTER column of the same materialization when the source leaf is non-PASS
# (FILTER=LowQual, DP=7).  GATK writes `.`: GenotypingEngine builds the output
# record with `new VariantContextBuilder(callSourceString(), vc.getContig(),
# vc.getStart(), vc.getEnd(), outputAlleles)` (GenotypingEngine.java:181), which
# copies NOTHING from the source record, so `filtersWereApplied` stays false and
# htsjdk's VCFEncoder renders the unfiltered token `.`
# (VCFEncoder.addFilterString: `!vc.filtersWereApplied()` -> ".").  The source
# FILTER is only ever replaced by the engine's own decision at :184-186, which
# is a threshold test on the *recomputed* phred confidence -- Infinity here, so
# no filter is applied.
GATK_STAR_ONLY_DENSE_NONPASS_ROW = ("chr1\t3\t.\tA\t.\tInfinity\t.\t"
                                    "DP=7;MLEAC=.;MLEAF=.\tGT\t./.")
GATK_STAR_ONLY_DENSE_NONPASS_TWO_SAMPLES_ROW = (
    "chr1\t3\t.\tA\t.\tInfinity\t.\tDP=7;MLEAC=.;MLEAF=.\tGT\t./.\t./.")
# The ordinary variant path with a non-PASS source FILTER is byte-identical to
# GATK_G_ROW: the source FILTER never reaches the output, so only the QUAL
# decides, and 82.26 passes the call threshold.
GATK_G_NONPASS_ROW = GATK_G_ROW
# The positive half of GATK's FILTER rule, which native does not implement:
# phredScaledConfidence 23.14 < 30, so :184-186 applies LowQual.  Pinned as
# REPORTED ONLY (see the case's `why`).
GATK_WEAK_LOCUS_LOWQUAL_ROW = ("chr1\t2\t.\tA\t.\t23.14\tLowQual\t"
                               "DP=20;MLEAC=.;MLEAF=.\tGT\t./.")
# The finite sibling of the same branch, measured on the STAR_RECORD fixture
# (`GATK_DEFAULT_ROW`, 127.78): log10ProbVariantPresent() is a small NEGATIVE
# number there rather than -Infinity, so the row is finite.  Together the two
# rows pin that the token follows the value and is not a blanket rule for
# monomorphic sites.
# The same publication one locus further down (the record at 4 whose own
# deletion does not cover itself).
GATK_DEL_DOWNSTREAM_ROW = ("chr1\t4\t.\tAA\tA\t82.19\t.\t"
                           "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11\t"
                           "GT:AD:DP:PL\t0/1:0,20:20:0,0,0")
GATK_DEL_LONG_DELETION_ROW = ("chr1\t2\t.\tAAA\tA\t92.60\t.\t"
                              "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63\t"
                              "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")
# A '*' that a strictly upstream covering deletion owns is preserved, and the
# site is genotyped as a two-ALT call (unchanged from before this round).
GATK_STAR_AND_DEL_ROW = ("chr1\t4\t.\tAA\t*,A\t82.19\t.\t"
                         "AC=1,1;AF=0.500,0.500;AN=2;DP=20;ExcessHet=0.0000;"
                         "MLEAC=1,1;MLEAF=0.500,0.500;QD=4.11\t"
                         "GT:AD:DP:GQ:PL\t1/2:0,0,20:20:99:100,100,100,100,0,100")

# A '*' that a strictly upstream covering deletion owns, published together with
# a concrete ALT at the same locus.  Measured GATK 4.6.2.0 keeps both.
GATK_STAR_PLUS_CONCRETE_ROW = (
    "chr1\t3\t.\tA\t*,G\t82.26\t.\t"
    "AC=1,1;AF=0.500,0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1,1;MLEAF=0.500,0.500;QD=4.11\t"
    "GT:AD:DP:GQ:PL\t1/2:0,0,20:20:99:100,100,100,100,0,100")
# The same locus with the '*' implausible: only the concrete ALT is published.
GATK_CONCRETE_ONLY_ROW = (
    "chr1\t3\t.\tA\tG\t92.63\t.\t"
    "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63\t"
    "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")

# The same concrete-ALT publication one locus further down (the locus-3 form of
# GATK_G_ROW, used by the residual input-vs-emitted case).
GATK_G_DOWNSTREAM_ROW = ("chr1\t3\t.\tA\tG\t82.26\t.\t"
                         "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.11\t"
                         "GT:AD:DP:PL\t0/1:0,20:20:0,0,0")

# ---------------------------------------------------------------------------
# --annotate-with-num-discovered-alleles (NDA): fixtures and measured rows
# ---------------------------------------------------------------------------
# The flag is spelled identically on both sides, so it travels in the case's
# shared `args` and is passed to pinned GATK and native alike.  Native binds it
# at genotype_gvcf_tool.cpp:643-648, GATK at
# GenotypeCalculationArgumentCollection.java:73-74.

NDA_FLAG = "--annotate-with-num-discovered-alleles"

# One concrete ALT besides the symbolic <NON_REF>, called 0/1: the minimal
# shape on which GATK's value is 1 rather than 0 or 2.  It is the shape
# native's `record.alleles.size() < 3` guard silently skipped.
NDA_SINGLE_ALT_RECORD = (
    "chr1\t2\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\n"
)

# The same locus with two samples, so the count is pinned as a property of the
# locus and not of the sample count (GATK_NDA_TWO_SAMPLES_ROW has NDA=1, not 2).
NDA_TWO_SAMPLES_RECORD = (
    "chr1\t2\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\t"
    "0/1:20:0,20,0:100,0,100,100,100,100\n"
)

# Measured with pinned GATK 4.6.2.0 (.diag/nda-gate-probe-before.log).  The
# `NDA` key sits between `MLEAF` and `QD`: GenotypeGVCFsEngine.addGenotyping-
# Annotations() copies the attributes into a LinkedHashMap in the order MLEAC,
# MLEAF, NDA, AS_QUAL (GenotypeGVCFsEngine.java:233-250, the NDA carry at
# :238-239) before htsjdk's
# VCFEncoder writes them, which is the same order native's compatibility
# annotation writer already produces.
GATK_NDA_SINGLE_ALT_ROW = (
    "chr1\t2\t.\tA\tG\t92.64\t.\t"
    "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;NDA=1;QD=4.63\t"
    "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")
GATK_NDA_TWO_SAMPLES_ROW = (
    "chr1\t2\t.\tA\tG\t190.50\t.\t"
    "AC=2;AF=0.500;AN=4;DP=20;ExcessHet=1.7609;MLEAC=2;MLEAF=0.500;NDA=1;QD=4.76\t"
    "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100\t0/1:0,20:20:99:100,0,100")
# The control: the very same fixture and invocation WITHOUT the flag.  Neither
# the `NDA` key nor its `##INFO` declaration may appear on either side
# (GenotypingEngine.java:464-465 and :90-94 are both behind the same boolean).
GATK_NDA_CONTROL_ROW = GATK_NDA_SINGLE_ALT_ROW.replace("NDA=1;", "")
# G_PLAUSIBLE_RECORD ('*,G,<NON_REF>') publishes only `G` but has TWO input ALT
# alleles: NDA=2 proves the value follows the input set, not the published one.
GATK_NDA_REDUCED_ALT_SET_ROW = (
    "chr1\t2\t.\tA\tG\t82.26\t.\t"
    "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;NDA=2;QD=4.11\t"
    "GT:AD:DP:PL\t0/1:0,20:20:0,0,0")
# COVERED_STAR_PLUS_CONCRETE_RECORD: two rows, NDA=1 for the upstream deletion
# and NDA=2 for the two-ALT row ('*' counts as a discovered allele).
GATK_NDA_COVERED_STAR_PLUS_CONCRETE_ROWS = [
    ("chr1\t2\t.\tAA\tA\t92.60\t.\t"
     "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;NDA=1;QD=4.63\t"
     "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100"),
    ("chr1\t3\t.\tA\t*,G\t82.26\t.\t"
     "AC=1,1;AF=0.500,0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1,1;MLEAF=0.500,0.500;"
     "NDA=2;QD=4.11\t"
     "GT:AD:DP:GQ:PL\t1/2:0,0,20:20:99:100,100,100,100,0,100"),
]
# STAR_ONLY_COVERED_RECORD in dense mode: the materialized REF-only row keeps
# the count of the input record whose only ALT was the covered '*'.
GATK_NDA_COVERED_STAR_ONLY_DENSE_ROWS = [
    ("chr1\t2\t.\tAA\tA\t92.60\t.\t"
     "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;NDA=1;QD=4.63\t"
     "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100"),
    "chr1\t3\t.\tA\t.\tInfinity\t.\tDP=20;MLEAC=.;MLEAF=.;NDA=1\tGT\t./.",
]

# ---------------------------------------------------------------------------
# Cross-sample ALT sets at one locus: one record vs. two records
# ---------------------------------------------------------------------------
# The SUPPORTED shape: ONE record per locus whose ALT set already carries every
# sample's allele -- what CombineGVCFs writes (measured: two single-sample GVCFs
# through CombineGVCFs produce exactly one record `A  T,G,<NON_REF>`).  Both
# tools agree on it byte for byte; this is the gated half.
CROSS_SAMPLE_ALT_UNION_RECORD = (
    "chr1\t2\t.\tA\tG,T,<NON_REF>\t.\tPASS\tDP=40\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\t"
    "2/2:20:0,0,20:200,200,200,200,200,0\n"
)
GATK_CROSS_SAMPLE_ALT_UNION_ROW = (
    "chr1\t2\t.\tA\tG,T\t277.88\t.\t"
    "AC=1,2;AF=0.250,0.500;AN=4;DP=40;ExcessHet=0.0000;"
    "MLEAC=1,2;MLEAF=0.250,0.500;QD=6.95\t"
    "GT:AD:DP:GQ:PL\t0/1:0,20,0:20:99:100,0,100,100,100,100\t"
    "2/2:0,0,20:20:99:200,200,200,200,200,0")

# The UNSUPPORTED shape: TWO records starting at the same locus, each carrying
# one sample.  GATK's default traversal is by-variant (GenotypeGVCFs.java:284-285
# -> VariantLocusWalker.java:132-142), so it never merges them: it genotyped each
# record on its own and wrote one row per record, with the other sample rendered
# as a bare `./.`.  Native coalesces the group by (rid, pos, REF)
# (genotype_gvcf_tool.cpp:1180-1187) and writes one merged row, so these cases
# are REPORTED ONLY -- see the docstring section "Two records at one locus".
SAME_POSITION_TWO_RECORDS_RECORD = (
    "chr1\t2\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\t./.:.:.:.\n"
    "chr1\t2\t.\tA\tT,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t./.:.:.:.\t0/1:20:0,20,0:100,0,100,100,100,100\n"
)
GATK_SAME_POSITION_TWO_RECORDS_ROWS = [
    ("chr1\t2\t.\tA\tG\t92.64\t.\t"
     "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63\t"
     "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100\t./."),
    ("chr1\t2\t.\tA\tT\t92.64\t.\t"
     "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;MLEAC=1;MLEAF=0.500;QD=4.63\t"
     "GT:AD:DP:GQ:PL\t./.\t0/1:0,20:20:99:100,0,100"),
]
# The same two rows with --annotate-with-num-discovered-alleles: NDA=1 on each,
# because each row is genotyped from ONE input record whose merged ALT set holds
# exactly one allele (GenotypingEngine.java:464-465 measures `vc`, and `vc` here
# is that single record -- the merger never runs in this traversal mode).
GATK_SAME_POSITION_TWO_RECORDS_NDA_ROWS = [
    row.replace("MLEAF=0.500;", "MLEAF=0.500;NDA=1;")
    for row in GATK_SAME_POSITION_TWO_RECORDS_ROWS
]

CASES = [
    {
        "case": "default",
        "why": "the exact invocation of verify_genotype_gvcf.py:383-386 (no "
               "option flags): GATK's GenotypingEngine returns null for this "
               "monomorphic-after-pruning locus and writes no record, so the "
               "stale assertion at line 391 pinned native-only output",
        "body": STAR_RECORD,
        "args": [],
        "gated": True,
        "expect": [],
    },
    {
        "case": "include-non-variant-sites",
        "why": "same fixture with OutputMode.EMIT_ALL_ACTIVE_SITES: the "
               "REF-only no-call row survives and is the measured GATK contract "
               "for the audit's '-all-sites emits ALT=. GT=./.' claim",
        "body": STAR_RECORD,
        "args": ["--include-non-variant-sites"],
        "gated": True,
        "expect": [GATK_DEFAULT_ROW],
    },
    {
        "case": "stand-call-conf-50-include-non-variant-sites",
        "why": "the standard-confidence option cannot change the outcome here: "
               "the ALT set is already empty, so the confidence value is "
               "irrelevant to allele pruning and the row is unchanged",
        "body": STAR_RECORD,
        "args": ["--standard-min-confidence-threshold-for-calling", "50",
                 "--include-non-variant-sites"],
        "gated": True,
        "expect": [GATK_DEFAULT_ROW],
    },
    {
        "case": "surviving-concrete-alt",
        "why": "GATK keeps ALT='G' and calls 0/1.  The subset projection onto the "
               "kept alleles [A,G] collapses the source PL row [100,100,100] to a "
               "constant, so GATKVariantContextUtils.makeGenotypeCall() takes its "
               "PREFER_PLS fallback: isInformative is false "
               "(GATKVariantContextUtils.java:54-58; :333-338) and the call is "
               "projected from the SOURCE genotype with bestMatchToOriginalGT() "
               "(:397-403), which keeps every source allele that survived the "
               "subset and replaces a pruned one by the reference.  No GQ is "
               "assigned on that branch, and the published PL row is the "
               "min-shifted projection (:90-95, GenotypeLikelihoods.GLsToPLs) = "
               "0,0,0.  AC/AF/AN then follow from that genotype through the "
               "StandardAnnotation pass: ChromosomeCounts counts the called "
               "alleles of the published genotypes "
               "(ChromosomeCounts.java:43-53 -> "
               "VariantContextUtils.calculateChromosomeCounts) after "
               "GenotypeGVCFsEngine.regenotypeVC() has finished them "
               "(GenotypeGVCFsEngine.java:187-190)",
        "body": G_PLAUSIBLE_RECORD,
        "args": [],
        "gated": True,
        "expect": [GATK_G_ROW],
    },
    {
        "case": "surviving-concrete-alt-star-not-first",
        "why": "the same PREFER_PLS fallback with the pruned '*' in the middle of "
               "the ALT list: the projection is by allele identity, so the row is "
               "byte-identical to surviving-concrete-alt and pins that GATK does "
               "not depend on the position of the dropped allele",
        "body": G_STAR_SECOND_RECORD,
        "args": [],
        "gated": True,
        "expect": [GATK_G_ROW],
    },
    {
        "case": "surviving-concrete-alt-reversed-source-gt",
        "why": "the measured contract for the copy order: bestMatchToOriginalGT() "
               "maps the source allele LIST position by position "
               "(GATKVariantContextUtils.java:397-403), so a source call written "
               "2/0 yields 1/0 and not a canonicalized 0/1.  This case separates "
               "'project the source genotype' from 'derive the genotype from the "
               "PLs' (whose argmax genotype would be (*,G) here)",
        "body": G_REVERSED_GT_RECORD,
        "args": [],
        "gated": True,
        "expect": [GATK_G_REVERSED_ROW],
    },
    {
        "case": "surviving-concrete-alt-phased-source-gt",
        "why": "GenotypeBuilder retains the source genotype's phase flag "
               "(GATKVariantContextUtils.java:397-403 rebuilds the allele list "
               "over the same builder), so a phased source call stays phased "
               "after the projection",
        "body": G_PHASED_GT_RECORD,
        "args": [],
        "gated": True,
        "expect": [GATK_G_PHASED_ROW],
    },
    {
        "case": "surviving-deletion-alt",
        "why": "the locus' OWN deletion allele must not own its own '*'.  GATK "
               "evaluates the spanning-deletion test with "
               "isVcCoveredByDeletion() (GenotypingEngine.java:365-371), which "
               "requires loc.getStart() < vc.getStart(), and it only calls "
               "recordDeletions() for the emitted alleles AFTER the output "
               "subset has been computed (:178-179), so a deletion allele of "
               "the record being genotyped can never cover it.  The '*' is "
               "spurious (:314), is dropped from the ALT subset (:316), and "
               "the surviving 'A' is published with the source call projected "
               "by bestMatchToOriginalGT() (0/2 -> 0/1) and the min-shifted PL "
               "row 0,0,0",
        "body": G_DELETION_ALT_RECORD,
        "args": [],
        "gated": True,
        "expect": [GATK_DEL_ROW],
    },
    {
        "case": "surviving-deletion-alt-noncovering-upstream",
        "why": "the same decision with an explicit upstream deletion record: "
               "the deletion emitted at 2 spans 2-3 (recordDeletions, "
               "GenotypingEngine.java:343-357) and therefore does NOT cover "
               "locus 4 (:365-371), so '*' stays spurious and the locus-4 row "
               "is byte-identical to surviving-deletion-alt",
        "body": DEL_ALT_NONCOVERING_UPSTREAM,
        "args": [],
        "gated": True,
        "expect": [GATK_DEL_UPSTREAM_LOCUS_ROW, GATK_DEL_DOWNSTREAM_ROW],
    },
    {
        "case": "surviving-deletion-alt-covering-upstream",
        "why": "positive control for the same boundary: the upstream deletion "
               "emitted at 2 now spans 2-4 (REF AAA -> ALT A, deletionSize 2, "
               "GenotypingEngine.java:348-355), so locus 4 IS covered, '*' "
               "survives the subset and the site is genotyped 1/2.  Without "
               "this case a fix that pruned every '*' would look correct",
        "body": DEL_ALT_COVERING_UPSTREAM,
        "args": [],
        "gated": True,
        "expect": [GATK_DEL_LONG_DELETION_ROW, GATK_STAR_AND_DEL_ROW],
    },
    {
        "case": "surviving-deletion-alt-only-alt-implausible",
        "why": "pruning the orphan '*' leaves 'A' as the only ALT and it fails "
               "the standard-confidence threshold, so the site is monomorphic "
               "(:311-318) and the default traversal writes nothing "
               "(:167-169, GenotypeGVCFsEngine.java:188-190)",
        "body": SELF_DELETION_ONLY_RECORD,
        "args": [],
        "gated": True,
        "expect": [],
    },
    {
        "case": "covered-star-only-record",
        "why": "the record's ONLY surviving allele is the symbolic spanning "
               "deletion.  Ownership is correct -- the deletion emitted at 2 "
               "spans 2-3 and covers locus 3 (GenotypingEngine.java:365-371) -- "
               "but GATK still publishes nothing, because when the output "
               "allele set is exactly [SPAN_DEL] the engine returns null unless "
               "the traversal is EMIT_ALL_ACTIVE_SITES "
               "(GenotypingEngine.java:173-175: 'return a null call if we "
               "aren't forcing site emission and the only alt allele is a "
               "spanning deletion'), and GenotypeGVCFs.apply() additionally "
               "refuses such a record (GenotypeGVCFs.java:327-328 with "
               "GATKVariantContextUtils.isSpanningDeletionOnly() at "
               "GATKVariantContextUtils.java:2089-2091).  Only the upstream "
               "deletion row is written",
        "body": STAR_ONLY_COVERED_RECORD,
        "args": [],
        "gated": True,
        "expect": [GATK_DEL_UPSTREAM_LOCUS_ROW],
    },
    {
        "case": "covered-star-only-record-non-pass",
        "why": "the same rule with the dropped record non-PASS (FILTER=LowQual) "
               "and shallower (DP=7): the measured GATK output is identical "
               "apart from the surviving locus, which pins that the "
               "star-only drop is a property of the output allele set and not "
               "of the record's FILTER/INFO/FORMAT content",
        "body": STAR_ONLY_COVERED_NONPASS_RECORD,
        "args": [],
        "gated": True,
        "expect": [GATK_DEL_UPSTREAM_LOCUS_ROW],
    },
    {
        "case": "covered-star-plus-concrete-alt",
        "why": "positive control for the same rule: the locus also carries the "
               "concrete ALT G and the source row's best cell is (*,G), so the "
               "'*' is plausible AND covered.  Two ALTs survive the output "
               "subset, isSpanningDeletionOnly() is false "
               "(GATKVariantContextUtils.java:2089-2091), and GATK publishes "
               "ALT='*,G' with the PL-argmax call 1/2.  A fix that dropped "
               "every record carrying a '*' would break this case",
        "body": COVERED_STAR_PLUS_CONCRETE_RECORD,
        "args": [],
        "gated": True,
        "expect": [GATK_DEL_UPSTREAM_LOCUS_ROW, GATK_STAR_PLUS_CONCRETE_ROW],
    },
    {
        "case": "covered-star-implausible-plus-concrete-alt",
        "why": "the same locus with the source row's best cell moved to (A,G): "
               "the '*' is still covered but is individually implausible, and "
               "GenotypingEngine.calculateOutputAlleleSubset() prunes every "
               "ALT unless it passesThreshold() (GenotypingEngine.java:312-316) "
               "-- the AF threshold applies to '*' exactly as it does to a "
               "concrete ALT.  Only the concrete ALT is published, so the rule "
               "is 'drop the unqualified ALLELE', not 'drop the record'",
        "body": COVERED_STAR_IMPLAUSIBLE_PLUS_CONCRETE_RECORD,
        "args": [],
        "gated": True,
        "expect": [GATK_DEL_UPSTREAM_LOCUS_ROW, GATK_CONCRETE_ONLY_ROW],
    },
    {
        "case": "covered-star-only-record-dense",
        "why": "dense mode (OutputMode.EMIT_ALL_ACTIVE_SITES, "
               "GenotypeGVCFsEngine.java:381-383) skips both null returns, so "
               "the locus is materialized as GATK's REF-only no-call row "
               "(subsetToRefOnly at GenotypingEngine.java:190 and "
               "cleanupGenotypeAnnotations at GenotypeGVCFsEngine.java:191-194) "
               "and apply() writes it because forceOutput is true "
               "(GenotypeGVCFs.java:326-330; the ALT set is empty, so the "
               "record is not spanning-deletion-only).  Its QUAL is the literal "
               "token Infinity because the site is monomorphic, so the "
               "confidence branch is GenotypingEngine.java:158-163 with "
               "AFresult.log10ProbVariantPresent(), and "
               "MathUtils.log10OneMinusPow10(0.0) is "
               "Double.NEGATIVE_INFINITY -- a genuine positive-infinity QUAL "
               "that htsjdk renders with Java's `%.2f` (see the "
               "GATK_STAR_ONLY_DENSE_ROW comment)",
        "body": STAR_ONLY_COVERED_RECORD,
        "args": ["--include-non-variant-sites"],
        "gated": True,
        "expect": [GATK_DEL_UPSTREAM_LOCUS_ROW, GATK_STAR_ONLY_DENSE_ROW],
    },
    {
        "case": "covered-star-only-record-dense-two-samples",
        "why": "the same locus with two samples in one gVCF: the allele rule "
               "and the QUAL branch are per-locus, the upstream row is "
               "genotyped from both samples, and the materialized REF-only row "
               "repeats GATK's Infinity token once per sample.  This separates "
               "'the token belongs to the monomorphic-confidence branch' from "
               "'the token is an artefact of a single-sample path'",
        "body": STAR_ONLY_COVERED_TWO_SAMPLES_RECORD,
        "header": HEADER_TWO_SAMPLES,
        "args": ["--include-non-variant-sites"],
        "gated": True,
        "expect": [GATK_DEL_UPSTREAM_TWO_SAMPLES_ROW,
                   GATK_STAR_ONLY_DENSE_TWO_SAMPLES_ROW],
    },
    {
        "case": "covered-star-only-record-dense-non-pass",
        "why": "THE FILTER CASE.  Dense mode over the non-PASS variant of the "
               "same fixture: the materialized REF-only row is "
               "`chr1 3 . A . Infinity . DP=7;MLEAC=.;MLEAF=. GT ./.` -- GATK "
               "writes `.`, not the source record's LowQual.  "
               "GenotypingEngine.calculateGenotypes() rebuilds the call with "
               "`new VariantContextBuilder(callSourceString(), vc.getContig(), "
               "vc.getStart(), vc.getEnd(), outputAlleles)` "
               "(GenotypingEngine.java:181), which copies no filter state from "
               "the source VariantContext, so `filtersWereApplied` is false and "
               "htsjdk's VCFEncoder emits the unfiltered token `.`; the ONLY "
               "filter GATK can apply is its own threshold decision at "
               "`if (!passesCallThreshold(phredScaledConfidence)) "
               "builder.filter(LOW_QUAL_FILTER_NAME)` (:184-186), and this "
               "locus' recomputed phred confidence is the infinite QUAL, which "
               "passes.  A source FILTER must therefore never be inherited by "
               "the output record",
        "body": STAR_ONLY_COVERED_NONPASS_RECORD,
        "args": ["--include-non-variant-sites"],
        "gated": True,
        "expect": [GATK_DEL_UPSTREAM_LOCUS_ROW,
                   GATK_STAR_ONLY_DENSE_NONPASS_ROW],
    },
    {
        "case": "covered-star-only-record-dense-non-pass-two-samples",
        "why": "the same FILTER rule with two samples, so the gate does not "
               "depend on a single-sample merge path: the upstream row is "
               "genotyped from both samples and the materialized REF-only row "
               "carries one './.' per sample with GATK's `.` FILTER",
        "body": STAR_ONLY_COVERED_NONPASS_TWO_SAMPLES_RECORD,
        "header": HEADER_TWO_SAMPLES,
        "args": ["--include-non-variant-sites"],
        "gated": True,
        "expect": [GATK_DEL_UPSTREAM_TWO_SAMPLES_ROW,
                   GATK_STAR_ONLY_DENSE_NONPASS_TWO_SAMPLES_ROW],
    },
    {
        "case": "covered-star-only-record-dense-non-pass-stream-by-locus",
        "why": "the same fixture through native's OTHER writer "
               "(--stream-by-locus).  Both writers funnel through the same "
               "GATK-compatibility transform, so a FILTER fix that only covers "
               "the aggregate path would leave this case red",
        "body": STAR_ONLY_COVERED_NONPASS_RECORD,
        "args": ["--include-non-variant-sites"],
        "native_args": ["--stream-by-locus"],
        "gated": True,
        "expect": [GATK_DEL_UPSTREAM_LOCUS_ROW,
                   GATK_STAR_ONLY_DENSE_NONPASS_ROW],
    },
    {
        "case": "header-filter-line-with-extra-declared-filter",
        "why": "THE HEADER CASE.  The input declares `##FILTER=<ID=q10,...>` and "
               "NOT LowQual.  GATK's output header still carries the standard "
               "`##FILTER=<ID=LowQual,Description=\"Low quality\">` because "
               "setupVCFWriter() adds it unconditionally as its last header "
               "line (GenotypeGVCFsEngine.java:416, text from "
               "GATKVCFHeaderLines.java:89 + GATKVCFConstants.java:179), and "
               "htsjdk's sorted VCFHeader emission places it BEFORE q10.  The "
               "gate compares the whole ##FILTER group byte for byte, so this "
               "case fails both when the LowQual declaration is missing and "
               "when it is appended at the wrong place in the group",
        "body": STAR_ONLY_COVERED_RECORD,
        "header": HEADER_WITH_Q10_FILTER,
        "args": [],
        "gated": True,
        "expect": [GATK_DEL_UPSTREAM_LOCUS_ROW],
    },
    {
        "case": "header-filter-line-with-input-pass-declaration",
        "why": "THE PASS CASE.  The input declares the reserved "
               "`##FILTER=<ID=PASS,Description=\"All filters passed\">`.  GATK "
               "KEEPS it: the writer's header starts from the input's own "
               "metadata (`new LinkedHashSet<>(inputVCFHeader."
               "getMetaDataInInputOrder())`, GenotypeGVCFsEngine.java:395) and "
               "the only removal is the GVCFBlock strip at :398-399, so the line "
               "is written back verbatim; the only filter GATK adds is LowQual "
               "(:416).  htsjdk does not synthesize a PASS line either (the byte "
               "string `All filters passed` occurs in neither the pinned "
               "gatk-package jar nor htsjdk 4.2.0), so the text proves "
               "propagation rather than generation.  Measured group, in order: "
               "[LowQual, PASS].",
        "body": STAR_ONLY_COVERED_RECORD,
        "header": HEADER_WITH_PASS_FILTER,
        "args": [],
        "gated": True,
        "expect": [GATK_DEL_UPSTREAM_LOCUS_ROW],
    },
    {
        "case": "header-filter-line-with-nonstandard-pass-description",
        "why": "the same case with a non-standard Description.  GATK keeps the "
               "input's text byte for byte -- it does not replace or normalize a "
               "PASS declaration, because it never generates one at all "
               "(GenotypeGVCFsEngine.java:395 + :416).  This case separates "
               "'propagate the input's PASS line' from 'emit a canonical PASS "
               "line', which the standard-description case alone could not "
               "distinguish.  Measured group, in order: "
               "[LowQual, PASS(\\\"Some other PASS text\\\")].",
        "body": STAR_ONLY_COVERED_RECORD,
        "header": HEADER_WITH_NONSTANDARD_PASS_FILTER,
        "args": [],
        "gated": True,
        "expect": [GATK_DEL_UPSTREAM_LOCUS_ROW],
    },
    {
        "case": "header-filter-line-with-pass-and-other-filter",
        "why": "PASS together with an unrelated filter.  The whole group is "
               "propagated and htsjdk emits its header in sorted order, so the "
               "measured GATK group is [LowQual, PASS, q10] at indices 2-4: "
               "LowQual first because the tool adds it (:416) and 'L' < 'P' < "
               "'q', and PASS before q10 for the same reason.  The order is the "
               "assertion here, not just the set -- dropping only the PASS line "
               "leaves [LowQual, q10], which this case must reject.",
        "body": STAR_ONLY_COVERED_RECORD,
        "header": HEADER_WITH_PASS_AND_Q10_FILTER,
        "args": [],
        "gated": True,
        "expect": [GATK_DEL_UPSTREAM_LOCUS_ROW],
    },
    {
        "case": "header-content-with-unquoted-descriptions",
        "why": "THE HEADER CONTENT CASE.  The input header writes every "
               "Description unquoted, the shape verify_genotype_gvcf.py's "
               "fixture has.  htsjdk parses each `##INFO`/`##FORMAT`/`##ALT` "
               "line into a VCFCompoundHeaderLine and writes it back out of its "
               "own fields, so GATK's output quotes every one of them "
               "(measured: `Description=\"Read depth\"` for an input that said "
               "`Description=Read depth`), and on top of the input's lines the "
               "tool declares its own annotation lines "
               "(GenotypeGVCFsEngine.java:395-416): the three rank-sum INFO "
               "lines, the standard INFO/DP line (:407, *\"needed for gVCFs "
               "without DP tags\"*) and the standard FORMAT/AD line "
               "(DepthPerAlleleBySample via VariantAnnotation.java:22-33).  Both "
               "halves are asserted as one multiset comparison, so this case "
               "fails while either half is wrong",
        "body": STAR_ONLY_COVERED_RECORD,
        "args": [],
        "gated": True,
        "expect": [GATK_DEL_UPSTREAM_LOCUS_ROW],
    },
    {
        "case": "header-content-with-quoted-descriptions",
        "why": "the same case with every Description ALREADY quoted, i.e. the "
               "shape a GATK-produced gVCF has.  Nothing needs re-quoting here, "
               "so only the missing declarations can make this case fail -- it "
               "separates 'quote the input's lines' from 'declare GATK's own "
               "lines', which the unquoted case alone could not distinguish",
        "body": STAR_ONLY_COVERED_RECORD,
        "header": HEADER_QUOTED_DESCRIPTIONS,
        "args": [],
        "gated": True,
        "expect": [GATK_DEL_UPSTREAM_LOCUS_ROW],
    },
    {
        "case": "header-content-with-gatk-mle-wording",
        "why": "the input already declares MLEAC/MLEAF with GATK's own long "
               "wording, as a HaplotypeCaller gVCF does.  GATK's header set "
               "de-duplicates by the whole line, so its identical pair is not "
               "added twice and native must produce exactly one of each: this "
               "case fails if native over-declares them (and it also pins the "
               "wording half, because native used to declare these two keys with "
               "a short description of its own, "
               "`Description=Maximum likelihood allele count`)",
        "body": STAR_ONLY_COVERED_RECORD,
        "header": HEADER_WITH_GATK_MLE_LINES,
        "args": [],
        "gated": True,
        "expect": [GATK_DEL_UPSTREAM_LOCUS_ROW],
    },
    {
        "case": "non-pass-variant-input-undeclared-filter",
        "why": "the ORDINARY (non-dense) variant path with a non-PASS source "
               "FILTER, and the source FILTER declared nowhere in the header.  "
               "GATK writes the byte-identical GATK_G_ROW (FILTER `.`): the "
               "engine rebuilds the record (GenotypingEngine.java:181) and the "
               "recomputed confidence 82.26 passes :430-431, so no filter is "
               "applied.  This case also pins the wrong-NAME half of the "
               "defect: the undeclared source FILTER is registered by HTSlib "
               "while a record is parsed (vcf_parse_filter appends "
               "`##FILTER=<ID=...,Description=\"Dummy\">`), i.e. AFTER the "
               "output header was duplicated from the header-only first pass, "
               "so the inherited index names an unrelated id of the output "
               "dictionary",
        "body": G_PLAUSIBLE_NONPASS_RECORD,
        "args": [],
        "gated": True,
        "expect": [GATK_G_NONPASS_ROW],
    },
    {
        "case": "non-pass-variant-input-declared-filter",
        "why": "the same record with `##FILTER=<ID=LowQual,...>` declared in "
               "the input header, so the source index resolves to the right "
               "name.  GATK STILL writes `.` (GATK_G_ROW), which separates the "
               "two layers of the defect: the undeclared header only decides "
               "whether the inherited filter is spelled `LowQual` or something "
               "unrelated; carrying the source FILTER at all is wrong either "
               "way",
        "body": G_PLAUSIBLE_NONPASS_RECORD,
        "header": HEADER_WITH_LOWQUAL_FILTER,
        "args": [],
        "gated": True,
        "expect": [GATK_G_NONPASS_ROW],
    },
    {
        "case": "weak-locus-lowqual-filter-not-implemented",
        "why": "REPORTED ONLY, different root cause.  A dense-mode locus whose "
               "recomputed phred confidence is 23.14 < "
               "--standard-min-confidence-threshold-for-calling (30): GATK "
               "applies its OWN filter at GenotypingEngine.java:184-186 "
               "(`builder.filter(LOW_QUAL_FILTER_NAME)`) and writes "
               "`FILTER=LowQual`.  Native does not implement that decision at "
               "all, so it writes `.`.  This is the positive half of the "
               "FILTER contract; it is not caused by the inherited-filter "
               "defect this gate was added for, and it stays reported-only "
               "until the threshold rule is implemented",
        "body": WEAK_LOCUS_RECORD,
        "args": ["--include-non-variant-sites"],
        "gated": False,
        "expect": [GATK_WEAK_LOCUS_LOWQUAL_ROW],
    },
    {
        "case": "unemitted-upstream-deletion-star-plus-concrete-alt",
        "why": "reported only: the upstream deletion record is IMPLAUSIBLE, so "
               "GATK drops it (GenotypingEngine.java:167-169) and "
               "recordDeletions() is never called for it (:178-179); the "
               "downstream '*' is therefore a spurious spanning deletion "
               "(:314) and is pruned even though it is plausible, leaving "
               "ALT='G' and the projected 0/1 call.  Native builds its deletion "
               "spans from the INPUT records, so it still treats the '*' as "
               "owned and publishes '*,G' with 1/2.  This is the residual "
               "'input records vs emitted alleles' structural difference "
               "flagged in .diag/round-star-ownership.md section 8, now visible "
               "in the default output mode because the star-only record rule no "
               "longer hides it; fixing it requires the ordered "
               "emitted-deletions state and is out of scope for this round",
        "body": IMPLAUSIBLE_UPSTREAM_STAR_PLUS_CONCRETE_RECORD,
        "args": [],
        "gated": False,
        "expect": [GATK_G_DOWNSTREAM_ROW],
    },
    {
        "case": "nda-annotation-single-alt",
        "why": "the --annotate-with-num-discovered-alleles data-row contract on "
               "the shape native's `record.alleles.size() < 3` guard skipped: "
               "one concrete ALT besides <NON_REF>.  GATK's value comes from "
               "`vc.getAlternateAlleles().size()` on the MERGED input record "
               "(GenotypingEngine.java:464-465), i.e. 1 here because the merger "
               "has already dropped <NON_REF> "
               "(ReferenceConfidenceVariantContextMerger.java:339-345 via "
               "GenotypeGVCFsEngine.java:136), while native wrote no NDA at all",
        "body": NDA_SINGLE_ALT_RECORD,
        "args": [NDA_FLAG],
        "gated": True,
        "expect": [GATK_NDA_SINGLE_ALT_ROW],
    },
    {
        "case": "nda-annotation-control-without-flag",
        "why": "the control for the same fixture with NO flag: both the `NDA` "
               "key and its `##INFO=<ID=NDA,...>` declaration are produced only "
               "under the flag (GenotypingEngine.java:464-465 for the key, "
               ":90-94 for the header line), so this case fails if native (or "
               "GATK) annotates or declares it unconditionally",
        "body": NDA_SINGLE_ALT_RECORD,
        "args": [],
        "gated": True,
        "expect": [GATK_NDA_CONTROL_ROW],
    },
    {
        "case": "nda-annotation-reduced-alt-set",
        "why": "the count follows the INPUT allele set, not the published one: "
               "G_PLAUSIBLE_RECORD ('*,G,<NON_REF>') is pruned to ALT='G' but "
               "GATK writes NDA=2, because composeCallAttributes() is handed the "
               "pre-pruning `vc` (GenotypingEngine.java:464-465) rather than "
               "outputAlternativeAlleles (:155)",
        "body": G_PLAUSIBLE_RECORD,
        "args": [NDA_FLAG],
        "gated": True,
        "expect": [GATK_NDA_REDUCED_ALT_SET_ROW],
    },
    {
        "case": "nda-annotation-covered-star-plus-concrete-alt",
        "why": "two published ALTs including a symbolic '*' that an upstream "
               "deletion owns: the '*' is a discovered allele, so the two-ALT "
               "row carries NDA=2 while the upstream one-ALT deletion row "
               "carries NDA=1 (ReferenceConfidenceVariantContextMerger.java:"
               "340-342 re-adds Allele.SPAN_DEL to the merged allele set)",
        "body": COVERED_STAR_PLUS_CONCRETE_RECORD,
        "args": [NDA_FLAG],
        "gated": True,
        "expect": GATK_NDA_COVERED_STAR_PLUS_CONCRETE_ROWS,
    },
    {
        "case": "nda-annotation-covered-star-only-dense",
        "why": "the REF-only dense materialization of a covered-'*' locus keeps "
               "the INPUT record's count: the merged set is [A, *] so GATK "
               "writes NDA=1 on the row whose ALT is '.' and whose FORMAT is "
               "GT only (GenotypeGVCFsEngine.java:191-194).  A locus whose "
               "merged ALT set is EMPTY writes no NDA at all -- that path is "
               "`regenotypeVC`'s non-variant branch (:148/:154) and never reaches "
               "composeCallAttributes()",
        "body": STAR_ONLY_COVERED_RECORD,
        "args": ["--include-non-variant-sites", NDA_FLAG],
        "gated": True,
        "expect": GATK_NDA_COVERED_STAR_ONLY_DENSE_ROWS,
    },
    {
        "case": "nda-annotation-two-samples",
        "why": "the count is a property of the locus, not of the sample count: "
               "the same single-ALT fixture with two samples still writes "
               "NDA=1, so a fix that summed per-sample discoveries would fail",
        "body": NDA_TWO_SAMPLES_RECORD,
        "header": HEADER_TWO_SAMPLES,
        "args": [NDA_FLAG],
        "gated": True,
        "expect": [GATK_NDA_TWO_SAMPLES_ROW],
    },
    {
        "case": "cross-sample-alt-union-single-record",
        "why": "the SUPPORTED cross-sample shape, gated: one record at chr1:2 "
               "whose ALT set already holds both samples' alleles (G for the "
               "0/1 sample, T for the 2/2 one), i.e. exactly what CombineGVCFs "
               "writes for two single-sample GVCFs (measured: `chr1 2 . A "
               "T,G,<NON_REF>`).  Both tools emit one row with ALT=G,T and the "
               "per-sample genotypes in it, byte for byte -- this is the "
               "boundary case that shows native's locus-level coalescing "
               "(genotype_gvcf_tool.cpp:1180-1187) is exact whenever the input "
               "has one record per locus, which is the only shape "
               "GenotypeGVCFs.java:65-66 documents and the one its group-by-"
               "locus traversal asserts (:359-364)",
        "body": CROSS_SAMPLE_ALT_UNION_RECORD,
        "header": HEADER_TWO_SAMPLES,
        "args": [],
        "gated": True,
        "expect": [GATK_CROSS_SAMPLE_ALT_UNION_ROW],
    },
    {
        "case": "same-position-two-records-cross-sample-alt",
        "why": "REPORTED ONLY -- the input shape GATK does not support.  Two "
               "records starting at chr1:2, one per sample, with DIFFERENT ALT "
               "sets.  GATK's default traversal is by-variant "
               "(GenotypeGVCFs.java:284-285 -> VariantLocusWalker.java:132-142), "
               "so nothing is merged across records: each record is genotyped "
               "alone and written as its own row, and the record that does not "
               "carry the other sample leaves it as a bare `./.`.  Native "
               "coalesces the two records (record_key = rid:pos:REF, "
               "genotype_gvcf_tool.cpp:1180-1187; aggregate group loop :6538-"
               "6650) and writes ONE row (ALT=G, second record's sample "
               "dropped by the overlapping-sample rule at :6594-6600).  The "
               "input is not producible by the documented pipeline: "
               "CombineGVCFs collapses these two records into one "
               "(`A T,G,<NON_REF>`, pinned by the gated case above) and "
               "GenotypeGVCFs takes a single input track (GenotypeGVCFs.java:"
               "65-66).  See .diag/round-crosssample-merge.md",
        "body": SAME_POSITION_TWO_RECORDS_RECORD,
        "header": HEADER_TWO_SAMPLES,
        "args": [],
        "gated": False,
        "expect": GATK_SAME_POSITION_TWO_RECORDS_ROWS,
    },
    {
        "case": "same-position-two-records-cross-sample-alt-nda",
        "why": "REPORTED ONLY -- the same unsupported shape with "
               "--annotate-with-num-discovered-alleles, which is how the shape "
               "was first measured (.diag/round-nda.md section 6, first "
               "bullet): GATK writes TWO rows, each carrying NDA=1, because "
               "the value GenotypingEngine.java:464-465 reads is the ALT count "
               "of the single record that row was genotyped from.  Native "
               "writes one row with NDA=2 (the coalesced union).  Same root "
               "cause and same scope as the case above",
        "body": SAME_POSITION_TWO_RECORDS_RECORD,
        "header": HEADER_TWO_SAMPLES,
        "args": [NDA_FLAG],
        "gated": False,
        "expect": GATK_SAME_POSITION_TWO_RECORDS_NDA_ROWS,
    },
    {
        "case": "same-position-two-records-dense-unsupported",
        "why": "REPORTED ONLY -- GATK's OTHER traversal mode rejects the shape "
               "outright.  --include-non-variant-sites (or any "
               "--force-output-intervals) turns off the by-variant switch "
               "(GenotypeGVCFs.java:284-285), so apply() receives every variant "
               "overlapping the locus; GenotypeGVCFsEngine.getVariantSubsetTo"
               "Process() then insists on at most ONE variant starting at the "
               "locus and throws `IllegalStateException: Variant input contains "
               "more than one variant starting at location: chr1:2-2` "
               "(:349-368, the comment at :359-360: 'since this tool only "
               "accepts a single input source, there should never be more than "
               "one variant at a given starting locus').  GATK exits 3 and "
               "writes no data row, so there is no GATK output contract for "
               "this input to match; native exits 0 with its merged row",
        "body": SAME_POSITION_TWO_RECORDS_RECORD,
        "header": HEADER_TWO_SAMPLES,
        "args": ["--include-non-variant-sites"],
        "gated": False,
        "gatk_expected_exit": 3,
        "expect": [],
    },
]


def data_rows(path: pathlib.Path) -> list[str]:
    opener = gzip.open if path.name.endswith(".gz") else open
    rows: list[str] = []
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line and not line.startswith("#"):
                rows.append(line.rstrip("\n"))
    return rows


def header_lines(path: pathlib.Path) -> list[str]:
    """Every header line of an output, in order, '#CHROM' included."""
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return [line.rstrip("\n") for line in stream if line.startswith("#")]


def filter_header_lines(lines: list[str]) -> list[str]:
    """The ##FILTER group, in order -- the surface this gate asserts."""
    return [line for line in lines if line.startswith("##FILTER=")]


def content_header_lines(lines: list[str]) -> list[str]:
    """Every header content line, in order, minus the producer-identity ones.

    ``##GATKCommandLine`` records the command line and version of the GATK
    process that produced the file and ``##source`` names the program that wrote
    it; native is a different program and is deliberately not required to
    reproduce either.  Everything else -- ``##fileformat``, ``##ALT``,
    ``##FILTER``, ``##FORMAT``, ``##INFO``, ``##contig`` -- is a declaration
    about the file's contents and must match byte for byte (order aside; the
    order-only residual is reported, not gated).
    """
    return [line for line in lines
            if line.startswith("##")
            and not line.startswith(("##GATKCommandLine=", "##source="))]


def write_reference(work: pathlib.Path) -> pathlib.Path:
    """A 100 bp chr1 whose bases match every fixture REF allele."""
    reference = work / "reference.fa"
    reference.write_text(">chr1\n" + "A" * 100 + "\n", encoding="utf-8")
    reference.with_name(reference.name + ".fai").write_text(
        "chr1\t100\t6\t100\t101\n", encoding="utf-8")
    # GATK resolves the sequence dictionary as <stem>.dict.
    reference.with_suffix(".dict").write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n", encoding="utf-8")
    return reference


def invoke(command: list[str], label: str, timeout: int) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True,
                            check=False, timeout=timeout)
    if result.returncode != 0:
        print(f"# {label} exited {result.returncode}")
        print(result.stderr[-4000:])
    return result


def run_case(case: dict, work: pathlib.Path, reference: pathlib.Path,
             java: pathlib.Path, jar: pathlib.Path, native: pathlib.Path,
             timeout: int) -> dict:
    source = work / f"{case['case']}.g.vcf"
    source.write_text(case.get("header", HEADER) + case["body"], encoding="utf-8")
    # GATK's --include-non-variant-sites switches to a group-by-locus traversal
    # that requires an index; a plain .vcf gets a tribble .idx.  Creating it for
    # every case keeps the two tools' inputs byte-identical.
    index_result = invoke([str(java), "-Xmx1g", "-jar", str(jar),
                           "IndexFeatureFile", "-I", str(source)],
                          f"GATK IndexFeatureFile [{case['case']}]", timeout)
    gatk_out = work / f"{case['case']}-gatk.vcf"
    native_out = work / f"{case['case']}-native.vcf"
    common = ["-R", str(reference), "-V", str(source), *case["args"]]
    gatk_result = invoke([str(java), "-Xmx1g", "-jar", str(jar),
                          "GenotypeGVCFs", *common, "-O", str(gatk_out),
                          "--create-output-variant-index", "false"],
                         f"GATK GenotypeGVCFs [{case['case']}]", timeout)
    # Native's default profile is a non-GATK diagnostic output, so the
    # GATK-compatibility writer is the only byte-comparable surface here; the
    # semantic arguments above are identical on both sides.
    native_result = invoke([str(native), *common,
                            "--gatk-compatible-annotations",
                            *case.get("native_args", []),
                            "-O", str(native_out)],
                           f"native GenotypeGVCFs [{case['case']}]", timeout)

    gatk_rows = data_rows(gatk_out) if gatk_out.exists() else []
    native_rows = data_rows(native_out) if native_out.exists() else []
    gatk_headers = header_lines(gatk_out) if gatk_out.exists() else []
    native_headers = header_lines(native_out) if native_out.exists() else []
    gatk_filters = filter_header_lines(gatk_headers)
    native_filters = filter_header_lines(native_headers)
    # The declarations themselves: every ## line but the two producer-identity
    # ones, compared as a multiset (duplicates included -- GATK really does
    # declare two ##INFO=<ID=DP,...> and two ##FORMAT=<ID=AD,...> lines here).
    gatk_content = content_header_lines(gatk_headers)
    native_content = content_header_lines(native_headers)
    gatk_content_counter = collections.Counter(gatk_content)
    native_content_counter = collections.Counter(native_content)
    # Everything that is not the asserted ##FILTER group is recorded, not gated.
    # ##GATKCommandLine is excluded from the diff because it embeds the run's
    # absolute paths and timestamp and native never emits it (native has no
    # --add-output-vcf-command-line equivalent on this surface).
    def observed(lines: list[str]) -> list[str]:
        return [line for line in lines if not line.startswith(("##GATKCommandLine=",))]

    gatk_observed = observed(gatk_headers)
    native_observed = observed(native_headers)

    result = {
        "case": case["case"],
        "why": case["why"],
        "gated": case["gated"],
        "args": common,
        "native_args": case.get("native_args", []),
        "gatk_exit": gatk_result.returncode,
        "native_exit": native_result.returncode,
        "expect": case["expect"],
        "gatk_expected_exit": case.get("gatk_expected_exit", 0),
        "gatk_rows": gatk_rows,
        "native_rows": native_rows,
        "gatk_filter_header_lines": gatk_filters,
        "native_filter_header_lines": native_filters,
        "gatk_content_header_lines": gatk_content,
        "native_content_header_lines": native_content,
        "header_content_missing_in_native":
            sorted((gatk_content_counter - native_content_counter).elements()),
        "header_content_extra_in_native":
            sorted((native_content_counter - gatk_content_counter).elements()),
        "header_observations": {
            "only_in_gatk": [line for line in gatk_observed
                             if line not in native_observed],
            "only_in_native": [line for line in native_observed
                               if line not in gatk_observed],
            "order_differs": (gatk_observed != native_observed
                              and [line for line in gatk_observed
                                   if line not in native_observed] == []
                              and [line for line in native_observed
                                   if line not in gatk_observed] == []),
        },
        "violations": [],
    }
    if index_result.returncode != 0:
        result["violations"].append(f"GATK IndexFeatureFile exited {index_result.returncode}")
    # A case may pin a GATK run that is *expected* to fail (the group-by-locus
    # rejection of an input with two records starting at one locus); the exit
    # code is then an observation, not a violation.
    if (gatk_result.returncode != 0
            and gatk_result.returncode != case.get("gatk_expected_exit")):
        result["violations"].append(f"GATK exited {gatk_result.returncode}")
    if native_result.returncode != 0:
        result["violations"].append(f"native exited {native_result.returncode}")
    if case["expect"] is not None and gatk_rows != case["expect"]:
        result["violations"].append(
            "GATK moved away from the measured truth: "
            f"GATK={gatk_rows} expected={case['expect']}")
    if case["gated"]:
        if len(gatk_rows) != len(native_rows):
            result["violations"].append(
                f"record count differs: GATK={len(gatk_rows)} native={len(native_rows)}")
        for index, (gatk_row, native_row) in enumerate(zip(gatk_rows, native_rows)):
            if gatk_row != native_row:
                result["violations"].append(
                    f"row {index} is not byte-identical: GATK={gatk_row!r} "
                    f"NATIVE={native_row!r}")
        # The header contract: GATK declares the filter it may apply
        # (GenotypeGVCFsEngine.java:416 + GATKVCFHeaderLines.java:89); the
        # ##FILTER group must match byte for byte, in order.
        if gatk_filters != native_filters:
            result["violations"].append(
                "##FILTER header lines are not byte-identical: "
                f"GATK={gatk_filters} NATIVE={native_filters}")
        # GATK must still declare what this gate was written from; if the
        # pinned GATK stops emitting one of these lines, the comparison below
        # could pass for the wrong reason.
        missing_from_gatk = [line for line in GATK_DECLARED_HEADER_LINES
                             if line not in gatk_content_counter]
        if missing_from_gatk:
            result["violations"].append(
                "GATK moved away from the measured header truth: it no longer "
                f"declares {missing_from_gatk}")
        # The header's content lines: same multiset (order aside, and excluding
        # the producer-identity lines).  GATK re-quotes every Description and
        # declares its own annotations on top of the input's lines
        # (GenotypeGVCFsEngine.java:395-416), so both the re-quoted input lines
        # and the added declarations are asserted here.
        if gatk_content_counter != native_content_counter:
            result["violations"].append(
                "header content lines are not byte-identical (as a set): "
                "missing in native="
                f"{sorted((gatk_content_counter - native_content_counter).elements())} "
                "extra in native="
                f"{sorted((native_content_counter - gatk_content_counter).elements())}")
        # The NDA declaration is flag-conditional on BOTH sides
        # (GenotypingEngine.java:90-94 for GATK, the
        # `options.annotate_with_num_discovered_alleles` guard in
        # genotype_gvcf_tool.cpp for native), so it is asserted in both
        # directions rather than left to the set comparison above.
        if NDA_FLAG in case["args"]:
            if GATK_NDA_HEADER_LINE not in gatk_content_counter:
                result["violations"].append(
                    "GATK moved away from the measured header truth: it no "
                    f"longer declares {GATK_NDA_HEADER_LINE!r} under {NDA_FLAG}")
            if GATK_NDA_HEADER_LINE not in native_content_counter:
                result["violations"].append(
                    f"native does not declare {GATK_NDA_HEADER_LINE!r} under "
                    f"{NDA_FLAG}")
        elif GATK_NDA_HEADER_LINE in native_content_counter:
            result["violations"].append(
                f"native declares {GATK_NDA_HEADER_LINE!r} without {NDA_FLAG}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict oracle for GenotypeGVCFs on a gVCF record whose ALT "
                    "list carries the symbolic spanning deletion '*', against "
                    "pinned GATK 4.6.2.0.")
    parser.add_argument(
        "--native", default=os.environ.get("FASTGATK_GENOTYPE_BINARY"),
        help="native GenotypeGVCFs binary (default: $FASTGATK_GENOTYPE_BINARY "
             "or fastgatk-native/build/fastgatk-genotype-gvcf)")
    parser.add_argument(
        "--expect-divergence", action="store_true",
        help="diagnostic mode: report the measured GATK/native rows and exit 0 "
             "instead of gating parity (use this while the defect is unfixed)")
    parser.add_argument("--case", action="append", default=None,
                        help="run only the named case(s)")
    parser.add_argument("--timeout", type=int, default=300,
                        help="per-process timeout in seconds (default 300)")
    arguments = parser.parse_args()

    root = pathlib.Path(__file__).resolve().parents[2]
    native = pathlib.Path(arguments.native) if arguments.native else (
        pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD",
                                    root / "fastgatk-native" / "build"))
        / "fastgatk-genotype-gvcf")
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = pathlib.Path(os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java")))

    assets = [native, java, jar]
    if not all(path.is_file() for path in assets):
        oracle_guard.oracle_not_verified(
            "verify_genotype_gvcf_spandel_gatk_oracle.py", java, jar)
        missing = sorted(str(path) for path in assets if not path.is_file())
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing oracle assets: {missing}")
        print(json.dumps({"status": "skipped",
                          "reason": "missing native binary or pinned GATK assets",
                          "missing": missing}, sort_keys=True))
        return 0

    selected = [case for case in CASES
                if arguments.case is None or case["case"] in arguments.case]
    if not selected:
        raise SystemExit(f"no such case: {arguments.case}")
    strict_mode = not arguments.expect_divergence

    results: list[dict] = []
    with tempfile.TemporaryDirectory(
            prefix="fastgatk-genotype-spandel-oracle-") as directory:
        work = pathlib.Path(directory)
        reference = write_reference(work)
        for case in selected:
            results.append(run_case(case, work, reference, java, jar, native,
                                    arguments.timeout))

    violations: list[str] = []
    for result in results:
        if result["gated"]:
            violations.extend(f"[{result['case']}] {item}" for item in result["violations"])

    print(f"# {pathlib.Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# GATK rule under test: an ALT survives only when it passes "
          "standardConfidenceForCalling and, for '*', when an emitted deletion "
          "covers the locus (GenotypingEngine.java:311-327); with an empty ALT "
          "set the default traversal writes nothing (:167-169) while "
          "--include-non-variant-sites forces the REF-only no-call row "
          "(GenotypeGVCFsEngine.java:191-194, :381-383; cleanup at :479-491).")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; cases={[case['case'] for case in selected]}")
    for result in results:
        marker = "gated" if result["gated"] else "REPORTED ONLY (not gated)"
        print(f"[{result['case']}] {marker}")
        print(f"    why: {result['why']}")
        print(f"    args: {' '.join(result['args']) if result['args'] else '(none)'}")
        if result["native_args"]:
            print(f"    native-only args: {' '.join(result['native_args'])}")
        print(f"    gatk_exit={result['gatk_exit']} native_exit={result['native_exit']}")
        if result["expect"] is not None:
            print(f"    expected rows: {result['expect']}")
        print(f"    GATK   rows ({len(result['gatk_rows'])}):")
        for row in result["gatk_rows"]:
            print(f"        {row}")
        print(f"    NATIVE rows ({len(result['native_rows'])}):")
        for row in result["native_rows"]:
            print(f"        {row}")
        print(f"    GATK   ##FILTER header lines: {result['gatk_filter_header_lines']}")
        print(f"    NATIVE ##FILTER header lines: {result['native_filter_header_lines']}")
        print(f"    header content lines: GATK={len(result['gatk_content_header_lines'])} "
              f"NATIVE={len(result['native_content_header_lines'])}")
        for line in result["header_content_missing_in_native"]:
            print(f"        missing in native: {line}")
        for line in result["header_content_extra_in_native"]:
            print(f"        extra in native  : {line}")
        observation = result["header_observations"]
        print("    header observations (NOT gated): "
              f"only_in_gatk={len(observation['only_in_gatk'])} "
              f"only_in_native={len(observation['only_in_native'])} "
              f"order_differs={observation['order_differs']}")
        for line in observation["only_in_gatk"]:
            print(f"        GATK-only  : {line}")
        for line in observation["only_in_native"]:
            print(f"        NATIVE-only: {line}")
        for violation in result["violations"]:
            print(f"    VIOLATION: {violation}")

    payload = {
        "status": "pass" if not violations else "divergence",
        "gatk_version": "4.6.2.0",
        "mode": "strict" if strict_mode else "expect-divergence",
        "cases": [
            {key: value for key, value in result.items() if key != "gated"}
            for result in results
        ],
        "violations": violations,
        "rows_compared": "data rows byte-identical (CHROM..sample columns)",
        "header_compared": "##FILTER header lines byte-identical and in order; "
                           "all other header content lines byte-identical as a "
                           "set (order-only differences reported, not gated); "
                           "##GATKCommandLine and ##source excluded as "
                           "producer identity",
        "header_observations": "order-only header differences are reported, not gated",
    }
    print(json.dumps(payload, sort_keys=True))

    if not strict_mode:
        return 0
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
