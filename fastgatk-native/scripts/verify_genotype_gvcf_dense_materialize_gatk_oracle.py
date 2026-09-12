#!/usr/bin/env python3
"""Strict pinned-GATK oracle for the two dense-mode
(``--include-non-variant-sites``) record-materialization behaviours of
GenotypeGVCFs that native does not implement:

1. **spanning-deletion row synthesis** -- for a position that an input record
   SPANS (a record whose span covers the position but which does not start
   there), GATK's per-locus traversal still visits the position, the
   reference-confidence merger turns that record's deletion allele into the
   symbolic ``*``, and dense mode forces the locus out as a record.  Native
   never visits the position at all, so the row is missing.
2. **reverse trimming** of the merged record -- ``reverseTrimAlleles()`` clips
   the trailing bases every non-symbolic allele shares (and, when the whole
   allele would be clipped away, keeps exactly one base).  Native emits the
   untrimmed alleles.

Measured rules (pinned GATK 4.6.2.0, see ``.diag/round-dense-materialize.md``)
-----------------------------------------------------------------------------
*Traversal* -- ``GenotypeGVCFs`` extends ``VariantLocusWalker`` and uses the
default ``getDrivingVariantCacheLookAheadBases()`` of 100_000 bases
(``VariantWalkerBase.java:39,146``).  ``VariantLocusWalker.traverse()`` walks
**every reference locus** of a shard that contains a variant
(``VariantLocusWalker.java:150-176`` via ``IntervalLocusIterator``) and calls
``apply(locus, overlappingVariants, ...)`` whenever
``drivingVariants.query(locus)`` is non-empty.  A locus is therefore visited if
and only if some input record's span covers it -- not only where a record
starts.

*Merging* -- in ``ReferenceConfidenceVariantContextMerger.merge()`` every
record that does not start at the locus is a "spanning event"
(``loc.getStart() != vc.getStart()``, ``:150-151``) and its alleles go through
``replaceWithNoCallsAndDels()`` (``:222-245``): the reference becomes NO_CALL,
``<NON_REF>`` stays, and **every ALT shorter than the record's own reference
becomes ``*``**; everything else becomes NO_CALL, which
``filterAllelesForFinalSet()`` (``:297-303``) then drops.  So a spanning record
contributes exactly ``*`` (plus ``<NON_REF>``, removed by
``removeNonRefSymbolicAllele=true`` at ``GenotypeGVCFsEngine.java:136``).

*Dense only* -- with the record starting at the locus and
``includeNonVariants`` set, ``getVariantSubsetToProcess()`` keeps ONLY that
record (``GenotypeGVCFsEngine.java:339-354``); without the flag the merged
``*``-only locus cannot pass ``calculateGenotypes()``
(``GenotypingEngine.java:173-175``) and nothing is written.  Hence the
synthesised rows are dense-only and never appear where a record starts.

*The synthesised row* -- it is an ordinary regenotyped record: QUAL ``0``,
``FILTER=LowQual`` (GenotypeGVCFsEngine.java:416), the full site-annotation set
``AC/AF/AN/DP/ExcessHet/MLEAC/MLEAF/QD`` and
``FORMAT=GT:AD:DP:GQ:PL`` with the spanning record's projected call.  Its
``QD`` renders as the literal ``-0.00`` (negative zero) -- GATK's ``QUAL`` is
``-0.0`` there; recorded, not gated, because native renders ``0.00`` for the
same row shape even when the record is supplied explicitly
(``owned-star-only-locus-dense-negative-zero-qual``, REPORTED ONLY).

*Reverse trim* -- ``GenotypeGVCFsEngine.java:167`` calls
``GATKVariantContextUtils.reverseTrimAlleles()`` (``:1443-1445`` ->
``trimAlleles(vc,false,true)``, ``:1455-1478``) on **every** record that reaches
the polymorphism branch, i.e. in **both** modes, before the site annotations
(``:189``).  The guard at ``:1458`` short-circuits when the record has one
allele or when any allele other than ``*`` is a single base -- which is why
HaplotypeCaller's minimally represented indels are never trimmed.  The number
of clipped bases is ``normalizeAlleles(...).getRight()`` (``:1465-1467``),
reduced by one when the clip would consume an allele completely and no
leading bases are clipped (``:1469-1475``, so an allele can never be emptied --
measured in ``.diag/dense-materialize-probe3.log`` case ``R``:
``ACGTACGT/ACGT`` -> ``ACGTA/A``).  ``*`` and symbolic alleles pass through
untouched (``:1497-1499``) and the stop is recomputed from the trimmed
reference (``:1512``).

Controls pinned here (green before and after the gap is closed)
--------------------------------------------------------------
* ``default-mode-no-synthesis`` -- no dense flag, so the spanning record's locus
  is never published.
* ``starting-record-wins-dense`` -- dense mode with a record starting at the
  covered position: GATK uses only that record, so a fix that materializes
  every covered position must not add a row here.
* ``no-deletion-alt-no-star-rows`` / ``pruned-deletion-no-star-rows`` -- a
  spanning record with no deletion allele, and one whose deletion allele is
  pruned: the materialized rows must not carry ``*`` (GATK publishes them as
  REF-only no-call rows).  These two cases also show the gap is broader than
  ``*``: GATK materializes a row for **every** covered position, whatever the
  spanning allele is, while native materializes none of them; the control
  therefore asserts the ``*``-free property only, and records the full GATK
  rows as the measured truth.
* ``no-trim-needed-multi-base-alt`` -- a multi-base substitution with no shared
  trailing base must be emitted untrimmed.

Exit status: 0 when every gated case passes, non-zero otherwise.
``--expect-divergence`` turns the run into a diagnostic that always exits 0.
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import oracle_guard  # noqa: E402

# Copied verbatim from verify_genotype_gvcf.py:328-337 (via
# verify_genotype_gvcf_spandel_gatk_oracle.py) so this gate pins the same
# fixture bytes as the registered contract test.
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

DENSE = "--include-non-variant-sites"

# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

# One emitted deletion at 2 whose span is 2-4 (REF AAA), so GATK's traversal
# also visits 3 and 4 -- positions no record starts at.
DELETION_ALONE = (
    "chr1\t2\t.\tAAA\tA,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\n"
)

# The same deletion plus a downstream orphan-'*' locus with a concrete ALT.
DELETION_PLUS_DOWNSTREAM_STAR = DELETION_ALONE + (
    "chr1\t5\t.\tA\t*,G,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/2:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100\n"
)

# The shape the previous round left ungated: a '*' record with a two-base REF
# at a locus an emitted deletion also spans.
MERGED_STAR_RECORD = DELETION_ALONE + (
    "chr1\t4\t.\tAA\t*,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t1/1:20:0,20:100,100,0\n"
    "chr1\t5\t.\tA\t*,G,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/2:20:0,0,20,0:100,100,100,100,0,100,100,100,100,100\n"
)

# Reverse trim with no '*' anywhere: every non-symbolic allele is two or more
# bases and all of them end in 'A', so one trailing base is clipped.
SUFFIX_SUBSTITUTION = (
    "chr1\t2\t.\tAAAA\tAACA,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\n"
)

# Control: a spanning record with NO deletion allele (ALT as long as REF), so no
# '*' can be synthesised -- but GATK still materializes the covered positions
# (3, 4 and 5) as REF-only no-call rows.  The ALT is chosen so the two alleles
# share no trailing base, keeping the reverse trim out of this control.
SPANNING_SUBSTITUTION = (
    "chr1\t2\t.\tAAAA\tAACC,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\n"
)

# Control: the same deletion, but its likelihoods are homozygous-reference, so
# the deletion allele is pruned at locus 2; the covered loci are still visited
# and published as REF-only rows without any '*'.
PRUNED_DELETION = (
    "chr1\t2\t.\tAAA\tA,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/0:20:20,0,0:0,100,100,100,100,100\n"
)

# Control: dense mode where a record DOES start at the only covered position
# (2 AA A spans 2-3), so nothing may be synthesised.
STARTING_RECORD_WINS = (
    "chr1\t2\t.\tAA\tA,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\n"
    "chr1\t3\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\n"
)

# Control: the same 'starting record wins' fixture with the starting record
# carrying the symbolic '*'; the row renders, but GATK's QD is the literal
# '-0.00' while native writes '0.00' (REPORTED ONLY).
OWNED_STAR_ONLY_LOCUS = (
    "chr1\t2\t.\tAA\tA,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\n"
    "chr1\t3\t.\tA\t*,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\n"
)

# REPORTED ONLY (same fixture, default mode): GATK refuses the locus outright
# when the surviving ALT set is exactly [SPAN_DEL] and the traversal is not
# EMIT_ALL_ACTIVE_SITES (GenotypingEngine.java:173-175), while native publishes
# the '*'-only row.  A different call site from this round's two behaviours.
OWNED_STAR_ONLY_LOCUS_DEFAULT = OWNED_STAR_ONLY_LOCUS

# Control: a multi-base ALT that shares no trailing base with its REF -- the
# guard at GATKVariantContextUtils.java:1495 must return the record unchanged.
NO_TRIM_NEEDED = (
    "chr1\t2\t.\tAAAA\tAAC,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\n"
)

# ---------------------------------------------------------------------------
# measured rows (GRCh37-free 100 bp poly-A chr1 reference, see write_reference)
# ---------------------------------------------------------------------------

ROW_DELETION_SITE = ("chr1\t2\t.\tAAA\tA\t92.60\t.\t"
                     "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;"
                     "MLEAC=1;MLEAF=0.500;QD=4.63\t"
                     "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")
# The synthesised row: QUAL 0 -> FILTER LowQual, the spanning record's projected
# call, and QD as the literal negative zero.
ROW_SYNTHESISED_STAR = ("chr1\t3\t.\tA\t*\t0\tLowQual\t"
                        "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;"
                        "MLEAC=1;MLEAF=0.500;QD=-0.00\t"
                        "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")
ROW_SYNTHESISED_STAR_AT_4 = ROW_SYNTHESISED_STAR.replace("\t3\t", "\t4\t", 1)
# Locus 4 of the previous round's shape: its own '*' record (1/1) wins over the
# spanning deletion, and the merged REF 'AA' is reverse-trimmed to 'A'.
ROW_TRIMMED_STAR_RECORD = ("chr1\t4\t.\tA\t*\t0\tLowQual\t"
                           "AC=2;AF=1.00;AN=2;DP=20;ExcessHet=0.0000;"
                           "MLEAC=2;MLEAF=1.00;QD=0.00\t"
                           "GT:AD:DP:GQ:PL\t1/1:0,20:20:99:100,100,0")
ROW_STAR_PLUS_CONCRETE_AT_5 = (
    "chr1\t5\t.\tA\t*,G\t82.26\t.\t"
    "AC=1,1;AF=0.500,0.500;AN=2;DP=20;ExcessHet=0.0000;"
    "MLEAC=1,1;MLEAF=0.500,0.500;QD=4.11\t"
    "GT:AD:DP:GQ:PL\t1/2:0,0,20:20:99:100,100,100,100,0,100")
ROW_CONCRETE_G_ONLY_AT_5 = ("chr1\t5\t.\tA\tG\t82.26\t.\t"
                            "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;"
                            "MLEAC=1;MLEAF=0.500;QD=4.11\t"
                            "GT:AD:DP:PL\t0/1:0,20:20:0,0,0")
ROW_TRIMMED_SUFFIX_SITE = ("chr1\t2\t.\tAAA\tAAC\t92.64\t.\t"
                           "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;"
                           "MLEAC=1;MLEAF=0.500;QD=4.63\t"
                           "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")
ROW_SPANNING_SUBSTITUTION_SITE = (
    "chr1\t2\t.\tAAAA\tAACC\t92.64\t.\t"
    "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;"
    "MLEAC=1;MLEAF=0.500;QD=4.63\t"
    "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")


def _ref_only_row(position: int, qual: str) -> str:
    """The dense REF-only materialization of a covered locus, as measured.

    Two shapes appear in this gate: the monomorphic-recovery shape
    (``QUAL 136.80/192.21``, ``FORMAT=GT``, ``GT ./.``) for a locus whose only
    surviving allele was a pruned deletion, and the merge-null shape
    (``QUAL .``, ``FORMAT=GT:AD``, ``GT:AD ./. :0``) for a locus all of whose
    spanning records replaced every ALT with NO_CALL.
    """
    assert qual in (".", "136.80", "192.21")
    if qual == ".":
        return (f"chr1\t{position}\t.\tA\t.\t.\t.\tDP=20\t"
                "GT:AD\t./.:0")
    return (f"chr1\t{position}\t.\tA\t.\t{qual}\t.\t"
            "DP=20;MLEAC=.;MLEAF=.\tGT\t./.")


ROW_PRUNED_DELETION_SITE = ("chr1\t2\t.\tAAA\t.\t136.80\t.\t"
                            "DP=20;MLEAC=.;MLEAF=.\tGT\t./.")
ROW_STARTING_RECORD_DELETION = ("chr1\t2\t.\tAA\tA\t92.60\t.\t"
                                "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;"
                                "MLEAC=1;MLEAF=0.500;QD=4.63\t"
                                "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")
ROW_STARTING_RECORD_SNP = ("chr1\t3\t.\tA\tG\t92.64\t.\t"
                           "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;"
                           "MLEAC=1;MLEAF=0.500;QD=4.63\t"
                           "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")
ROW_OWNED_STAR_ONLY_LOCUS = ("chr1\t3\t.\tA\t*\t0\tLowQual\t"
                             "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;"
                             "MLEAC=1;MLEAF=0.500;QD=-0.00\t"
                             "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")
ROW_NO_TRIM_NEEDED = ("chr1\t2\t.\tAAAA\tAAC\t92.60\t.\t"
                      "AC=1;AF=0.500;AN=2;DP=20;ExcessHet=0.0000;"
                      "MLEAC=1;MLEAF=0.500;QD=4.63\t"
                      "GT:AD:DP:GQ:PL\t0/1:0,20:20:99:100,0,100")

# ---------------------------------------------------------------------------
# cases
# ---------------------------------------------------------------------------

# ``mode`` is the assertion applied to a gated case:
#   "rows"    -- GATK and native rows must be equal, count included;
#   "no-star" -- GATK rows must equal the measured truth and no native row may
#                carry the symbolic '*' ALT; every native row must still be
#                byte-identical to the GATK row at its own position (a control
#                against a fix that over-materializes).
CASES = [
    {
        "case": "deletion-alone-dense-materializes-covered-positions",
        "why": "one emitted deletion (REF AAA at 2, span 2-4) and nothing else. "
               "GATK's locus traversal visits 3 and 4 because the record spans "
               "them (VariantLocusWalker.java:150-176), the merger replaces the "
               "deletion allele of that spanning event with '*' "
               "(ReferenceConfidenceVariantContextMerger.java:150-151,222-245) "
               "and dense mode publishes the locus "
               "(GenotypeGVCFsEngine.java:339-354,371-384).  Native's locus list "
               "is built from record starts only, so both materialized rows are "
               "missing",
        "body": DELETION_ALONE,
        "args": [DENSE],
        "mode": "rows",
        "gated": True,
        "expect": [ROW_DELETION_SITE, ROW_SYNTHESISED_STAR,
                   ROW_SYNTHESISED_STAR_AT_4],
    },
    {
        "case": "deletion-plus-downstream-star-dense",
        "why": "the same deletion plus a downstream '*'-with-concrete-ALT "
               "record: the two synthesized rows stay, and the locus at 5 is "
               "unchanged, so the row count is 4 against native's 2",
        "body": DELETION_PLUS_DOWNSTREAM_STAR,
        "args": [DENSE],
        "mode": "rows",
        "gated": True,
        "expect": [ROW_DELETION_SITE, ROW_SYNTHESISED_STAR,
                   ROW_SYNTHESISED_STAR_AT_4, ROW_CONCRETE_G_ONLY_AT_5],
    },
    {
        "case": "merged-star-record-dense-reverse-trim",
        "why": "the shape the previous round left ungated: locus 4 carries its "
               "own two-base-REF '*' record, so dense mode keeps only that "
               "record (GenotypeGVCFsEngine.java:339-354) and the merged "
               "'AA/*' record is reverse-trimmed to 'A/*' "
               "(GenotypeGVCFsEngine.java:167; GATKVariantContextUtils.java:"
               "1443-1478) -- while locus 3, covered only by the deletion, is "
               "materialized.  Native emits the untrimmed REF and no row at 3",
        "body": MERGED_STAR_RECORD,
        "args": [DENSE],
        "mode": "rows",
        "gated": True,
        "expect": [ROW_DELETION_SITE, ROW_SYNTHESISED_STAR,
                   ROW_TRIMMED_STAR_RECORD, ROW_STAR_PLUS_CONCRETE_AT_5],
    },
    {
        "case": "suffix-substitution-default-reverse-trim",
        "why": "the reverse trim is NOT dense-only: this is the default "
               "traversal, with no '*' anywhere, and GATK still clips the "
               "trailing base the two non-symbolic alleles share "
               "(GenotypeGVCFsEngine.java:167, reached from :160-169).  Native "
               "emits 'AAAA/AACA' untrimmed",
        "body": SUFFIX_SUBSTITUTION,
        "args": [],
        "mode": "rows",
        "gated": True,
        "expect": [ROW_TRIMMED_SUFFIX_SITE],
    },
    {
        "case": "no-deletion-alt-no-star-rows",
        "why": "control: a spanning record whose ALT is as long as its REF "
               "contributes no allele at all (replaceWithNoCallsAndDels maps it "
               "to NO_CALL, :222-244), so the materialized rows at 3, 4 and 5 "
               "must not carry '*'.  It also shows the gap is wider than '*': "
               "GATK materializes those positions (as REF-only no-call rows) "
               "and native emits nothing",
        "body": SPANNING_SUBSTITUTION,
        "args": [DENSE],
        "mode": "no-star",
        "gated": True,
        "expect": [ROW_SPANNING_SUBSTITUTION_SITE,
                   _ref_only_row(3, "."), _ref_only_row(4, "."),
                   _ref_only_row(5, ".")],
    },
    {
        "case": "pruned-deletion-no-star-rows",
        "why": "control: the deletion allele is pruned at locus 2 (hom-ref "
               "likelihoods), so the covered loci 3 and 4 are published as "
               "monomorphic REF-only rows and must not carry a synthesized '*'",
        "body": PRUNED_DELETION,
        "args": [DENSE],
        "mode": "no-star",
        "gated": True,
        "expect": [ROW_PRUNED_DELETION_SITE, _ref_only_row(3, "192.21"),
                   _ref_only_row(4, "192.21")],
    },
    {
        "case": "starting-record-wins-dense",
        "why": "control: dense mode with a record starting at every covered "
               "position (2 AA A spans 2-3 and 3 starts there).  GATK keeps "
               "only the starting record (:339-354), so a fix must not add a "
               "synthesized row or let the spanning deletion contribute '*'",
        "body": STARTING_RECORD_WINS,
        "args": [DENSE],
        "mode": "rows",
        "gated": True,
        "expect": [ROW_STARTING_RECORD_DELETION, ROW_STARTING_RECORD_SNP],
    },
    {
        "case": "default-mode-no-synthesis",
        "why": "control: the same deletion without --include-non-variant-sites "
               "produces exactly one row.  Without the flag the spanning locus "
               "never reaches the writer (GenotypingEngine.java:173-175) and "
               "the row materialization must be dense-only",
        "body": DELETION_ALONE,
        "args": [],
        "mode": "rows",
        "gated": True,
        "expect": [ROW_DELETION_SITE],
    },
    {
        "case": "no-trim-needed-multi-base-alt",
        "why": "control: a multi-base substitution whose alleles share no "
               "trailing base -- the guard at GATKVariantContextUtils.java:1495 "
               "returns the record unchanged, so native's untrimmed row is "
               "already correct and a fix must leave it alone",
        "body": NO_TRIM_NEEDED,
        "args": [],
        "mode": "rows",
        "gated": True,
        "expect": [ROW_NO_TRIM_NEEDED],
    },
    {
        "case": "owned-star-only-locus-dense-negative-zero-qual",
        "why": "REPORTED ONLY, different root cause: the same row shape with "
               "the record supplied explicitly in the input diverges by the "
               "rendering of QUAL's negative zero.  GATK's QUAL is -0.0 and "
               "QD prints as '-0.00'; native prints '0.00'.  Not gated here "
               "because it is a number-formatting defect reachable without the "
               "materialization gap, not a record-supply one",
        "body": OWNED_STAR_ONLY_LOCUS,
        "args": [DENSE],
        "mode": "rows",
        "gated": False,
        "expect": [ROW_STARTING_RECORD_DELETION, ROW_OWNED_STAR_ONLY_LOCUS],
    },
    {
        "case": "owned-star-only-locus-default-mode-refused",
        "why": "GATED since the round that implemented it: without the dense flag "
               "GATK returns no record at all for a locus whose only surviving "
               "ALT is the symbolic '*' (GenotypingEngine.java:173-175: "
               "'!emitAllActiveSites() && alleles.size()==1 && "
               "SPAN_DEL.equals(alleles.get(0))'), and native used to publish the "
               "'*'-only row with the full annotation set.  Found by this gate's "
               "round and fixed later; the registered spandel oracle's "
               "covered-star-only-record case reaches native's *empty*-ALT-set "
               "refusal instead, so this shape was previously unmeasured",
        "body": OWNED_STAR_ONLY_LOCUS_DEFAULT,
        "args": [],
        "mode": "rows",
        "gated": True,
        "expect": [ROW_STARTING_RECORD_DELETION],
    },
]


def data_rows(path: pathlib.Path) -> list[str]:
    rows: list[str] = []
    if not path.exists():
        return rows
    with open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line and not line.startswith("#"):
                rows.append(line.rstrip("\r\n"))
    return rows


def write_reference(work: pathlib.Path) -> pathlib.Path:
    """A 100 bp chr1 whose bases match every fixture REF allele."""
    reference = work / "reference.fa"
    reference.write_text(">chr1\n" + "A" * 100 + "\n", encoding="utf-8")
    reference.with_name(reference.name + ".fai").write_text(
        "chr1\t100\t6\t100\t101\n", encoding="utf-8")
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


def alt_column(row: str) -> str:
    fields = row.split("\t")
    return fields[4] if len(fields) > 4 else ""


def position_of(row: str) -> str:
    fields = row.split("\t")
    return f"{fields[0]}:{fields[1]}" if len(fields) > 1 else row


def run_case(case: dict, work: pathlib.Path, reference: pathlib.Path,
             java: pathlib.Path, jar: pathlib.Path, native: pathlib.Path,
             timeout: int) -> dict:
    source = work / f"{case['case']}.in.vcf"
    source.write_text(HEADER + case["body"], encoding="utf-8")
    # GATK's group-by-locus traversal needs an index; creating it for every case
    # keeps the two tools' inputs byte-identical.
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
    native_result = invoke([str(native), *common,
                            "--gatk-compatible-annotations",
                            *case.get("native_args", []),
                            "-O", str(native_out)],
                           f"native GenotypeGVCFs [{case['case']}]", timeout)

    gatk_rows = data_rows(gatk_out)
    native_rows = data_rows(native_out)

    result = {
        "case": case["case"],
        "why": case["why"],
        "gated": case["gated"],
        "mode": case["mode"],
        "args": common,
        "native_args": case.get("native_args", []),
        "gatk_exit": gatk_result.returncode,
        "native_exit": native_result.returncode,
        "expect": case["expect"],
        "gatk_rows": gatk_rows,
        "native_rows": native_rows,
        "violations": [],
    }
    if index_result.returncode != 0:
        result["violations"].append(
            f"GATK IndexFeatureFile exited {index_result.returncode}")
    if gatk_result.returncode != 0:
        result["violations"].append(f"GATK exited {gatk_result.returncode}")
    if native_result.returncode != 0:
        result["violations"].append(f"native exited {native_result.returncode}")

    # Always: the pinned GATK must still emit exactly what this gate was
    # written from, or the comparison below could pass for the wrong reason.
    if gatk_rows != case["expect"]:
        result["violations"].append(
            "GATK moved away from the measured truth: "
            f"GATK={gatk_rows} expected={case['expect']}")

    if case["gated"]:
        if case["mode"] == "rows":
            if len(gatk_rows) != len(native_rows):
                result["violations"].append(
                    f"record count differs: GATK={len(gatk_rows)} "
                    f"native={len(native_rows)}; missing native positions="
                    f"{sorted(set(map(position_of, gatk_rows)) - set(map(position_of, native_rows)))}")
            for index, (gatk_row, native_row) in enumerate(zip(gatk_rows, native_rows)):
                if gatk_row != native_row:
                    result["violations"].append(
                        f"row {index} is not byte-identical: GATK={gatk_row!r} "
                        f"NATIVE={native_row!r}")
        elif case["mode"] == "no-star":
            starred = [row for row in native_rows if alt_column(row).find("*") >= 0]
            if starred:
                result["violations"].append(
                    f"native synthesised a '*' row where no deletion allele "
                    f"reaches the locus: {starred}")
            expected_by_position = {position_of(row): row for row in gatk_rows}
            for row in native_rows:
                reference_row = expected_by_position.get(position_of(row))
                if reference_row is None:
                    result["violations"].append(
                        f"native emitted a row at a position GATK does not: {row!r}")
                elif reference_row != row:
                    result["violations"].append(
                        f"native row at {position_of(row)} is not byte-identical: "
                        f"GATK={reference_row!r} NATIVE={row!r}")
        else:
            raise AssertionError(f"unknown mode {case['mode']!r}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Strict oracle for GenotypeGVCFs dense-mode record "
                    "materialization (spanning-deletion row synthesis and "
                    "reverse trimming of merged records) against pinned GATK "
                    "4.6.2.0.")
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
            "verify_genotype_gvcf_dense_materialize_gatk_oracle.py", java, jar)
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
            prefix="fastgatk-genotype-dense-materialize-oracle-") as directory:
        work = pathlib.Path(directory)
        reference = write_reference(work)
        for case in selected:
            results.append(run_case(case, work, reference, java, jar, native,
                                    arguments.timeout))

    violations: list[str] = []
    for result in results:
        if result["gated"]:
            violations.extend(f"[{result['case']}] {item}"
                              for item in result["violations"])

    print(f"# {pathlib.Path(__file__).name}: pinned GATK 4.6.2.0 vs native ({native})")
    print("# GATK rules under test: (1) every reference locus covered by an input "
          "record is visited (VariantLocusWalker.java:150-176), a spanning "
          "record's shorter ALT becomes '*' (ReferenceConfidenceVariantContext"
          "Merger.java:150-151,222-245) and dense mode publishes the locus "
          "(GenotypeGVCFsEngine.java:339-354,371-384); (2) every regenotyped "
          "record is reverse-trimmed (GenotypeGVCFsEngine.java:167; "
          "GATKVariantContextUtils.java:1443-1478).")
    print(f"# mode: {'expect-divergence (diagnostic)' if not strict_mode else 'strict assertion'}"
          f"; cases={[case['case'] for case in selected]}")
    for result in results:
        marker = f"gated/{result['mode']}" if result["gated"] else "REPORTED ONLY (not gated)"
        print(f"[{result['case']}] {marker}")
        print(f"    why: {result['why']}")
        print(f"    args: {' '.join(result['args']) if result['args'] else '(none)'}")
        print(f"    gatk_exit={result['gatk_exit']} native_exit={result['native_exit']}")
        print(f"    expected rows ({len(result['expect'])}):")
        for row in result["expect"]:
            print(f"        {row}")
        print(f"    GATK   rows ({len(result['gatk_rows'])}):")
        for row in result["gatk_rows"]:
            print(f"        {row}")
        print(f"    NATIVE rows ({len(result['native_rows'])}):")
        for row in result["native_rows"]:
            print(f"        {row}")
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
        "control_modes": {
            "rows": "full byte-identical row comparison, count included",
            "no-star": "GATK truth pinned; no native row may carry '*'; every "
                       "native row must equal the GATK row at its position",
        },
    }
    print(json.dumps(payload, sort_keys=True))

    if not strict_mode:
        return 0
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
