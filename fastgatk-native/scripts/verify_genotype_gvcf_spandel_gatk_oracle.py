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
# SURVIVES the subset -- but the record has no other ALT, so GATK drops it at
# GenotypeGVCFs.apply() (isSpanningDeletionOnly).  Reported only: native writes
# the '*' record, which is a separate record-emission rule, not ownership.
STAR_ONLY_COVERED_RECORD = (
    "chr1\t2\t.\tAA\tA,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,0,100,100,100,100\n"
    "chr1\t3\t.\tA\t*,<NON_REF>\t.\tPASS\tDP=20\t"
    "GT:DP:AD:PL\t0/1:20:0,20,0:100,100,100,0,100,100\n"
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
        "why": "reported only: this '*' IS covered by the upstream deletion, so "
               "ownership is right, but the record has no other ALT and GATK "
               "drops it at GenotypeGVCFs.apply() because "
               "GATKVariantContextUtils.isSpanningDeletionOnly() is true "
               "(GenotypeGVCFs.java:327-328).  Native writes the '*' record. "
               "That is a separate record-emission rule, deliberately not "
               "gated in this round",
        "body": STAR_ONLY_COVERED_RECORD,
        "args": [],
        "gated": False,
        "expect": [GATK_DEL_UPSTREAM_LOCUS_ROW],
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
    source.write_text(HEADER + case["body"], encoding="utf-8")
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
                            "-O", str(native_out)],
                           f"native GenotypeGVCFs [{case['case']}]", timeout)

    gatk_rows = data_rows(gatk_out) if gatk_out.exists() else []
    native_rows = data_rows(native_out) if native_out.exists() else []

    result = {
        "case": case["case"],
        "why": case["why"],
        "gated": case["gated"],
        "args": common,
        "gatk_exit": gatk_result.returncode,
        "native_exit": native_result.returncode,
        "expect": case["expect"],
        "gatk_rows": gatk_rows,
        "native_rows": native_rows,
        "violations": [],
    }
    if index_result.returncode != 0:
        result["violations"].append(f"GATK IndexFeatureFile exited {index_result.returncode}")
    if gatk_result.returncode != 0:
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
        print(f"    gatk_exit={result['gatk_exit']} native_exit={result['native_exit']}")
        if result["expect"] is not None:
            print(f"    expected rows: {result['expect']}")
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
    }
    print(json.dumps(payload, sort_keys=True))

    if not strict_mode:
        return 0
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
