#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace fastgatk::kernels {

// A half-open coordinate range traversed by the source AssemblyRegionWalker.
// GATK's IntervalAlignmentContextIterator emits an empty pileup for every
// uncovered coordinate inside this range.  The native profile keeps that
// traversal information separate from the compact Kokkos evidence arrays so
// a large uncovered interval need not become a device allocation.
struct ActivityProfileSpan {
    std::int32_t tid = -1;
    std::int32_t start = 0;
    std::int32_t end = 0;
};

struct ActivityProfileInput {
    std::vector<std::int32_t> tids;
    std::vector<std::int32_t> positions;
    // Four deterministic base counts per locus in A,C,G,T order.
    std::vector<std::uint32_t> counts;
    // Optional reference base per locus, encoded as A,C,G,T -> 0,1,2,3 and
    // any other value for an unknown/N reference.  When present, activity is
    // the non-reference fraction (rather than merely 1-max(count)/depth),
    // matching the Host reference projection used by HC's active-region
    // model.  An empty vector retains the legacy no-reference fallback.
    std::vector<std::uint8_t> reference_bases;
    // Number of explicit CIGAR insertion/deletion observations anchored at
    // each locus.  GATK's ReadThreading/AssemblyRegion path seeds activity
    // from these events even when the aligned bases themselves are all
    // reference-matching; keeping this parallel to `tids` lets the portable
    // Kokkos signal consume the same evidence for HC and Mutect2.
    std::vector<std::uint32_t> indel_counts;
    // Optional per-locus pileup observations for the GATK activity models.
    // `quality_offsets` has loci+1 entries and indexes the parallel
    // `quality_bases`/`quality_values`/`quality_alt_flags` arrays.  Mutect2
    // uses the bases and qualities for its somatic log-likelihood ratio;
    // a set flag means that the aligned base was classified as indel evidence
    // (a current/following indel or useful adjacent soft clip) and therefore
    // must not also enter the substitution bucket.  HaplotypeCaller consumes
    // the same flags as CIGAR-adjacent ALT evidence in Ref-vs-Any.
    std::vector<std::uint32_t> quality_offsets;
    std::vector<std::uint8_t> quality_bases;
    std::vector<std::uint8_t> quality_values;
    std::vector<std::uint8_t> quality_alt_flags;
    // Mutect2 PileupQualBuffer keeps insertion/deletion (including useful
    // terminal soft-clip) evidence in one separate ALT bucket.  These
    // per-locus quality lists retain its LOD calculation on Kokkos instead
    // of treating every CIGAR indel as unconditionally active on Host.
    std::vector<std::uint32_t> indel_quality_offsets;
    std::vector<std::uint8_t> indel_quality_values;
    // Optional Host-derived Mutect2 normal-pileup gate, one byte per locus.
    // `isActive()` first evaluates the tumor LOD, then deactivates a locus
    // when the normal supports the same dominant allele above its source
    // thresholds.  Read/sample ownership stays on Host; the Kokkos kernel
    // consumes only this compact, deterministic suppression mask.
    std::vector<std::uint8_t> somatic_normal_suppressed;
    // Optional Host-derived Mutect2 germline-resource/PON gate, one byte per
    // locus. Resource VCF parsing and dominant-allele selection remain Host
    // work; Kokkos receives the resulting fixed-shape activity veto only.
    std::vector<std::uint8_t> somatic_feature_suppressed;
    // Optional Host-derived HaplotypeCaller GenotypeGivenAlleles trigger,
    // one byte per locus.  A concrete, eligible --alleles VCF event makes
    // HaplotypeCallerEngine.isActive() return probability one before it
    // evaluates pileup evidence.  VCF parsing/filter policy stays on Host;
    // this Kokkos mask preserves the resulting per-locus activity state
    // without broadening a feature trigger into --force-active.
    std::vector<std::uint8_t> forced_allele_active;
    // Per-locus running-average inputs for HaplotypeCaller's
    // HIGH_QUALITY_SOFT_CLIPS ActivityProfileState.  A nonzero average over
    // the matching event count causes the Host profile state machine to copy
    // the Ref-vs-Any probability to neighbouring loci before band-pass
    // filtering, exactly as GATK's ActivityProfile.processState does.
    std::vector<std::uint32_t> high_quality_softclip_bases;
    std::vector<std::uint32_t> high_quality_softclip_events;
};

struct ActivityProfileOptions {
    std::uint32_t min_depth = 4;
    // GATK's BandPassActivityProfile default active probability threshold.
    double min_activity = 0.002;
    std::uint32_t halo = 128;
    std::uint32_t max_region_size = 300;
    // Maximum distance over which active probability is propagated into the
    // same AssemblyRegion.  This is intentionally separate from the hard
    // maximum region length: GATK uses both controls, and conflating them can
    // create unbounded graph/PairHMM windows on dense active runs.
    std::uint32_t max_prob_propagation_distance = 300;
    // BandPassActivityProfile controls.  The Gaussian filter is deliberately
    // kept in the Host state machine because it creates a variable-length,
    // accumulated profile; the raw per-locus activity remains a Kokkos
    // kernel.  Defaults match GATK's AssemblyRegionArgumentCollection and
    // BandPassActivityProfile constants.
    std::uint32_t min_region_size = 50;
    std::uint32_t bandpass_max_filter_size = 50;
    double bandpass_sigma = 17.0;
    bool bandpass_adaptive_filter = true;
    // When enabled and quality observations are present, use Mutect2's
    // pileup log-likelihood-ratio activity signal instead of a raw mismatch
    // fraction.  The threshold is in log10 odds, matching the CLI spelling;
    // the Kokkos kernel converts it to natural-log space for the GATK formula.
    bool quality_aware_somatic = false;
    // HaplotypeCaller isActive() does not threshold raw mismatch fraction.
    // It evaluates a diploid Ref-vs-Any pileup likelihood, then combines it
    // with the same SNP heterozygosity/standard-deviation pseudocounts used
    // by GATK's minimal active-region genotyper.  Keep this distinct from
    // Mutect2's quality-aware somatic LOD model so the shared kernel remains
    // explicit about which Java engine it is matching.
    bool genotype_aware_hc = false;
    std::uint8_t hc_min_base_quality = 10;
    double hc_snp_heterozygosity = 1.0e-3;
    double hc_heterozygosity_stdev = 0.01;
    // GATK Mutect2 --force-active marks every emitted profile segment as an
    // active AssemblyRegion while preserving the ActivityProfile boundaries.
    // Keep the switch at the profile/region boundary so callers can force
    // graph/PairHMM materialization without changing the Kokkos activity
    // signal itself.
    bool force_active = false;
    double initial_tumor_log10_odds = 2.0;
    std::uint8_t pcr_snv_quality = 40;
    std::uint8_t multiple_substitution_quality_correction = 5;
    // Optional source traversal ranges.  Each range is half-open and is
    // interpreted as zero activity at coordinates absent from
    // ActivityProfileInput.  Host-side sparse traversal preserves GATK's
    // band-pass boundary semantics without staging empty pileups on Kokkos.
    std::vector<ActivityProfileSpan> traversal_spans;
};

struct ActivityRegion {
    std::int32_t tid = -1;
    std::int32_t start = -1;
    std::int32_t end = -1;
    std::int32_t active_start = -1;
    std::int32_t active_end = -1;
    std::size_t active_loci = 0;
    // True for an active profile segment and false for the intervening
    // inactive segment.  The calling path still consumes `regions` (active
    // padded windows); `profile_regions` retains both states for the GATK
    // AssemblyRegion IGV/debug contract.
    bool active = true;
};

struct ActivityProfileResult {
    bool used = false;
    bool reference_aware = false;
    std::size_t loci = 0;
    std::size_t active_loci = 0;
    std::vector<std::uint8_t> active;
    std::vector<double> activity;
    std::vector<ActivityRegion> regions;
    // Raw, unpadded profile segments in traversal order.  This is deliberately
    // separate from `regions`, whose active entries are halo-expanded for
    // graph/PairHMM ownership.
    std::vector<ActivityRegion> profile_regions;
    std::uint32_t filter_size = 0;
    std::uint32_t effective_max_prob_propagation_distance = 0;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

// Per-locus output of Mutect2's SomaticReferenceConfidenceModel.  This is
// intentionally separate from ActivityProfileResult: activity detection uses
// PileupQualBuffer's *best* alternative allele (and its PCR/error-quality
// policy), whereas reference confidence treats every post-finalization
// non-reference pileup element as evidence for <NON_REF>.  The Host owns
// CIGAR/finalizeRegion projection and supplies the compact primitive arrays;
// Kokkos owns the quality/LOD reduction so reference-confidence output keeps
// the same portable numerical boundary as the rest of the caller.
struct SomaticReferenceConfidenceResult {
    std::vector<double> log10_lods;
    std::vector<std::uint32_t> depths;
    std::vector<std::uint32_t> reference_depths;
    std::vector<std::uint32_t> non_reference_depths;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

// Compute a deterministic pileup activity signal on Kokkos and construct
// sorted, halo-expanded AssemblyRegions on Host.  The region builder is kept
// outside the kernel because interval merging and stable ordering are
// variable-length control flow, while the per-locus signal is portable.
ActivityProfileResult compute_activity_profile_kokkos(
    const ActivityProfileInput& input, ActivityProfileOptions options = {});

// Port SomaticReferenceConfidenceModel::calcGenotypeLikelihoodsOfRefVsAny
// for the default flat Beta(1,1) allele-fraction prior.  `minimum_base_quality`
// is the source model's strict lower bound (elements with quality <= it are
// excluded).  The result is in the log10 space used by the TLOD FORMAT field
// and SomaticGVCFWriter's LOD partitions.
SomaticReferenceConfidenceResult calculate_somatic_reference_confidence_kokkos(
    const ActivityProfileInput& input, std::uint8_t minimum_base_quality = 6);

}  // namespace fastgatk::kernels
