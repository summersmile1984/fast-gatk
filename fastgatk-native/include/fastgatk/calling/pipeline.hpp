#pragma once

#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/kernels/activity_profile.hpp"
#include "fastgatk/kernels/kmer_graph.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace fastgatk::calling {

// Result of the flow-space long-homopolymer restoration pass. GATK builds
// flow likelihoods from a collapsed haplotype view, then restores deleted
// reference HMER bases before the haplotypes re-enter genotyping. The
// alignment score/CIGAR recurrence is shared with the Kokkos SW kernel; this
// small Host result keeps the variable-length sequence transformation
// testable without leaking HTSlib state into a device kernel.
struct FlowHmerUncollapseResult {
    std::vector<std::uint8_t> bases;
    std::int32_t alignment_offset = 0;
    bool expanded = false;
};

FlowHmerUncollapseResult uncollapse_flow_hmers(
    const std::vector<std::uint8_t>& bases,
    const std::vector<std::uint8_t>& reference,
    std::size_t hmer_size_threshold,
    bool partial_mode = false,
    bool limit_to_hmer_size_threshold = false);

// A caller may score this concrete allele even when its own assembly did not
// discover the branch. Mutect2 uses this for matched normals: tumor assembly
// supplies the allele/haplotype set and both samples receive likelihoods on
// that shared set.
struct ForcedAllele {
    std::int32_t tid = -1;
    std::int32_t position = -1;  // 0-based VCF anchor
    std::string reference;
    std::string alternate;
    // A recovered low-quality tumor SNP can require a nearby insertion on
    // its assembled alternate haplotype. Preserve that path context when the
    // allele is replayed against matched-normal reads.
    bool somatic_low_quality_assembly_recovery = false;
    std::int32_t linked_insertion_reference_position = -1;
    std::string linked_insertion_bases;
    // Preserve the source AssemblyResultSet's branch support when an allele
    // is replayed against another sample. This is graph topology evidence,
    // not that sample's AD: ChainPruner must not reinterpret a shared ALT as
    // a zero-support normal-only branch before PairHMM has scored its reads.
    std::uint32_t reference_count = 0;
    std::uint32_t alternate_count = 0;
    // Retained EventMap provenance.  A shared AssemblyResultSet may also
    // carry an explicitly read-supported compatibility event while compact
    // graph projection is being completed, but only a graph-derived allele
    // is a direct haplotype-to-reference EventMap projection.
    bool graph_derived = false;
};

// GATK assigns an EventMap to every retained Haplotype, not only to the
// eventual union of candidate alleles.  Keep that irregular, per-haplotype
// state on the C++ Host: it has variable-length alleles and CIGAR-dependent
// compound-event rules, so it is deliberately outside the Kokkos ABI.
struct HaplotypeEventMap {
    std::vector<ForcedAllele> events;
};

// The sequence interval retained by GATK Haplotype.trim() for a reference
// half-open span. The offsets index the original haplotype sequence. A
// missing value means the target span begins or ends in a deletion (or the
// CIGAR cannot describe exactly the supplied haplotype/reference lengths),
// which GATK discards for a non-reference haplotype. In particular, an
// insertion immediately before the left boundary and immediately after the
// right boundary is not retained.
struct HaplotypeCigarTrim {
    std::size_t sequence_begin = 0;
    std::size_t sequence_end = 0;
};

// Host-side coordinate component of Haplotype.trim(). The graph, SW score,
// and PairHMM arithmetic remain Kokkos work; this function only preserves the
// variable-length CIGAR boundary semantics before a trimmed haplotype enters
// PairHMM.
std::optional<HaplotypeCigarTrim> trim_haplotype_by_reference_cigar(
    const std::string& cigar, std::size_t reference_length,
    std::size_t haplotype_length, std::size_t trim_start,
    std::size_t trim_end);

// The identity and reference flag of one already-trimmed haplotype. This is
// deliberately Host state: AssemblyResultSet.trimTo() removes equal-base
// haplotypes, lets the reference source win that collision, and then orders
// the retained set before it reaches the likelihood engine.
struct TrimmedHaplotypePopulationEntry {
    std::string bases;
    bool is_reference = false;
};

// Returns indices into `haplotypes` for the exact AssemblyResultSet.trimTo()
// population order: one representative per equal base string (the reference
// representative wins), sorted by length and then lexicographic bases. The
// caller owns the parallel Host metadata and must permute it by these indices.
std::vector<std::size_t> normalized_trimmed_haplotype_indices(
    const std::vector<TrimmedHaplotypePopulationEntry>& haplotypes);

// Host-only EventMap construction from one haplotype/reference CIGAR.
// `reference_start` is the CIGAR offset within the supplied complete reference
// string (not the genomic origin of a sliced reference span). The full graph
// caller obtains the CIGAR from a Kokkos-scored SW request and passes it
// through this same state transition; keeping this focused entry point makes
// GATK's compound-event rules directly regression-testable.
HaplotypeEventMap event_map_from_haplotype_cigar(
    const std::string& reference, const std::string& haplotype,
    std::int32_t tid, std::size_t reference_start,
    const std::string& cigar, std::uint32_t max_mnp_distance);

// The Host representation of one GATK AssemblyResult.  GATK does not merge
// arbitrary local graphs into a synthetic contig-wide graph: each
// AssemblyRegion owns its immutable haplotype population and its EventMap,
// then each sample scores reads against that exact regional population.  The
// graph and variable-length alleles remain Host data; Kokkos is still used
// for the graph and PairHMM numerical kernels that consume them.
struct AssemblyRegionAssembly {
    fastgatk::kernels::ActivityRegion region;
    fastgatk::kernels::KmerGraphResult graph;
    // Preserves the graph→haplotype→EventMap ownership GATK keeps inside one
    // AssemblyResultSet. `event_map` below is the region's later union used
    // for allele replay; it must not erase the original per-haplotype maps.
    std::vector<HaplotypeEventMap> haplotype_event_maps;
    std::vector<ForcedAllele> event_map;
};

// A resolved Mutect2 feature predicate for one 0-based pileup locus. VCF
// parsing and contig-name resolution remain in the C++ Host tool; the shared
// calling pipeline receives only the compact facts that Mutect2Engine.isActive
// needs before it dispatches the Kokkos ActivityProfile kernel.
struct SomaticActivityFeatureMask {
    std::int32_t tid = -1;
    std::int32_t position = -1;
    // Bit i is set when a common germline substitution has reference base
    // A/C/G/T encoded by i. The source compares that base against the
    // dominant tumor PileupQualBuffer substitution bucket.
    std::uint8_t common_germline_substitution_reference_bases = 0;
    bool common_germline_indel = false;
    // PON membership has no allele-specific condition in Mutect2Engine.
    bool panel_of_normals = false;
};

// This is deliberately a small, composable caller contract.  It represents a
// single assembly region and is independent of HTSlib file lifetime, so it can
// later be used by a dispatcher that feeds one AssemblyRegion at a time.
struct Options {
    // The shared low-level API keeps its historical bounded-smoke defaults;
    // the HC CLI overrides these to one because GATK itself has no hard
    // pileup depth/ALT-support gate and delegates confidence to AF genotyping.
    std::size_t min_depth = 4;
    std::size_t min_alt_support = 3;
    std::uint8_t min_base_quality = 10;
    // ReferenceConfidenceModel uses this quality for a CIGAR deletion at the
    // Ref-vs-Any boundary.  It is distinct from the fixed Q45 indel-error
    // cache used to choose the conservative no-indel confidence model.
    std::uint8_t reference_model_deletion_quality = 30;
    // PairHMM input-quality controls mirror GATK's
    // LikelihoodEngineArgumentCollection and are separate from the pileup
    // candidate support threshold above.
    std::uint8_t pairhmm_base_quality_score_threshold = 18;
    bool pairhmm_disable_cap_base_qualities_to_mapq = false;
    // Strict compatibility callers use the direct Float64 recurrence.  The
    // HC CLI keeps this true unless the native-only optional boolean is
    // explicitly supplied as false; that opt-in path selects GKL-compatible
    // Float32 plus low-likelihood fallback for the regular base-space Kokkos
    // PairHMM.  GATK's FlowBased likelihood engine does not consume this
    // native PairHMM setting, so flow-tagged batches retain their calibrated
    // double flow-space kernel in either mode.
    bool pairhmm_use_double_precision = true;
    // GATK's default PCR_INDEL_MODEL=CONSERVATIVE turns the no-tag Q45
    // insertion/deletion default into an effective Q40 penalty.  NONE is
    // retained as an explicit compatibility mode for callers that request
    // the legacy Q45 behavior.
    bool pairhmm_conservative_indel_model = true;
    // Full GATK PCRErrorModel rate-factor selection.  The legacy boolean above
    // remains part of the public Host contract; when it is false the effective
    // model is NONE regardless of this factor.  Values are NONE=0, HOSTILE=1,
    // AGGRESSIVE=2 and CONSERVATIVE=3.
    std::string pairhmm_pcr_indel_model = "CONSERVATIVE";
    double pairhmm_pcr_error_rate_factor = 3.0;
    // GATK ReadLikelihoodCalculationEngine disqualifies a read when the
    // best haplotype likelihood is below ceil(read_length * rate) * -4.0,
    // capped at two expected errors.  Keep the policy in Options so HC and
    // Mutect2 share the same Host-side evidence contract.
    double pairhmm_expected_error_rate_per_base = 0.02;
    bool pairhmm_filter_poorly_modeled_reads = true;
    // GATK's `--phred-scaled-global-read-mismapping-rate` controls the
    // per-read likelihood floor applied after PairHMM.  A non-negative value
    // is phred-scaled (45 -> a 4.5 log10 cap); a negative value disables the
    // cap while retaining the best-likelihood disqualification calculation.
    double phred_scaled_global_read_mismapping_rate = 45.0;
    // GATK flow assembly's long-homopolymer collapse threshold.  Zero keeps
    // the historical no-collapse path; positive values apply the deterministic
    // Host collapse before flow-haplotype key encoding; -1 derives the limit
    // from the maximum read-group mc/flow_max_hmer in the current batch.
    std::int32_t flow_assembly_collapse_hmer_size = 0;
    // GATK's partial-mode restoration stops a deleted run at the first
    // reference base that no longer matches the restored homopolymer.
    bool flow_assembly_collapse_partial_mode = false;
    // AssemblyRegionTrimmer's default genotyping padding.  The Host chooses
    // the indel padding when a local candidate changes reference span.
    std::uint32_t snp_padding_for_genotyping = 20;
    // AssemblyRegionArgumentCollection.indelPaddingForGenotyping defaults to
    // 75 in GATK 4.6.2.0.  This is deliberately distinct from the 150 bp
    // historical native window: using that larger value retains extra reads
    // at a trimmed indel region and changes the PairHMM request ownership.
    std::uint32_t indel_padding_for_genotyping = 75;
    // AssemblyRegionTrimmer uses this base padding, plus the length of the
    // longest inferred repeat, for an indel that is a tandem repeat.  Keep it
    // separate from ordinary indel padding: the variable extension belongs to
    // Host-side AssemblyRegion/CIGAR ownership, not the Kokkos PairHMM kernel.
    std::uint32_t str_padding_for_genotyping = 75;
    // AssemblyBasedCallerArgumentCollection.informativeReadOverlapMargin
    // controls the VariantContext interval used by retainEvidence() after
    // PairHMM marginalization.  It is deliberately a Host-side interval
    // ownership rule: changing it must not widen the Kokkos read×haplotype
    // matrix or change graph construction.  GATK 4.6.2.0 defaults to 2 bp.
    std::uint32_t informative_read_overlap_margin = 2;
    // ActivityProfile/AssemblyRegion controls.  Keep the probability signal,
    // padding and propagation distance separate so a caller can bound the
    // graph/PairHMM window without changing which active loci are grouped.
    double active_probability_threshold = 0.002;
    // Mutect2 --force-active makes every ActivityProfile segment eligible for
    // AssemblyRegion graph/PairHMM processing while retaining profile
    // boundaries.  HC and direct callers retain the normal activity mask.
    bool force_active = false;
    // Mutect2's pileup active-region detector uses a quality-aware
    // log-likelihood ratio rather than HC's raw mismatch fraction.  The
    // native CLI enables this path and passes the GATK log10 threshold here;
    // direct callers retain the count-only default for API compatibility.
    bool activity_quality_aware_somatic = false;
    double activity_initial_tumor_log10_odds = 2.0;
    // Mutect2's --pcr-snv-qual is used both by the quality-aware activity
    // profile and by cleanOverlappingReadPairs, which receives half this
    // Phred value.  HC leaves the explicit optional cap absent and therefore
    // retains FragmentUtils' Q20 default.
    std::uint8_t activity_pcr_snv_quality = 40;
    // Mutect2 likewise passes --pcr-indel-qual / 2 to
    // cleanOverlappingReadPairs for BI/BD.  Keep it distinct from the
    // PairHMM PCR repeat error model: they are separate GATK arguments and
    // operate at different Host boundaries.
    std::uint8_t overlapping_pcr_indel_quality = 40;
    std::uint8_t activity_multiple_substitution_quality_correction = 5;
    // A joint Mutect2 AssemblyRegion traversal labels tumor and matched-normal
    // records through ReadBatch::sample_ids.  Negative values retain the
    // ordinary single-sample activity path.  When present, Host builds the
    // tumor PileupQualBuffer input and normal suppression gate before the
    // Kokkos LOD reduction.
    std::int32_t activity_tumor_sample_id = -1;
    std::int32_t activity_normal_sample_id = -1;
    // Multi-sample Mutect2 keeps sample identity through read-pair cleanup
    // and PairHMM, but isActive() applies roles: every non-normal sample is
    // tumor and every selected normal contributes to the normal gate.  These
    // Host-side role sets supersede the single-ID compatibility fields above
    // when non-empty.
    std::vector<std::int32_t> activity_tumor_sample_ids;
    std::vector<std::int32_t> activity_normal_sample_ids;
    // GATK's --genotype-germline-sites also disables the matched-normal
    // early activity veto; the regular somatic genotyper later decides which
    // alleles to write. Keep that source switch at the Host activity boundary.
    bool activity_genotype_germline_sites = false;
    // Host-resolved GATK feature inputs for Mutect2's pre-assembly activity
    // veto. Empty for HC and for Mutect2 invocations without an applicable
    // germline resource/PON predicate.
    std::vector<SomaticActivityFeatureMask> activity_feature_masks;
    // GATK 4.6.2.0 AssemblyRegionArgumentCollection defaults to a 100-base
    // halo. Keep the native default identical so small
    // --max-assembly-region-size windows have the same read ownership.
    std::uint32_t assembly_region_padding = 100;
    std::uint32_t min_assembly_region_size = 50;
    std::uint32_t max_assembly_region_size = 300;
    // BandPassActivityProfile's release default is 50 bases. A larger native
    // default changes active-region boundaries when the AssemblyRegion cap is
    // exercised, which in turn changes PairHMM evidence at shard edges.
    std::uint32_t max_probability_propagation_distance = 50;
    // Candidate discovery is owned by each GATK AssemblyRegion; it has no
    // caller-wide 64-site ceiling.  A global default cap silently discards
    // later regions in a real interval before they can be assembled (for
    // example the DREAM chr20 somatic calls beyond the first 12 kb).  Keep an
    // explicit native caller limit available as a resource guard, but make
    // the compatibility default unbounded and rely on the per-region
    // haplotype/link budgets below for bounded PairHMM work.
    std::size_t max_candidates = std::numeric_limits<std::size_t>::max();
    // GATK HaplotypeCaller --max-mnp-distance.  A value of zero keeps
    // substitutions as independent SNPs (the GATK default); positive values
    // allow graph-haplotype substitution runs whose consecutive mismatches
    // are at most this many reference bases apart to be emitted as one MNP.
    std::uint32_t max_mnp_distance = 0;
    // Bounded graph construction/traversal controls.  They are carried in the
    // Host Options object so HC/Mutect2 and future callers share one graph
    // policy instead of hard-coding device-independent constants.
    std::uint32_t graph_kmer_size = 5;
    std::uint32_t graph_min_kmer_count = 1;
    std::uint32_t graph_min_pruning = 2;
    // ReadThreadingAssemblerArgumentCollection.maxNumHaplotypesInPopulation
    // defaults to 128 in GATK 4.6.2.0.  This is particularly important for a
    // joint tumor/normal graph: a lower cap can silently drop the normal
    // competing path before either sample reaches PairHMM.
    std::size_t graph_max_paths = 128;
    // GATK K-best traversal ends at reference sink vertices, with no
    // fixed-bases default.  Zero carries that source-compatible behaviour;
    // a positive --max-haplotype-depth is an explicit resource guard.
    std::size_t graph_max_depth = 0;
    // GATK ReadThreadingAssembler-shaped dangling-branch controls.  Short
    // branches that do not rejoin the reference are pruned by default;
    // recovery keeps them as unanchored graph telemetry only.
    std::uint32_t graph_min_dangling_branch_length = 4;
    bool graph_recover_all_dangling_branches = false;
    bool graph_dont_increase_kmer_sizes_for_cycles = false;
    // GATK's default -1 selects the legacy mismatch-bounded dangling-head
    // recovery; non-negative values select the exact-prefix gate.
    std::int32_t graph_min_dangling_matching_bases = -1;
    // Match ReadThreadingAssembler's default: reject a k-mer attempt when the
    // reference window has repeated k-mers, then let Host try a larger k.
    bool graph_allow_non_unique_kmers_in_ref = false;
    std::uint32_t graph_num_pruning_samples = 1;
    bool graph_use_adaptive_pruning = false;
    double graph_initial_error_rate_for_pruning = 0.001;
    // GATK ReadThreadingAssembler uses natural-log thresholds
    // log10ToLog(1) and log10ToLog(4), not literal 2/3.  In particular the
    // latter is ln(10^4); lowering it incorrectly preserves singleton error
    // chains and can turn them into spurious dangling-head recoveries.
    double graph_pruning_log_odds_threshold = 2.30258509299404568402;
    double graph_pruning_seeding_log_odds_threshold = 9.21034037197618273607;
    std::uint32_t graph_max_unpruned_variants = 100;
    bool graph_linked_de_bruijn = false;
    bool graph_disable_artificial_haplotype_recovery = false;
    bool graph_enable_legacy_cycle_detection = false;
    // GATK --dont-use-soft-clipped-bases controls assembly evidence only.
    // Aligned bases and explicit CIGAR I/D events remain eligible; terminal
    // and internal soft-clip-derived haplotypes are excluded before graph
    // construction so HC and Mutect2 share the same policy.
    bool use_soft_clipped_bases = true;
    // GATK HaplotypeCaller --soft-clip-low-quality-ends changes only the
    // finalizeRegion representation.  Terminal Q<=minTailQuality bases stay
    // in the AssemblyRegion read as CIGAR soft clips, so the read-threading
    // graph still observes their sequence where its independent Q10 gate
    // permits it.  PairHMM subsequently hard-clips those new S operations;
    // this flag must therefore never bypass the numeric PairHMM read span.
    bool soft_clip_low_quality_ends = false;
    // Opt-in reference-confidence read-to-haplotype realignment.  The Host
    // GATK realigns every best-haplotype read before reference-confidence
    // calculation.  Preserve that as the native default: CIGAR composition
    // and sparse coordinate projection remain Host work, while PairHMM and
    // RCM numerical reductions remain Kokkos kernels.
    bool use_haplotype_realignment_for_rcm = true;
    // Optional GATK NearbyKmerErrorCorrector-shaped preprocessing.  It is
    // deliberately disabled by default and affects assembly/graph reads only;
    // original reads remain the evidence source for PairHMM and RCM.
    bool error_correct_reads = false;
    std::uint32_t error_correction_kmer_length = 25;
    std::uint32_t error_correction_min_solid_observations = 20;
    // Finite values select GATK PileupReadErrorCorrector and take precedence
    // over NearbyKmerErrorCorrector, matching AssemblyBasedCallerUtils.  The
    // default -infinity keeps both correction paths disabled.
    double pileup_error_correction_log_odds = -std::numeric_limits<double>::infinity();
    // Maximum number of nearby candidate alleles for which the Host builds
    // one complete 2^N local haplotype combination set.  Larger groups are
    // partitioned into deterministic coordinate blocks of this size; each
    // block still receives its complete combination set instead of silently
    // collapsing the whole region to singleton alternatives.
    std::size_t max_haplotype_combination_alleles = 6;
    // After PairHMM scoring, graph-derived haplotypes more than this many
    // log10 likelihood units below the best local haplotype are pruned from
    // genotype marginalization.  A best graph path is always retained so a
    // graph-derived alternate cannot silently disappear from telemetry.
    double haplotype_pruning_log10 = 20.0;
    // GATK genotype priors used by AF/model calculations.  Default HC
    // genotype assignment remains USE_PLS_TO_ASSIGN; prior-aware GT/GQ is an
    // explicit USE_POSTERIOR_PROBABILITIES compatibility mode.
    bool use_genotype_priors = true;
    bool use_posterior_genotype_assignment = false;
    double heterozygosity = 1.0e-3;
    // HaplotypeCaller inherits HomoSapiensConstants.INDEL_HETEROZYGOSITY:
    // 1 / 8000.  This feeds the AF Dirichlet pseudocount for both concrete
    // indels and EventMap spanning-deletion (`*`) alleles.
    double indel_heterozygosity = 1.0 / 8000.0;
    // GATK's AlleleFrequencyCalculator uses this spread when converting
    // heterozygosity into Dirichlet pseudocounts.  Keep it at the Host
    // boundary so HaplotypeCaller call-confidence filtering and
    // GenotypeGVCFs use the same prior model on every execution space.
    double heterozygosity_stdev = 0.01;
    // Standard phred-scaled call threshold.  GATK defaults to 30 for a
    // normal VCF and forces it to zero for -ERC GVCF/BP_RESOLUTION; the HC
    // driver applies that mode-specific override after parsing.
    double standard_confidence_for_calling = 30.0;
    // Mutect2 emits candidate loci for downstream somatic filtering even when
    // the HC-oriented genotype confidence/read-haplotype gate has no retained
    // ALT owner.  This flag keeps that policy explicit instead of inferring it
    // from a zero confidence threshold.
    bool somatic_mode = false;
    // Mutect2's --independent-mates changes only the Host evidence identity
    // passed to AlleleLikelihoods.groupEvidence: each read becomes its own
    // one-read Fragment instead of sharing the read-name fragment. Graph
    // construction and the Kokkos PairHMM/somatic reductions are unchanged.
    bool somatic_independent_mates = false;
    // Sample ploidy used by both VCF and reference-confidence genotype rows.
    // The Host writer bounds this value to keep Number=G workspace finite.
    int sample_ploidy = 2;
    // GenotypeCalculationArgumentCollection.maxGenotypeCount (default 1024).
    // HaplotypeCaller applies this to each EventMap allele map before
    // read-likelihood marginalization; the Host owns that variable-length
    // structural subsetting while the remaining numeric likelihood work stays
    // in Kokkos.
    std::size_t max_genotype_count = 1024;
    // Maximum indel size considered by ReferenceConfidenceModel when deciding
    // whether a read is informative for a hom-ref call. GATK's default is 10.
    int indel_size_to_eliminate_in_ref_model = 10;
    // Apply the deterministic Kokkos read-filter mask to all downstream
    // stages.  The Host still retains the input count for telemetry; callers
    // can opt out explicitly for compatibility/debug comparisons.
    bool apply_read_filters = true;
    // GATK walkers first apply their complete read-filter list and only then
    // run any tool-specific post-filter transformer.  The callback keeps
    // that variable-length Host transformation at the exact boundary while
    // the downstream graph, likelihood and reduction kernels stay unchanged.
    // It returns the number of records transformed for Result telemetry.
    std::function<std::size_t(io::ReadBatch&)> post_read_filter_transform;
    // GATK's --disable-tool-default-read-filters removes only the tool's
    // implicit filter set; explicitly requested --read-filter predicates must
    // remain active.  CLI adapters resolve that distinction before entering
    // the shared pipeline and retain this bit for manifest provenance.
    bool disable_tool_default_read_filters = false;
    // GATK HaplotypeCaller/Mutect2 default read-filter policy.  The native
    // filter kernel accepts these as data so the same mask is consumed by
    // assembly, graph, SW/PairHMM and RCM rather than silently applying a
    // different threshold in one stage.  Neither caller's standard filter
    // list includes NotSupplementaryAlignmentReadFilter; an explicit
    // --read-filter can still enable that predicate at the CLI boundary.
    std::uint8_t minimum_mapping_quality = 20;
    bool exclude_mapping_quality_unavailable = false;
    bool exclude_mapping_quality_zero = false;
    bool exclude_duplicates = true;
    bool exclude_unmapped = true;
    bool exclude_secondary = true;
    bool exclude_supplementary = false;
    bool exclude_qcfail = true;
    bool require_good_cigar = false;
    bool require_nonzero_reference_span = false;
    bool require_non_chimeric_original_alignment = false;
    // Components of GATK WellformedReadFilter used by production HC/Mutect2
    // defaults. Direct library callers retain permissive defaults.
    bool require_no_n_cigar = false;
    bool require_read_group = false;
    bool require_read_length = false;
    std::uint32_t min_read_length = 1;
    std::uint32_t max_read_length = std::numeric_limits<std::uint32_t>::max();
    // GATK HaplotypeCaller cleans overlapping paired fragments before
    // assembly.  The default caps concordant overlapping base qualities at
    // half the default PCR SNV quality (Phred 20) and zeros conflicts.  Keep
    // the inverse flag explicit for command-line compatibility/debug oracles.
    bool do_not_correct_overlapping_base_qualities = false;
    // AssemblyRegionWalker applies a positional reservoir downsampler after
    // its read filters. HaplotypeCaller and Mutect2 default to 50 reads at
    // one alignment start; zero explicitly disables that GATK control. This
    // remains Host ownership before assembly/activity/PairHMM consume reads.
    std::size_t max_reads_per_locus = 50;
    std::uint64_t downsampling_seed = 47382911ULL;
    // Mutect2's dedicated MutectDownsampler keeps only MAPQ > 50 reads when
    // a downsampling pool exceeds the coverage budget.  Its pool begins at
    // the first observed alignment start and spans `stride` bases; an
    // optional suspicious-read limit rejects the complete pool. HC retains
    // the ordinary AssemblyRegionWalker reservoir policy.
    bool mutect2_downsampling = false;
    std::uint32_t mutect2_downsampling_stride = 1;
    std::int32_t mutect2_max_suspicious_reads_per_alignment_start = 0;
    // Optional 0-based interval supplied by the Host reader.  When present,
    // gVCF reference confidence is emitted across the complete interval,
    // including zero-depth blocks.  Whole-contig GVCF mode uses
    // emit_reference_confidence with -1 interval coordinates.
    std::int32_t interval_tid = -1;
    std::int32_t interval_start = -1;
    std::int32_t interval_end = -1;  // exclusive
    // Repeatable -L/--intervals selectors normalized by the Host reader. When
    // non-empty these take precedence over the legacy single-interval fields.
    std::vector<io::HtsInterval> intervals;
    // Repeatable -XL/--exclude-intervals selectors. Exclusions are applied
    // after include selection at the projected-locus boundary, preserving
    // evidence from reads that cross an excluded and included locus.
    std::vector<io::HtsInterval> exclusion_intervals;
    // Optional wider decode/assembly domain used by region streaming.  The
    // public `intervals` remain the output core, while this halo domain keeps
    // reads crossing a tile boundary available to activity/graph/PairHMM.
    // Empty means use `intervals` for both domains.
    std::vector<io::HtsInterval> locus_intervals;
    // Set by the shared caller when this invocation is already one
    // AssemblyRegion unit.  A contiguous multi-region interval is split into
    // independent assemble→likelihood→genotype passes; the recursive unit
    // must not split again, even if its padded window still contains more
    // than one activity island.
    bool assembly_region_independent_pass = false;
    // Experimental Mutect2 parity path: execute each active AssemblyRegion
    // through its complete assemble→EventMap→PairHMM pipeline before merging
    // the resulting per-region calls.  This must travel with the complete
    // downstream pipeline; merely splitting the compact k-mer graph permits
    // cross-region semantics to leak back through candidate/Pairs ownership.
    bool partition_somatic_assembly_regions = false;
    // When set, this invocation is one GATK ActiveRegion: do not re-cluster
    // activity into multiple calling windows.  PairHMM/graph see exactly this
    // padded region so independent units report one real group each.
    bool force_single_calling_region = false;
    fastgatk::kernels::ActivityRegion forced_calling_region{};
    // Optional Host-provided ActiveRegion partition.  Mutect2 determines its
    // AssemblyRegions from the tumor traversal, then evaluates matched-normal
    // evidence in those same graph/PairHMM windows.  A non-empty vector
    // replaces the locally discovered calling windows while retaining this
    // invocation's own reads and likelihood calculations.
    std::vector<fastgatk::kernels::ActivityRegion> forced_calling_regions;
    // A Mutect2 AssemblyResultSet is built from the union of every sample's
    // reads, then the very same haplotype population is scored separately
    // for each sample.  Supplying this value replaces only the locally built
    // graph/haplotype population; read filtering, clipping and Kokkos
    // PairHMM scoring still run against this invocation's own ReadBatch.
    // Keeping the graph as an explicit Host input makes that ownership
    // auditable and avoids re-assembling normal-only paths.
    std::optional<fastgatk::kernels::KmerGraphResult> shared_assembly_graph;
    // Region-scoped form of shared_assembly_graph.  When populated, run()
    // executes one complete likelihood/genotyping pass per supplied
    // AssemblyResult, with that result's EventMap and graph kept together.
    // This is the Mutect2Engine.callRegion ownership boundary.
    std::vector<AssemblyRegionAssembly> shared_assembly_regions;
    // Preserve the local AssemblyResult objects when partitioning an assembly
    // traversal.  This is intentionally opt-in because normal HC callers
    // only need flattened telemetry/calls, whereas Mutect2 needs to replay
    // the joint region objects for each sample's PairHMM pass.
    bool retain_assembly_region_results = false;
    // Build and retain only the Host-owned AssemblyResult contract
    // (ActivityRegion, graph and EventMap), then stop before pileup
    // likelihoods, PairHMM and genotyping.  Mutect2 uses this for the joint
    // tumor/normal assembly traversal: GATK builds that shared graph once and
    // scores each sample separately, so a joint PairHMM matrix has no
    // downstream consumer and must not be computed merely to obtain the
    // AssemblyResultSet.  This is an internal scheduling mode, not a CLI
    // approximation; the actual numerical likelihood work remains on the
    // Kokkos path during the subsequent per-sample replays.
    bool assembly_only = false;
    // In addition to graph/EventMap snapshots, retain the completed local
    // likelihood Result for each replayed AssemblyResult.  This is required
    // only by the Mutect2 Host reducer: SomaticGenotypingEngine evaluates
    // TLOD before regional results are merged, while ordinary callers keep
    // the flattened Result fast path.
    bool retain_assembly_region_likelihood_results = false;
    // Additional tumor-derived alleles to score in this sample. They remain
    // distinct from locally assembled candidates, but are included in the
    // local PairHMM haplotype set even at zero local ALT support.
    std::vector<ForcedAllele> forced_alleles;
    // Use only forced alleles (and their local combination haplotypes),
    // discarding branches assembled from this sample's reads. This is useful
    // only once the caller has supplied the complete shared AssemblyResultSet.
    // GATK Mutect2 assembles all samples in one AssemblyRegion before it
    // calculates either sample's likelihoods; a separate normal invocation
    // must retain local paths until native represents that joint assembly.
    bool restrict_to_forced_alleles = false;
    // When set, extra assembled haplotype *sequences* stay out of PairHMM
    // (chr17 --dont-use-soft-clipped-bases QUAL).  Haplotype-vs-reference
    // variant sites still join the genotyping candidate set (GATK HC step 2).
    bool omit_graph_haplotype_injection = false;
    // HC Host sets this only for -ERC GVCF.  With no explicit interval, the
    // reference-backed path emits depth-0 blocks for uncovered contig spans;
    // without a reference it falls back to observed loci.
    bool emit_reference_confidence = false;
    // GATK -ERC BP_RESOLUTION emits one REF/<NON_REF> record per reference
    // base instead of merging hom-ref loci into GVCF quality bands.
    bool emit_reference_confidence_bp_resolution = false;
    // Optional GVCF genotype-quality partition upper bounds from GATK's
    // repeatable -GQB/--GVCF-GQ-bands option.  An empty vector preserves
    // HaplotypeCaller's default partitions (0-1, ..., 60-70, ..., 99-100).
    // When non-empty, the Host CLI appends the required terminal 100 bound.
    std::vector<int> gvcf_gq_bands;
    // GATK ReadThreadingAssembler accepts --kmer-size repeatedly (the
    // release default is [10,25]). Keep graph_kmer_size as the historical
    // scalar ABI field and append the optional ordered request list so
    // callers can preserve the full Java option semantics without changing
    // existing aggregate initialization of this public struct.
    std::vector<std::uint32_t> graph_kmer_sizes;
};

struct AssemblyCandidate {
    std::int32_t tid = -1;
    std::int32_t position = -1;  // 0-based reference coordinate
    std::uint8_t reference = 'N';
    std::uint8_t alternate = 'N';
    // Empty strings mean a single-base SNP represented by reference/alternate
    // above.  Non-empty strings carry an anchored indel allele pair for VCF,
    // SW and PairHMM; position is the left anchor coordinate.
    std::string reference_allele;
    std::string alternate_allele;
    std::uint32_t reference_count = 0;
    std::uint32_t alternate_count = 0;
    // FORMAT/DP is the retained genotyping evidence depth.  INFO/DP follows
    // GATK Coverage and counts all reads overlapping the variant's annotation
    // interval (including reads later disqualified by PairHMM).
    std::uint32_t depth = 0;
    std::uint32_t variant_depth = 0;
    // True when the allele was materialized directly from an input CIGAR I/D
    // event.  PairHMM's base-oriented evidence counter may be zero for these
    // alleles, so the writer can preserve the independent assembly support;
    // graph/soft-clip-only candidates deliberately leave this false.
    bool cigar_indel = false;
    // Mutect2-only recovery for a locally assembled SNP backed by at least
    // two independent raw observations below the pileup BQ gate. PairHMM
    // remains the authority for retaining and emitting this candidate.
    bool somatic_low_quality_assembly_recovery = false;
    // Optional assembled-path context for the recovery above. The insertion
    // is not an independently emitted allele, but belongs to the alternate
    // haplotype consumed by PairHMM in both tumor and matched-normal runs.
    std::int32_t linked_insertion_reference_position = -1;
    std::string linked_insertion_bases;
    // Set only for a tumor-derived allele injected into matched-normal
    // scoring. Such an allele must be retained in the local haplotype set
    // even when normal pileup support is below candidate discovery thresholds.
    bool forced_by_tumor = false;
    // Set for a concrete HaplotypeCaller --alleles FeatureInput event.  GATK
    // GenotypeGivenAlleles retains and emits this event even when PairHMM
    // correctly assigns every observation to the reference haplotype.
    // Keeping the provenance separate from the Mutect2 replay flag prevents
    // the two callers' output policies from leaking into one another.
    bool forced_by_alleles_feature = false;
    // True when this exact allele was projected from a retained local graph
    // haplotype. Default HC admits only this EventMap-derived set inside its
    // AssemblyRegion; separately-scoped somatic/compatibility paths may keep
    // their explicitly documented non-graph evidence.
    bool graph_derived = false;
    // GATK's LongHomopolymerHaplotypeCollapsingEngine marks an EventMap event
    // with INFO/XC when its source haplotype restored a long reference HMER.
    // This is Host-only EventMap provenance, never a PairHMM/kernel input.
    bool flow_hmer_collapsed = false;
};

struct Likelihoods {
    double hom_ref = 0.0;
    double het = 0.0;
    double hom_alt = 0.0;
};

struct GenotypeCall {
    AssemblyCandidate candidate;
    Likelihoods likelihoods;
    std::uint8_t genotype = 0;  // 0=0/0, 1=0/1, 2=1/1
    std::uint8_t gq = 0;
    // Site-level AF posterior quality used for call-thresholding and VCF
    // QUAL.  This is distinct from FORMAT/GQ (genotype-vs-next-genotype).
    double qual = 0.0;
    // Standard GATK VariantContext annotations computed at the Host evidence
    // boundary.  Missing values remain NaN and are omitted from the VCF;
    // device kernels never receive these string/report concerns.
    struct Annotations {
        double mq = std::numeric_limits<double>::quiet_NaN();
        double qd = std::numeric_limits<double>::quiet_NaN();
        double fs = std::numeric_limits<double>::quiet_NaN();
        double sor = std::numeric_limits<double>::quiet_NaN();
        // Mann-Whitney rank-sum annotations are retained as numeric Host
        // evidence.  They are intentionally separate from MQ/QD/FS/SOR:
        // GATK's rank tests use the informative-read likelihood partition,
        // while this bounded native path uses the CIGAR-projected allele
        // evidence available at the output boundary.
        double mq_rank_sum = std::numeric_limits<double>::quiet_NaN();
        double read_pos_rank_sum = std::numeric_limits<double>::quiet_NaN();
        double base_q_rank_sum = std::numeric_limits<double>::quiet_NaN();
        // GATK's gVCF candidate rows carry the raw MQ/DP accumulator and the
        // four-cell strand-bias table in FORMAT/SB.  These are Host evidence
        // counters (never device objects) retained alongside the scalar
        // annotations so the writer can reproduce that contract.
        std::uint64_t raw_mq_sum_square = 0;
        std::uint32_t raw_mq_depth = 0;
        std::uint32_t ref_forward = 0;
        std::uint32_t ref_reverse = 0;
        std::uint32_t alt_forward = 0;
        std::uint32_t alt_reverse = 0;
    } annotations;
};

// A multi-ALT FORMAT/AD row is materialized once at the Host read/projection
// boundary. Keeping the partition in Result lets the VCF writer use the same
// per-read allele ownership after region streaming or result merging, without
// passing HTSlib/Kokkos objects into the serializer.
struct MultiallelicDepth {
    std::int32_t tid = -1;
    std::int32_t position = -1;
    std::string reference;
    std::vector<std::string> alternates;
    std::vector<std::uint32_t> depths;
    std::uint32_t depth = 0;
};

struct ReferenceBlock {
    std::int32_t tid = -1;
    std::int32_t start = -1;  // inclusive, 0-based
    std::int32_t end = -1;    // inclusive, 0-based
    std::uint8_t reference = 'N';
    std::uint32_t depth = 0;
    // Reference-confidence fields are carried with the block representative.
    // They are computed from the same filtered/projected observations used by
    // assembly.  A block is merged only when its confidence band is stable;
    // the fields therefore remain valid after GVCF state-machine coalescing.
    std::uint32_t reference_count = 0;
    std::uint32_t non_ref_count = 0;
    std::uint8_t gq = 0;
    std::uint8_t gq_band = 0;
    std::array<std::int32_t, 3> pl{0, 0, 0};
    // Arbitrary-ploidy REF/<NON_REF> Number=G rows.  Diploid callers retain
    // `pl` for the established ABI and fill this vector with the same values.
    std::vector<std::int32_t> genotype_pl;
    // GATK HomRefBlock stores every per-site DP so the emitted FORMAT/DP is
    // the rounded median, while FORMAT/MIN_DP is the minimum.  Keep the
    // bounded block-local samples in Host memory; device kernels never see
    // this writer state.
    std::vector<std::uint32_t> depth_values;
    std::uint32_t min_depth = 0;
};

// Mutect2 reference confidence is not an HC GQ/PL block: the source
// SomaticReferenceConfidenceModel reduces every post-finalization pileup to
// a single <NON_REF> TLOD, and SomaticGVCFWriter subsequently bands those
// values.  Preserve the per-locus Host representation until that writer has
// made its TLOD-band decisions; no variable-length VCF state crosses into a
// Kokkos kernel.
struct SomaticReferenceConfidenceLocus {
    std::int32_t tid = -1;
    std::int32_t position = -1;  // 0-based
    std::uint32_t depth = 0;
    std::uint32_t reference_depth = 0;
    std::uint32_t non_reference_depth = 0;
    double log10_lod = 0.0;
};

// AssemblyRegion-scoped read realignment.  This is deliberately sparse: an
// AssemblyRegion only owns the reads it actually handed to PairHMM, so a
// dense candidate x all-input-read matrix both invents non-existent evidence
// slots and makes real genomic windows quadratic in memory.  Entries are
// sorted by (context_ordinal, source_record), where source_record is the
// compact likelihood-row ordinal (not Result::likelihood_read_source_records'
// outer-batch provenance ordinal).  An entry with a negative interval records
// that the local Host traceback fell back to the original read alignment.
struct CandidateReadRealignment {
    std::uint32_t context_ordinal = std::numeric_limits<std::uint32_t>::max();
    std::uint32_t source_record = std::numeric_limits<std::uint32_t>::max();
    std::int32_t tid = -1;
    std::int64_t start = -1;
    std::int64_t end = -1;
    // Inclusive source-read coordinate range that remains after the selected
    // haplotype-to-reference realignment.  Together with the original CIGAR
    // hard clips this recreates GATKRead's final CIGAR coordinate system for
    // ReadPosRankSum without retaining a variable-length CIGAR per read.
    std::uint32_t first_read_offset = std::numeric_limits<std::uint32_t>::max();
    std::uint32_t last_read_offset = std::numeric_limits<std::uint32_t>::max();
    bool qualified = false;
};

// One retained candidate/fragment likelihood row after
// AlleleLikelihoods.groupEvidence(...).marginalize(...).  Rows are kept in
// candidate-major CSR order because an AssemblyRegion usually overlaps only a
// small subset of all fragments in a genomic window.  The Host owns this
// irregular membership; Kokkos consumes the parallel numeric likelihood
// arrays for the somatic mixture reduction.
struct SomaticGroupLikelihood {
    std::uint32_t group = std::numeric_limits<std::uint32_t>::max();
    double reference = -1.0e299;
    double alternate = -1.0e299;
    // GATK's AssemblyBasedCallerUtils.createAlleleMapper(..., true) assigns
    // haplotypes whose EventMap deletion began before this locus to the
    // symbolic spanning-deletion allele (`*`).  It is only a Host-side
    // EventMap membership row; the somatic mixture reduction stays in the
    // existing Kokkos multiallelic kernel.
    double spanning_deletion = -1.0e299;
};

struct Result {
    std::size_t reads = 0;
    // `observations` preserves the historical smoke metric: all non-N input
    // bases.  `projected_observations` is the CIGAR-aware subset that maps to
    // one reference coordinate and is consumed by assembly kernels.
    std::size_t observations = 0;
    std::size_t projected_observations = 0;
    std::size_t loci = 0;
    std::size_t candidate_sites = 0;
    // Host-side assembly evidence gates.  These count pileup alleles that
    // were present in the primitive locus set but were rejected before
    // PairHMM because all support was adjacent to a soft-clip boundary or
    // collapsed to one overlapping read fragment.  Keeping the counters in
    // Result makes a compatibility audit distinguish assembly pruning from
    // genotype confidence filtering.
    std::size_t candidate_softclip_suppressed = 0;
    std::size_t candidate_fragment_suppressed = 0;
    // Pileup-only singleton alleles that do not meet the HC graph pruning
    // floor and are folded back after genotyping.  They remain in PairHMM so
    // neighboring candidate evidence cannot change with this output gate.
    std::size_t candidate_low_support_suppressed = 0;
    // Number of read bases remapped through a selected haplotype SW CIGAR
    // before reference-confidence calculation.  A zero value is expected
    // when the compatibility switch is disabled or no valid mapping exists.
    std::size_t rcm_realigned_observations = 0;
    std::size_t rcm_realignment_fallback_observations = 0;
    bool rcm_haplotype_realignment_used = false;
    // Number of somatic candidates whose pileup evidence was restored after
    // the reference-backed PairHMM ownership gate had no retained ALT owner.
    // This is an auditable compatibility fallback, not a claim of HC parity.
    std::size_t somatic_pileup_fallback_candidates = 0;
    std::size_t multiallelic_candidates_pruned = 0;
    std::size_t softclip_candidates = 0;
    std::size_t compound_cigar_candidates = 0;
    // Candidate assembly boundary accounting.  A candidate whose full REF
    // span crosses two disjoint AssemblyRegions is not safe to materialize in
    // either local graph and is suppressed; candidates outside all regions
    // remain an explicit contig-level fallback for compatibility.
    std::size_t assembly_unassigned_candidates = 0;
    std::size_t assembly_cross_region_candidates = 0;
    // Raw activity regions can overlap at their halo.  The calling path
    // coalesces such overlapping windows before graph/PairHMM materialization
    // so a compound allele whose reference span crosses the raw boundary is
    // not discarded merely because two local windows were emitted.
    std::size_t assembly_region_union_count = 0;
    std::size_t assembly_region_union_merges = 0;
    std::size_t assembly_cross_region_rescued = 0;
    bool assembly_region_partitioned = false;
    std::size_t graph_variant_candidates = 0;
    std::size_t graph_snp_candidates = 0;
    std::size_t graph_mnp_candidates = 0;
    std::size_t graph_deletion_candidates = 0;
    std::size_t variant_calls = 0;
    std::uint64_t signature = 0;
    double assembly_prepare_seconds = 0.0;
    double assembly_seconds = 0.0;
    double likelihood_prepare_seconds = 0.0;
    double likelihood_seconds = 0.0;
    double pairhmm_prepare_seconds = 0.0;
    double pairhmm_seconds = 0.0;
    std::size_t pairhmm_pairs = 0;
    // Candidate-to-haplotype associations before stable read×haplotype
    // request deduplication. `pairhmm_pairs` is the actual kernel request
    // count; this field exposes the avoided duplicate work.
    std::size_t pairhmm_request_links = 0;
    std::size_t pairhmm_haplotypes = 0;
    std::size_t pairhmm_graph_haplotypes = 0;
    std::size_t pairhmm_graph_snp_posterior_pairs = 0;
    std::size_t pairhmm_graph_haplotypes_considered = 0;
    std::size_t pairhmm_haplotypes_pruned = 0;
    std::size_t pairhmm_graph_haplotypes_pruned = 0;
    std::size_t pairhmm_haplotype_combination_blocks = 0;
    bool pairhmm_haplotype_combination_chunked = false;
    // PairHMM candidate/haplotype construction is partitioned by the Host
    // AssemblyRegion windows. These fields make the boundary explicit so a
    // caller can distinguish bounded region traversal from the contig-wide
    // fallback.
    std::size_t pairhmm_assembly_region_groups = 0;
    std::size_t pairhmm_unassigned_candidates = 0;
    bool pairhmm_assembly_region_partitioned = false;
    std::size_t pairhmm_marginalized_candidates = 0;
    std::size_t pairhmm_cached_shapes = 0;
    std::size_t pairhmm_cache_hits = 0;
    // GATK AlleleLikelihoods normalization (best-per-read plus global
    // mismapping cap) is a separate Kokkos plan from the PairHMM recurrence;
    // expose it so prepare/execute cost and backend are auditable.
    std::string pairhmm_normalization_execution_space;
    double pairhmm_normalization_prepare_seconds = 0.0;
    double pairhmm_normalization_seconds = 0.0;
    // Per-candidate/read REF-vs-ALT max-marginalization is a second shared
    // Kokkos plan; CIGAR and graph eligibility stay on Host, while the
    // numerical reduction is backend-independent.
    std::string pairhmm_marginalization_execution_space;
    double pairhmm_marginalization_prepare_seconds = 0.0;
    double pairhmm_marginalization_seconds = 0.0;
    // Best/second-best concrete allele reduction used for read uncertainty;
    // the Kokkos plan is shared by regular and flow PairHMM outputs.
    std::string pairhmm_uncertainty_execution_space;
    double pairhmm_uncertainty_prepare_seconds = 0.0;
    double pairhmm_uncertainty_seconds = 0.0;
    bool pairhmm_used = false;
    std::string pairhmm_precision = "double";
    // PairHMM model selected for the reference-backed path.  `regular` is
    // the default GATK-compatible base-space recurrence; `flow` is selected
    // only when every participating read carries a complete, consistent flow
    // tag set.  A run that executes both eligible flow and explicit regular
    // fallback partitions reports `mixed`; empty partitions do not affect the
    // executed-model telemetry.
    std::string pairhmm_error_model = "regular";
    // Zero denotes an uninitialized aggregate before the first batch is
    // merged; a completed run records the effective adjusted quality: 40
    // (Conservative) or 45 (NONE).
    std::uint8_t pairhmm_default_indel_quality = 0;
    // GATK's missing BI/BD fallback is Q45.  Conservative PCR adjustment
    // does not touch the terminal read position, so expose that raw terminal
    // quality separately from the effective adjusted-position scalar above.
    std::uint8_t pairhmm_default_indel_terminal_quality = 0;
    std::string pairhmm_pcr_indel_model = "CONSERVATIVE";
    double pairhmm_pcr_error_rate_factor = 0.0;
    // Number of read-position insertion/deletion penalties lowered by the
    // tandem-repeat PCR model before entering the Kokkos PairHMM kernel.
    std::size_t pairhmm_pcr_adjusted_positions = 0;
    std::size_t pairhmm_flow_reads = 0;
    std::size_t pairhmm_flow_haplotypes = 0;
    std::size_t pairhmm_flow_haplotypes_clipped = 0;
    // Retained output-schema compatibility counter.  GATK scores the raw
    // assembled FlowBasedHaplotype and does not collapse haplotypes before
    // flow-key encoding, so this remains zero; HMER restoration is recorded
    // by `pairhmm_flow_haplotypes_uncollapsed` below.
    std::size_t pairhmm_flow_haplotypes_collapsed = 0;
    // Number of flow haplotypes for which the post-likelihood reference
    // uncollapse pass restored one or more deleted reference HMER bases.
    // This transform occurs after Kokkos scoring, before downstream genotype
    // ownership is canonicalized.
    std::size_t pairhmm_flow_haplotypes_uncollapsed = 0;
    // Number of original flow haplotype entries remapped to a canonical
    // post-uncollapse sequence before AlleleLikelihoods marginalization.
    std::size_t pairhmm_flow_haplotype_remaps = 0;
    // Number of canonical sequence-identical groups containing at least two
    // original flow haplotype entries.
    std::size_t pairhmm_flow_identical_haplotype_groups = 0;
    std::size_t pairhmm_flow_reads_clipped = 0;
    std::size_t pairhmm_flow_clipping_fallbacks = 0;
    std::size_t pairhmm_insertion_quality_reads = 0;
    std::size_t pairhmm_deletion_quality_reads = 0;
    // AssemblyRegion hard-clipping is performed on Host using CIGAR/reference
    // coordinates before the read enters PairHMM.  Keep the counters explicit
    // so a caller can distinguish clipped-read semantics from the old
    // contig-wide read view.
    std::size_t pairhmm_reads_clipped = 0;
    std::size_t pairhmm_reads_dropped_after_clipping = 0;
    // GATK AssemblyRegion trimming drops fragments shorter than ten bases
    // before PairHMM. Expose the active contract for audit/fixture tests.
    std::size_t pairhmm_minimum_read_length_after_trimming = 10;
    std::size_t pairhmm_reads_disqualified = 0;
    double pairhmm_read_disqualification_threshold = 0.0;
    // Stable reason for a non-executed PairHMM path.  An empty request set is
    // a valid compatibility fallback (for example a pileup-only smoke input),
    // but it must never be indistinguishable from a kernel that was forgotten.
    std::string pairhmm_skip_reason;
    std::string pairhmm_execution_space;
    double sw_prepare_seconds = 0.0;
    double sw_seconds = 0.0;
    std::size_t sw_pairs = 0;
    std::size_t sw_simd_width = 1;
    std::size_t sw_simd_groups = 0;
    bool sw_used = false;
    std::string sw_execution_space;
    std::size_t sw_traceback_pairs = 0;
    std::size_t sw_indel_alignments = 0;
    std::uint64_t sw_cigar_signature = 0;
    double sw_traceback_seconds = 0.0;
    // ReadLikelihoodCalculationEngine provenance.  These counters are kept
    // separate from the generic SW fields so a manifest consumer can tell
    // that the CIGAR/score contract is for read x haplotype requests (the
    // exact matrix consumed by PairHMM), rather than a graph/reference-only
    // alignment.  Traceback is still Host-side; the score path is Kokkos.
    std::size_t read_haplotype_cigar_pairs = 0;
    std::size_t read_haplotype_cigar_valid_pairs = 0;
    std::size_t read_haplotype_cigar_indel_pairs = 0;
    // Candidate/read/haplotype associations whose local SW CIGAR does not
    // actually cover the candidate event are excluded from allele
    // marginalization.  The counters make the CIGAR-informed evidence gate
    // observable instead of presenting it as a hidden PairHMM heuristic.
    std::size_t read_haplotype_cigar_filtered_pairs = 0;
    std::size_t read_haplotype_softclip_filtered_pairs = 0;
    // GATK's BestAllele informative threshold is a 0.2 log10 likelihood
    // margin.  The native path reports the bounded top-two margin telemetry
    // and uses the same Kokkos owner reduction for one-read-one-allele
    // evidence depth/counts; full read-likelihood posterior parity remains
    // an explicit follow-up contract.
    std::size_t read_haplotype_uncertain_reads = 0;
    std::size_t read_haplotype_informative_reads = 0;
    double read_haplotype_margin_sum = 0.0;
    std::uint64_t read_haplotype_uncertainty_signature = 1469598103934665603ULL;
    bool graph_used = false;
    std::size_t graph_nodes = 0;
    std::size_t graph_edges = 0;
    std::size_t graph_branching_nodes = 0;
    std::size_t graph_reference_nodes = 0;
    std::size_t graph_reference_edges = 0;
    std::size_t graph_reference_connected_nodes = 0;
    std::size_t graph_dangling_nodes = 0;
    std::size_t graph_pruned_nodes = 0;
    std::size_t graph_dangling_branch_paths = 0;
    std::size_t graph_dangling_branch_bases = 0;
    std::size_t graph_dangling_recovered_paths = 0;
    std::size_t graph_dangling_recovered_bases = 0;
    std::size_t graph_artificial_haplotype_recovery_paths = 0;
    std::size_t graph_artificial_haplotype_recovery_bases = 0;
    std::size_t graph_seqgraph_nodes = 0;
    std::size_t graph_seqgraph_edges = 0;
    std::size_t graph_seqgraph_linear_chain_merges = 0;
    std::size_t graph_seqgraph_diamond_merges = 0;
    std::size_t graph_seqgraph_tail_merges = 0;
    std::size_t graph_seqgraph_suffix_splits = 0;
    std::size_t graph_seqgraph_suffix_merges = 0;
    std::size_t graph_adaptive_pruned_nodes = 0;
    std::size_t graph_non_unique_kmers = 0;
    std::size_t graph_reference_non_unique_kmers = 0;
    bool graph_reference_kmer_rejected = false;
    std::uint32_t graph_kmer_size_selected = 0;
    std::size_t graph_kmer_iterations = 1;
    bool graph_has_non_reference_cycles = false;
    std::size_t graph_reference_paths = 0;
    std::size_t graph_haplotype_paths = 0;
    std::vector<std::string> graph_haplotype_sequences;
    std::size_t graph_haplotype_sequence_count = 0;
    std::size_t graph_haplotype_sequence_max_length = 0;
    std::vector<std::int32_t> graph_haplotype_path_tids;
    std::vector<std::int32_t> graph_haplotype_path_starts;
    std::vector<std::int32_t> graph_haplotype_path_ends;
    std::vector<std::uint32_t> graph_haplotype_path_support;
    std::vector<double> graph_haplotype_path_scores;
    std::vector<std::uint8_t> graph_haplotype_path_has_non_reference_edge;
    std::vector<std::uint32_t> graph_haplotype_path_alt_read_starts;
    std::vector<std::uint32_t> graph_haplotype_path_alt_read_ends;
    // Reference-connected provenance for each materialized graph haplotype.
    // The alignment score is produced by the same SW scoring contract used by
    // assembly; variable-length CIGAR traceback remains Host-side by design.
    std::vector<std::string> graph_haplotype_path_cigars;
    std::vector<std::int32_t> graph_haplotype_path_alignment_scores;
    std::vector<std::int32_t> graph_haplotype_path_alignment_offsets;
    // One EventMap per retained graph haplotype.  This is populated from the
    // Kokkos-scored, Host-traced haplotype-to-reference CIGAR before the
    // candidate union is built, matching EventMap.buildEventMapsForHaplotypes.
    std::vector<HaplotypeEventMap> graph_haplotype_event_maps;
    std::size_t graph_haplotype_event_count = 0;
    std::size_t graph_haplotype_provenance_count = 0;
    std::uint64_t graph_haplotype_cigar_signature = 1469598103934665603ULL;
    // Haplotype-to-reference provenance score telemetry.  The affine score
    // recurrence is executed through the shared Kokkos SW API; only the
    // variable-length traceback remains Host-side.
    bool graph_haplotype_sw_used = false;
    std::string graph_haplotype_sw_execution_space;
    std::size_t graph_haplotype_sw_simd_width = 1;
    std::size_t graph_haplotype_sw_simd_groups = 0;
    double graph_haplotype_sw_prepare_seconds = 0.0;
    double graph_haplotype_sw_seconds = 0.0;
    double graph_prepare_seconds = 0.0;
    double graph_seconds = 0.0;
    std::string graph_execution_space;
    bool activity_used = false;
    bool activity_reference_aware = false;
    bool activity_quality_aware_somatic = false;
    std::size_t activity_loci = 0;
    std::size_t active_loci = 0;
    std::size_t assembly_regions = 0;
    // Raw active/inactive profile segments in traversal order.  The shared
    // calling path uses the padded active windows above; this Host metadata is
    // only for GATK's --assembly-region-out/debug boundary.
    std::vector<fastgatk::kernels::ActivityRegion> activity_profile_regions;
    // Sparse, Kokkos-classified raw active loci in the source traversal. The
    // profile segments above are band-pass/AssemblyRegion state; Mutect2's
    // callable statistic additionally credits raw active sites below its
    // configured depth threshold, so it needs this unexpanded Host trace.
    std::vector<std::pair<std::int32_t, std::int32_t>> activity_active_loci;
    // EventMap-trimmed HC reference-confidence transactions.  `start/end`
    // are the padded read/reference span and `active_start/active_end` are
    // the corresponding variant span whose RCM observations must use that
    // padding.  This keeps AssemblyRegionTrimmer ownership on the Host while
    // the Ref-vs-Any and genotype reductions remain Kokkos kernels.
    std::vector<fastgatk::kernels::ActivityRegion> reference_confidence_regions;
    // Final padded calling windows used for assembly and PairHMM.  This Host
    // planning metadata lets Mutect2 pass its tumor-driven region partition
    // to the matched-normal traversal.
    std::vector<fastgatk::kernels::ActivityRegion> calling_regions;
    // Optional region-level joint assembly snapshots.  These retain the
    // graph/EventMap pairing that is lost by the historical flattened graph
    // provenance vectors above; they are consumed by Mutect2's per-sample
    // likelihood passes and never cross the Kokkos API boundary.
    std::vector<AssemblyRegionAssembly> assembly_region_results;
    // Completed sample-local likelihood results in AssemblyRegion traversal
    // order.  A shared_ptr keeps Result self-referential without exposing
    // large read/Pairs matrices to the Kokkos ABI.
    std::vector<std::shared_ptr<Result>> assembly_region_likelihood_results;
    std::size_t activity_filter_size = 0;
    std::size_t activity_effective_max_probability_propagation_distance = 0;
    double activity_prepare_seconds = 0.0;
    double activity_seconds = 0.0;
    std::string activity_execution_space;
    bool read_filter_used = false;
    bool read_filter_applied = false;
    bool post_read_filter_transform_applied = false;
    std::size_t post_read_filter_transformed_reads = 0;
    std::size_t filtered_reads = 0;
    std::size_t downsampled_reads = 0;
    double read_filter_prepare_seconds = 0.0;
    double read_filter_seconds = 0.0;
    std::string read_filter_execution_space;
    bool overlapping_quality_correction_used = false;
    bool overlapping_quality_correction_metadata_available = false;
    std::size_t overlapping_pairs = 0;
    std::size_t overlapping_bases = 0;
    std::size_t overlapping_conflicting_bases = 0;
    std::size_t overlapping_quality_caps = 0;
    std::size_t overlapping_indel_quality_caps = 0;
    // Effective full Phred settings supplied to Mutect2's overlap cleanup.
    // The Host divides them by two while preparing each overlapping pair,
    // exactly as Mutect2Engine.callRegion does before graph/PairHMM work.
    std::uint8_t overlapping_pcr_snv_quality = 40;
    std::uint8_t overlapping_pcr_indel_quality = 40;
    bool error_correction_used = false;
    std::size_t error_correction_solid_kmers = 0;
    std::size_t error_correction_corrected_kmers = 0;
    std::size_t error_correction_uncorrectable_kmers = 0;
    std::size_t error_correction_pileup_loci = 0;
    std::size_t error_correction_pileup_skipped_indel_adjacent_bases = 0;
    std::size_t error_correction_corrected_reads = 0;
    std::size_t error_correction_corrected_bases = 0;
    double error_correction_prepare_seconds = 0.0;
    double error_correction_seconds = 0.0;
    std::string error_correction_execution_space;
    std::string error_correction_mode = "none";
    // ReferenceConfidenceModel pileup likelihoods are evaluated by a shared
    // Kokkos kernel; indel-informativeness and PL policy remain Host-side.
    double reference_confidence_prepare_seconds = 0.0;
    double reference_confidence_seconds = 0.0;
    std::string reference_confidence_execution_space;
    // Whole-reference GVCF emission uses an indexed coordinate lookup for
    // observed confidence rows and binary-searched candidate spans.  Keep the
    // policy in the manifest so resource/performance audits can distinguish
    // the bounded indexed path from a legacy linear fallback.
    bool reference_block_lookup_indexed = true;
    std::vector<ReferenceBlock> reference_blocks;
    std::vector<SomaticReferenceConfidenceLocus> somatic_reference_confidence_loci;
    double somatic_reference_confidence_prepare_seconds = 0.0;
    double somatic_reference_confidence_seconds = 0.0;
    std::string somatic_reference_confidence_execution_space;
    // Streaming callers can flush a contig immediately after writing it
    // instead of retaining every reference block in the aggregate result.
    // The vector remains the source of truth for the normal in-memory path;
    // this scalar carries its count across bounded stream reductions.
    std::size_t reference_block_count = 0;
    // PairHMM-derived per-read likelihoods for the concrete candidate alleles.
    // Layout is candidate x input-record; a non-finite value means that the
    // read did not cover that candidate.  The Host VCF writer uses this matrix
    // to form true diploid multi-allelic genotype likelihoods instead of
    // combining independently normalized biallelic PLs.
    std::vector<std::vector<double>> allele_read_likelihoods;
    std::vector<std::vector<double>> reference_read_likelihoods;
    // Optional candidate x input-record rows for EventMap's symbolic
    // spanning-deletion (`*`) allele. HC uses this row in the pre-subsetting
    // AF/QUAL posterior and, when GenotypingEngine's upstream-deletion state
    // says that the event is not spurious, in the emitted gVCF allele list.
    // Empty candidate rows mean that no retained haplotype spans that locus,
    // avoiding an extra dense matrix for ordinary calls.
    std::vector<std::vector<double>> spanning_deletion_read_likelihoods;
    std::vector<MultiallelicDepth> multiallelic_depths;
    // Source-record mask after the PairHMM poorly-modelled-read gate.  The
    // annotation interval also contains margin-only reads with no concrete
    // candidate likelihood row; retaining this mask lets MQ include those
    // reads while excluding disqualified evidence consistently.
    std::vector<std::uint8_t> annotation_read_qualified;
    // Read names corresponding to the second dimension above after the
    // shared read-filter/overlap-quality stage.  Mutect2 groups evidence by
    // GATKRead::getName before somatic genotyping; keeping this identity map
    // in Result lets that Host boundary reproduce the grouping without
    // leaking HTSlib objects into the Kokkos kernels.
    std::vector<std::string> likelihood_read_names;
    // Coordinate identity for the same post-filtered records. Mutect2's
    // SomaticGenotypingEngine retains a Fragment when its span overlaps the
    // candidate's informative interval; keeping compact tid/start/end
    // vectors here reproduces that Host-side rule after read filtering
    // without retaining a second ReadBatch or passing HTSlib objects into
    // Kokkos.
    std::vector<std::int32_t> likelihood_read_tids;
    std::vector<std::int64_t> likelihood_read_starts;
    std::vector<std::int64_t> likelihood_read_ends;
    // Base qualities for the same compact likelihood rows after the Host
    // overlap-mate correction.  PairHMM consumes this view, and Mutect2's
    // post-realignment per-allele annotations must use it as well; reading
    // the outer BAM qualities instead loses GATK's callRegion() Q20 cap.
    std::vector<std::uint32_t> likelihood_read_quality_offsets;
    std::vector<std::uint8_t> likelihood_read_qualities;
    // Stable outer-batch record identity for the same likelihood rows.  This
    // is Host-only provenance: AssemblyRegion slicing preserves it so the
    // Mutect2 annotation boundary never guesses an input record from an
    // ambiguous name/coordinate tuple.
    std::vector<std::uint32_t> likelihood_read_source_records;
    // Haplotype realignment metadata for the same likelihood rows. A negative
    // start marks a source record that fell back to its original alignment.
    // Mutect2 uses these intervals for fragment retention before its
    // annotation engine evaluates FORMAT/SB, F1R2 and F2R1.
    std::vector<std::int32_t> likelihood_read_realigned_tids;
    std::vector<std::int64_t> likelihood_read_realigned_starts;
    std::vector<std::int64_t> likelihood_read_realigned_ends;
    // Sparse base projections through the selected haplotype. The native
    // Mutect2 writer uses these only to reproduce GATK's per-read base-quality
    // lookup for OrientationBiasReadCounts after realignment.
    std::vector<std::uint32_t> likelihood_read_realigned_base_records;
    std::vector<std::uint32_t> likelihood_read_realigned_base_offsets;
    std::vector<std::int32_t> likelihood_read_realigned_base_positions;
    // Assembly regions can overlap. GATK performs its read realignment in the
    // candidate's own likelihood context, so a single source read may have a
    // different selected-haplotype projection for adjacent calls. These
    // Sparse AssemblyRegion-local contexts retain that distinction for somatic
    // annotations without allocating candidate x all-input-read matrices.
    // A negative start means the source record has no candidate-scoped
    // projection and consumers must use the original alignment.
    // candidate_read_context_ordinals maps each candidate to its owning
    // AssemblyRegion likelihood collection, allowing all EventMap entries in
    // that collection to share one GATK readRealignment.
    std::vector<std::uint32_t> likelihood_candidate_read_context_ordinals;
    std::vector<CandidateReadRealignment> likelihood_candidate_read_realignments;
    std::vector<std::uint32_t> likelihood_candidate_read_realigned_base_contexts;
    std::vector<std::uint32_t> likelihood_candidate_read_realigned_base_records;
    std::vector<std::uint32_t> likelihood_candidate_read_realigned_base_offsets;
    std::vector<std::int32_t> likelihood_candidate_read_realigned_base_positions;
    // A zero marks the synthetic offset used by ReadPosition when the
    // reference coordinate falls inside a realigned deletion.  It is not a
    // concrete read base and must not feed BaseQuality or F1R2/F2R1.
    std::vector<std::uint8_t> likelihood_candidate_read_realigned_base_is_concrete;
    // PairHMM runs on an AssemblyRegion EventMap, whereas this candidate list
    // also retains raw pileup/CIGAR observations for annotation bookkeeping.
    // Preserve exact EventMap membership at the Host likelihood boundary so
    // later regional somatic reductions cannot add raw sibling alleles.
    std::vector<std::uint8_t> likelihood_candidate_event_map_owned;
    // GATK's <NON_REF> likelihood is derived per read from the median of
    // qualified concrete-allele likelihoods.  This matrix is populated by
    // the reference-backed PairHMM path and consumed by the GVCF writer.
    std::vector<std::vector<double>> non_ref_read_likelihoods;
    // Mutect2 groups read likelihoods by fragment before marginalizing
    // haplotypes to alleles.  Candidate-major CSR preserves that boundary
    // without constructing an infeasible candidate x all-fragment matrix.
    // `somatic_candidate_group_offsets` has candidates.size() + 1 entries;
    // its final value equals `somatic_group_likelihoods.size()`.
    bool somatic_likelihoods_grouped = false;
    // One post-filtered source-record ordinal per fragment group; used only
    // to place a sparse Mutect2 pileup fallback on an overlapping fragment.
    std::vector<std::uint32_t> somatic_group_representative_sources;
    // AssemblyRegion likelihood collection that owns each fragment group.
    // Fragment identity is local to this collection in GATK, even when a
    // source read name recurs in an overlapping region.
    std::vector<std::uint32_t> somatic_group_context_ordinals;
    std::vector<std::uint32_t> somatic_candidate_group_offsets;
    std::vector<SomaticGroupLikelihood> somatic_group_likelihoods;
    // For every concrete candidate, the candidate indices carried by the
    // best read-supported haplotype that contains that candidate.  The
    // PairHMM likelihoods themselves stay in Kokkos; this compact Host-only
    // EventMap projection preserves SomaticGenotypingEngine's ECNTH boundary
    // (best supporting haplotype, then count potential somatic events on it).
    // It is region-local by construction and therefore never represents a
    // synthetic chromosome-wide haplotype.
    std::vector<std::vector<std::uint32_t>> somatic_best_haplotype_candidate_indices;
    // EventMap membership of every retained, concrete alternate haplotype for
    // each candidate. The IDs are local to this Result's AssemblyRegion and
    // let the Host reproduce AssemblyBasedCallerUtils.phaseCalls() without
    // inferring phase from shared read support.
    std::vector<std::vector<std::uint32_t>> somatic_candidate_haplotype_indices;
    // Mutect2's FORMAT/FAD is independently annotated from the retained
    // fragment allele matrix.  AD/DP remain the corresponding read-level
    // BestAllele counts on AssemblyCandidate; keeping FAD separate prevents
    // the VCF writer from conflating the two evidence granularities.
    std::vector<std::uint32_t> somatic_fragment_depth;
    std::vector<std::uint32_t> somatic_fragment_reference_count;
    std::vector<std::uint32_t> somatic_fragment_alternate_count;
    // Same reduction rule as reference_block_count: avoid retaining the
    // per-read matrix after a streamed contig has been materialized.
    bool non_ref_likelihoods_present = false;
    // True only when the explicit USE_POSTERIOR_PROBABILITIES assignment mode
    // supplied GATK GenotypePriorCalculator-style priors to final GT/GQ.
    // Keeping this bit with the result prevents writers from guessing whether
    // zero-valued priors mean flat PLS assignment or an unavailable path.
    bool genotype_priors_used = false;
    // Preserve the configured heterozygosity values with streamed results so
    // gVCF candidate writers can materialize the symbolic <NON_REF> allele's
    // GATK OTHER prior without reaching back into CLI options.
    double genotype_snp_heterozygosity = 1.0e-3;
    double genotype_indel_heterozygosity = 1.0 / 8000.0;
    double genotype_heterozygosity_stdev = 1.0e-2;
    // GenotypingEngine applies this same threshold twice: once to decide
    // whether a biallelic candidate becomes a call, then again to each ALT
    // in a joint multi-allelic VariantContext.  Retain it with streamed
    // results so the VCF writer can perform the latter GATK output-allele
    // subset without guessing from its call list.
    double genotype_standard_confidence_for_calling = 30.0;
    std::vector<double> candidate_prior_het;
    std::vector<double> candidate_prior_hom_alt;
    std::size_t genotype_prior_kernel_calls = 0;
    double genotype_prior_prepare_seconds = 0.0;
    double genotype_prior_seconds = 0.0;
    std::string genotype_prior_execution_space;
    double genotyping_seconds = 0.0;
    double call_confidence_prepare_seconds = 0.0;
    double call_confidence_seconds = 0.0;
    std::string call_confidence_execution_space;
    std::vector<AssemblyCandidate> candidates;
    std::vector<Likelihoods> likelihoods;
    std::vector<GenotypeCall> calls;
    // In reference-confidence mode GATK preserves an EventMap site that
    // reaches GenotypingEngine even when its final concrete GT is hom-ref.
    // These rows are emitted only by the gVCF writer and never enter the
    // ordinary VCF call set; keeping them separate makes that output
    // distinction explicit at the C++ Host boundary.
    std::vector<GenotypeCall> gvcf_hom_ref_calls;
};

// `reference_sequences` is indexed by BAM tid and may be empty.  If a
// reference is unavailable, assembly uses the deterministic majority base at
// each locus as the local reference.  This fallback is useful for smoke tests,
// but it is intentionally not advertised as GATK-compatible variant calling.
Result run(const io::ReadBatch& reads,
           const std::vector<std::string>& reference_sequences,
           const Options& options = {});

// Recreate HC's final VariantContext annotation boundary after joint allele
// subsetting.  `result` supplies the original EventMap likelihood rows, while
// `output_candidates` is the concrete REF/ALT set that will actually be
// serialized.  BestAllele is intentionally chosen from the complete original
// multi-allelic collection, then evidence assigned to an ALT that was removed
// by the joint AF step is ignored by StrandBias/RankSum exactly as in GATK.
// This is Host-side VCF semantics; PairHMM and all Kokkos likelihood/AF work
// remain unchanged.
std::optional<GenotypeCall::Annotations> calculate_output_variant_annotations(
    const io::ReadBatch& reads, const Result& result,
    const std::vector<const AssemblyCandidate*>& output_candidates,
    double qual, std::uint32_t informative_read_overlap_margin = 2,
    bool spanning_deletion_is_output = false,
    bool restrict_likelihoods_to_output_candidates = false);

}  // namespace fastgatk::calling
