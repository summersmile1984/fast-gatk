#include "fastgatk/calling/pipeline.hpp"

#include <Kokkos_Core.hpp>

#include <cstdint>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

// The first `length` bases of a quaternary de Bruijn cycle of order five.
// Every in-range 5-mer is unique, so a focused ReadThreadingGraph fixture
// does not accidentally exercise the non-unique-reference rejection path.
std::string unique_five_mer_reference(const std::size_t length) {
    constexpr std::size_t kOrder = 5;
    constexpr char kBases[] = {'A', 'C', 'G', 'T'};
    std::vector<std::size_t> symbols(kOrder * 2 + 1, 0);
    std::string cycle;
    cycle.reserve(1U << (kOrder * 2));
    std::function<void(std::size_t, std::size_t)> visit =
        [&](const std::size_t t, const std::size_t p) {
            if (t > kOrder) {
                if (kOrder % p == 0)
                    for (std::size_t index = 1; index <= p; ++index)
                        cycle.push_back(kBases[symbols[index]]);
                return;
            }
            symbols[t] = symbols[t - p];
            visit(t + 1, p);
            for (std::size_t symbol = symbols[t - p] + 1; symbol < 4; ++symbol) {
                symbols[t] = symbol;
                visit(t + 1, t);
            }
        };
    visit(1, 1);
    if (length > cycle.size()) throw std::invalid_argument("unique fixture reference is too long");
    return cycle.substr(0, length);
}

}  // namespace

int main() {
    Kokkos::initialize();
    try {
        // LongHomopolymerHaplotypeCollapsingEngine parity: a flow-limited
        // four-base A run is expanded back across the reference deletion by
        // the shared Kokkos-SW/Host-traceback path.  The reverse-direction
        // pass and the optional re-collapse mode are covered as well.
        const std::vector<std::uint8_t> collapsed{
            'C','C','C','C','C','A','A','A','A','G','G','G'};
        const std::vector<std::uint8_t> long_reference{
            'C','C','C','C','C','A','A','A','A','A','A','A','A','A','A','A','A','G','G','G'};
        const auto uncollapsed = fastgatk::calling::uncollapse_flow_hmers(
            collapsed, long_reference, 4, false, false);
        if (!uncollapsed.expanded || uncollapsed.bases != long_reference)
            throw std::runtime_error("flow HMER uncollapse did not restore reference run");
        const auto relimited = fastgatk::calling::uncollapse_flow_hmers(
            collapsed, long_reference, 4, false, true);
        const std::vector<std::uint8_t> relimited_expected{
            'C','C','C','C','C','A','A','A','A','G','G','G'};
        if (relimited.bases != relimited_expected)
            throw std::runtime_error("flow HMER limit-to-threshold mode is not deterministic");
        const std::vector<std::uint8_t> ordinary{'A','C','G','T'};
        const auto ordinary_result = fastgatk::calling::uncollapse_flow_hmers(
            ordinary, ordinary, 4, false, false);
        if (ordinary_result.expanded || ordinary_result.bases != ordinary)
            throw std::runtime_error("ordinary flow HMER region was changed");

        fastgatk::io::ReadBatch reads;
        reads.offsets.push_back(0);
        reads.cigar_offsets.push_back(0);
        reads.insertion_quality_offsets.push_back(0);
        reads.deletion_quality_offsets.push_back(0);
        const auto match = fastgatk::io::CigarOpCode::Match;
        const auto insertion = fastgatk::io::CigarOpCode::Insertion;
        // A non-repeating reference flank is essential here.  With every
        // k-mer equal to AAAAA, the read branch has neither a graph source
        // nor a graph sink and GATK's removePathsNotConnectedToRef() removes
        // it after dangling-end recovery.  Test the explicit CIGAR path with
        // an actual source-to-sink assembly branch instead of assuming that
        // a raw BAM CIGAR independently becomes an HC EventMap allele.
        const std::string cigar_reference =
            "CTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAA";
        const auto cigar_read = cigar_reference.substr(0, 10) + "C" +
            cigar_reference.substr(10, 10);
        for (std::size_t record = 0; record < 4; ++record) {
            reads.bases.insert(reads.bases.end(), cigar_read.begin(), cigar_read.end());
            reads.qualities.insert(reads.qualities.end(), cigar_read.size(), 30);
            // Exercise the GATK BI/BD path, including the Q<6 normalization
            // performed by PairHMMLikelihoodCalculationEngine.
            reads.insertion_qualities.insert(reads.insertion_qualities.end(), cigar_read.size(), 5);
            reads.deletion_qualities.insert(reads.deletion_qualities.end(), cigar_read.size(), 6);
            reads.offsets.push_back(static_cast<std::uint32_t>(reads.bases.size()));
            reads.insertion_quality_offsets.push_back(
                static_cast<std::uint32_t>(reads.insertion_qualities.size()));
            reads.deletion_quality_offsets.push_back(
                static_cast<std::uint32_t>(reads.deletion_qualities.size()));
            reads.positions.push_back(0);
            reads.tids.push_back(0);
            reads.mapq.push_back(60);
            reads.cigar_ops.push_back(fastgatk::io::CigarOp{10, match}.pack());
            reads.cigar_ops.push_back(fastgatk::io::CigarOp{1, insertion}.pack());
            reads.cigar_ops.push_back(fastgatk::io::CigarOp{10, match}.pack());
            reads.cigar_offsets.push_back(static_cast<std::uint32_t>(reads.cigar_ops.size()));
            reads.flags.push_back(0);
        }
        if (!reads.cigar_layout_valid()) throw std::runtime_error("synthetic indel CIGAR invalid");
        fastgatk::calling::Options options;
        options.min_depth = 1;
        options.min_alt_support = 2;
        options.min_base_quality = 10;
        options.max_candidates = 16;
        // These synthetic homopolymer references intentionally exercise the
        // low-level graph path. GATK's production default rejects repeated
        // reference k-mers; opt in explicitly for this focused graph oracle.
        options.graph_allow_non_unique_kmers_in_ref = true;
        const auto result = fastgatk::calling::run(
            reads, {cigar_reference}, options);
        bool found = false;
        for (const auto& candidate : result.candidates) {
            if (candidate.position == 9 && candidate.reference_allele == "T" &&
                candidate.alternate_allele == "TC") {
                found = true;
                break;
            }
        }
        if (!found) throw std::runtime_error("CIGAR insertion candidate was not assembled");

        // ReferenceConfidenceModel does not promote every CIGAR indel
        // anchor to an active site.  With one Q30 insertion against one
        // hundred Q30 reference reads, the hom-ref genotype remains the raw
        // likelihood maximum and GATK's active-region fast path returns
        // zero before applying priors.  This guards the Host/Kokkos activity
        // handoff against restoring the old unconditional-I/D activation
        // shortcut (deletions instead contribute their source Q30 pileup
        // elements at each deleted reference base).
        fastgatk::io::ReadBatch reference_dominant_indel_reads;
        reference_dominant_indel_reads.offsets.push_back(0);
        reference_dominant_indel_reads.cigar_offsets.push_back(0);
        const auto append_reference_dominant_read = [&](const std::string& sequence,
                                                         const bool with_insertion) {
            reference_dominant_indel_reads.bases.insert(
                reference_dominant_indel_reads.bases.end(), sequence.begin(), sequence.end());
            reference_dominant_indel_reads.qualities.insert(
                reference_dominant_indel_reads.qualities.end(), sequence.size(), 30U);
            reference_dominant_indel_reads.offsets.push_back(
                static_cast<std::uint32_t>(reference_dominant_indel_reads.bases.size()));
            reference_dominant_indel_reads.positions.push_back(0);
            reference_dominant_indel_reads.tids.push_back(0);
            reference_dominant_indel_reads.mapq.push_back(60);
            reference_dominant_indel_reads.flags.push_back(0);
            reference_dominant_indel_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{10, match}.pack());
            if (with_insertion)
                reference_dominant_indel_reads.cigar_ops.push_back(
                    fastgatk::io::CigarOp{1, insertion}.pack());
            reference_dominant_indel_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{10, match}.pack());
            reference_dominant_indel_reads.cigar_offsets.push_back(
                static_cast<std::uint32_t>(reference_dominant_indel_reads.cigar_ops.size()));
        };
        for (std::size_t record = 0; record < 100; ++record)
            append_reference_dominant_read(cigar_reference.substr(0, 20), false);
        append_reference_dominant_read(cigar_read, true);
        const auto reference_dominant_indel_result = fastgatk::calling::run(
            reference_dominant_indel_reads, {cigar_reference}, options);
        if (reference_dominant_indel_result.active_loci != 0 ||
            reference_dominant_indel_result.assembly_regions != 0 ||
            !reference_dominant_indel_result.candidates.empty())
            throw std::runtime_error(
                "reference-dominant single-indel activity ignored the GATK posterior gate");

        // A BAM hard clip is not present in SEQ.  Flow read alignments often
        // carry large terminal H operations while retaining a short local
        // segment (eg. 78H41M116H); the assembly region must thread that
        // retained segment instead of treating the hard-clip length as an
        // out-of-bounds read offset.  Keep the same graph-connected allele
        // as above but make both terminal hard clips larger than the stored
        // sequence, which catches exactly that Host CIGAR accounting error.
        fastgatk::io::ReadBatch hard_clipped_reads;
        hard_clipped_reads.offsets.push_back(0);
        hard_clipped_reads.cigar_offsets.push_back(0);
        for (std::size_t record = 0; record < 4; ++record) {
            hard_clipped_reads.bases.insert(
                hard_clipped_reads.bases.end(), cigar_read.begin(), cigar_read.end());
            hard_clipped_reads.qualities.insert(
                hard_clipped_reads.qualities.end(), cigar_read.size(), 30);
            hard_clipped_reads.offsets.push_back(
                static_cast<std::uint32_t>(hard_clipped_reads.bases.size()));
            hard_clipped_reads.positions.push_back(0);
            hard_clipped_reads.tids.push_back(0);
            hard_clipped_reads.mapq.push_back(60);
            hard_clipped_reads.cigar_ops.push_back(fastgatk::io::CigarOp{78,
                fastgatk::io::CigarOpCode::HardClip}.pack());
            hard_clipped_reads.cigar_ops.push_back(fastgatk::io::CigarOp{10, match}.pack());
            hard_clipped_reads.cigar_ops.push_back(fastgatk::io::CigarOp{1, insertion}.pack());
            hard_clipped_reads.cigar_ops.push_back(fastgatk::io::CigarOp{10, match}.pack());
            hard_clipped_reads.cigar_ops.push_back(fastgatk::io::CigarOp{116,
                fastgatk::io::CigarOpCode::HardClip}.pack());
            hard_clipped_reads.cigar_offsets.push_back(
                static_cast<std::uint32_t>(hard_clipped_reads.cigar_ops.size()));
            hard_clipped_reads.flags.push_back(0);
        }
        const auto hard_clipped_result = fastgatk::calling::run(
            hard_clipped_reads, {cigar_reference}, options);
        bool hard_clipped_found = false;
        for (const auto& candidate : hard_clipped_result.candidates) {
            hard_clipped_found = hard_clipped_found ||
                (candidate.position == 9 && candidate.reference_allele == "T" &&
                 candidate.alternate_allele == "TC" && candidate.graph_derived);
        }
        if (!hard_clipped_found || hard_clipped_result.graph_haplotype_paths < 2)
            throw std::runtime_error("hard-clipped assembly read was not threaded into graph");

        // A joint Mutect2 traversal owns only the immutable AssemblyResult
        // contract.  It must retain the same graph/EventMap candidate before
        // returning, but it must not launch a joint PairHMM matrix whose
        // likelihoods have no consumer; tumor and normal replays score that
        // snapshot separately.  Exercise the internal scheduling mode here
        // so a future refactor cannot accidentally restore the redundant
        // numerical pass while keeping the final VCF tests green.
        auto assembly_only_options = options;
        assembly_only_options.retain_assembly_region_results = true;
        assembly_only_options.assembly_only = true;
        const auto assembly_only_result = fastgatk::calling::run(
            reads, {cigar_reference}, assembly_only_options);
        if (assembly_only_result.pairhmm_used || assembly_only_result.pairhmm_pairs != 0 ||
            assembly_only_result.pairhmm_skip_reason != "assembly-only")
            throw std::runtime_error("AssemblyResult-only traversal executed PairHMM");
        if (assembly_only_result.assembly_region_results.empty())
            throw std::runtime_error("AssemblyResult-only traversal lost its region snapshot");
        bool snapshot_insertion_found = false;
        for (const auto& assembly : assembly_only_result.assembly_region_results) {
            for (const auto& event : assembly.event_map) {
                if (event.position == 9 && event.reference == "T" && event.alternate == "TC") {
                    snapshot_insertion_found = true;
                    break;
                }
            }
            if (snapshot_insertion_found) break;
        }
        if (!snapshot_insertion_found)
            throw std::runtime_error("AssemblyResult-only traversal lost the EventMap insertion");

        // HTSlib exposes BAM SEQ in alignment orientation.  A reverse-strand
        // read carrying reference-oriented insertion AG therefore also has
        // AG in its CIGAR-I slice.  The assembly candidate is anchored on the
        // forward reference and must emit TAG without a second complement.
        fastgatk::io::ReadBatch reverse_cigar_reads;
        reverse_cigar_reads.offsets.push_back(0);
        reverse_cigar_reads.cigar_offsets.push_back(0);
        const std::string reverse_cigar_reference =
            "TACATTTGCTTCGTTGACTAGCAACCCAGGGCTATAGCTA";
        const std::string reverse_cigar_sequence = "TACATTAGTGCTTC";
        for (std::size_t record = 0; record < 4; ++record) {
            reverse_cigar_reads.bases.insert(reverse_cigar_reads.bases.end(),
                                             reverse_cigar_sequence.begin(),
                                             reverse_cigar_sequence.end());
            reverse_cigar_reads.qualities.insert(
                reverse_cigar_reads.qualities.end(), reverse_cigar_sequence.size(), 30);
            reverse_cigar_reads.offsets.push_back(
                static_cast<std::uint32_t>(reverse_cigar_reads.bases.size()));
            reverse_cigar_reads.positions.push_back(0);
            reverse_cigar_reads.tids.push_back(0);
            reverse_cigar_reads.mapq.push_back(60);
            reverse_cigar_reads.flags.push_back(0x10U);
            reverse_cigar_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{6, match}.pack());
            reverse_cigar_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{2, insertion}.pack());
            reverse_cigar_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{6, match}.pack());
            reverse_cigar_reads.cigar_offsets.push_back(
                static_cast<std::uint32_t>(reverse_cigar_reads.cigar_ops.size()));
        }
        if (!reverse_cigar_reads.cigar_layout_valid())
            throw std::runtime_error("synthetic reverse CIGAR invalid");
        const auto reverse_cigar_result = fastgatk::calling::run(
            reverse_cigar_reads, {reverse_cigar_reference}, options);
        bool reverse_insertion_found = false;
        std::string reverse_insertion_alt;
        for (const auto& candidate : reverse_cigar_result.candidates) {
            if (candidate.position == 5 && candidate.reference_allele == "T" &&
                candidate.alternate_allele == "TAG") {
                reverse_insertion_found = true;
                reverse_insertion_alt = candidate.alternate_allele;
                break;
            }
        }
        if (!reverse_insertion_found) {
            std::string detail;
            for (const auto& candidate : reverse_cigar_result.candidates) {
                if (!detail.empty()) detail += ";";
                detail += std::to_string(candidate.position) + ":" +
                    candidate.reference_allele + ">" + candidate.alternate_allele;
            }
            throw std::runtime_error("reverse-strand CIGAR insertion was not reference-oriented: " + detail);
        }

        // A mixed-strand CIGAR event must retain one common, non-repeat
        // representation.  In particular, G>GA must retain its anchored
        // coordinate rather than being normalized into a different haplotype.
        // This mirrors the DREAM chr20 10080550 regression, where a previous
        // bad normalization split three G>GA reads into a non-emittable
        // singleton.
        fastgatk::io::ReadBatch mixed_strand_insertion_reads;
        mixed_strand_insertion_reads.offsets.push_back(0);
        mixed_strand_insertion_reads.cigar_offsets.push_back(0);
        const std::string mixed_strand_reference =
            "AACTGGCGAGTGGAGGACACATTAATAATTTGCTCACTCC";
        const std::string mixed_strand_sequence = "AACTGGACGAGTG";
        for (std::size_t record = 0; record < 3; ++record) {
            mixed_strand_insertion_reads.bases.insert(
                mixed_strand_insertion_reads.bases.end(), mixed_strand_sequence.begin(),
                mixed_strand_sequence.end());
            mixed_strand_insertion_reads.qualities.insert(
                mixed_strand_insertion_reads.qualities.end(), mixed_strand_sequence.size(), 30);
            mixed_strand_insertion_reads.offsets.push_back(
                static_cast<std::uint32_t>(mixed_strand_insertion_reads.bases.size()));
            mixed_strand_insertion_reads.positions.push_back(0);
            mixed_strand_insertion_reads.tids.push_back(0);
            mixed_strand_insertion_reads.mapq.push_back(60);
            mixed_strand_insertion_reads.flags.push_back(record == 2 ? 0x10U : 0U);
            mixed_strand_insertion_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{6, match}.pack());
            mixed_strand_insertion_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{1, insertion}.pack());
            mixed_strand_insertion_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{6, match}.pack());
            mixed_strand_insertion_reads.cigar_offsets.push_back(
                static_cast<std::uint32_t>(mixed_strand_insertion_reads.cigar_ops.size()));
        }
        const auto mixed_strand_insertion_result = fastgatk::calling::run(
            mixed_strand_insertion_reads, {mixed_strand_reference}, options);
        bool mixed_strand_insertion_found = false;
        for (const auto& candidate : mixed_strand_insertion_result.candidates) {
            if (candidate.position == 5 && candidate.reference_allele == "G" &&
                candidate.alternate_allele == "GA" && candidate.alternate_count == 3) {
                mixed_strand_insertion_found = true;
                break;
            }
        }
        if (!mixed_strand_insertion_found)
            throw std::runtime_error("mixed-strand G>GA insertion was split or left-shifted");
        if (result.pairhmm_insertion_quality_reads != 4 ||
            result.pairhmm_deletion_quality_reads != 4)
            throw std::runtime_error("PairHMM did not consume BI/BD per-base qualities");

        // GATK-shaped PairHMM Host controls must be accepted by the shared
        // calling API (HC and Mutect2 both route through this same Options).
        auto pairhmm_options = options;
        pairhmm_options.pairhmm_base_quality_score_threshold = 31;
        pairhmm_options.pairhmm_disable_cap_base_qualities_to_mapq = true;
        const auto pairhmm_controls = fastgatk::calling::run(
            reads, {cigar_reference}, pairhmm_options);
        if (pairhmm_controls.pairhmm_insertion_quality_reads != 4 ||
            pairhmm_controls.pairhmm_deletion_quality_reads != 4)
            throw std::runtime_error("PairHMM quality controls bypassed BI/BD consumption");
        bool indel_prior_checked = false;
        for (std::size_t index = 0; index < result.candidates.size(); ++index) {
            const auto& candidate = result.candidates[index];
            if (!candidate.reference_allele.empty() &&
                candidate.reference_allele.size() != candidate.alternate_allele.size() &&
                index < result.candidate_prior_het.size() && index < result.candidate_prior_hom_alt.size()) {
                // GATK HomoSapiensConstants.INDEL_HETEROZYGOSITY is 1/8000,
                // not the historical 1e-4 approximation.  Keep this derived
                // from the source value so the smoke oracle verifies the same
                // GenotypePriorCalculator.assumingHW contract as the caller.
                const double expected_indel_het = std::log10(1.0 / 8000.0);
                if (std::abs(result.candidate_prior_het[index] - expected_indel_het) > 1e-9 ||
                    std::abs(result.candidate_prior_hom_alt[index] - 2.0 * expected_indel_het) > 1e-9)
                    throw std::runtime_error("INDEL assumingHW prior is not GATK-compatible");
                indel_prior_checked = true;
                break;
            }
        }
        if (!indel_prior_checked) throw std::runtime_error("INDEL prior regression candidate missing");

        // The caller must consume the same deterministic Host downsampling
        // mask used by the kernel filter, and a repeated run with the same
        // seed must retain the same ordinals/signature.
        auto downsample_options = options;
        downsample_options.max_reads_per_locus = 2;
        downsample_options.downsampling_seed = 0x123456789abcdef0ULL;
        const auto downsample_a = fastgatk::calling::run(
            reads, {cigar_reference}, downsample_options);
        const auto downsample_b = fastgatk::calling::run(
            reads, {cigar_reference}, downsample_options);
        if (downsample_a.downsampled_reads != 2 ||
            downsample_b.downsampled_reads != 2 ||
            downsample_a.reads != downsample_b.reads ||
            downsample_a.signature != downsample_b.signature)
            throw std::runtime_error("deterministic per-locus downsampling regression failed");

        fastgatk::io::ReadBatch deletion_reads;
        deletion_reads.offsets.push_back(0);
        deletion_reads.cigar_offsets.push_back(0);
        for (std::size_t record = 0; record < 4; ++record) {
            const auto sequence = cigar_reference.substr(0, 10) +
                cigar_reference.substr(11, 10);
            deletion_reads.bases.insert(deletion_reads.bases.end(), sequence.begin(), sequence.end());
            deletion_reads.qualities.insert(deletion_reads.qualities.end(), sequence.size(), 30);
            deletion_reads.offsets.push_back(static_cast<std::uint32_t>(deletion_reads.bases.size()));
            deletion_reads.positions.push_back(0);
            deletion_reads.tids.push_back(0);
            deletion_reads.mapq.push_back(60);
            deletion_reads.cigar_ops.push_back(fastgatk::io::CigarOp{10, match}.pack());
            deletion_reads.cigar_ops.push_back(fastgatk::io::CigarOp{1,
                fastgatk::io::CigarOpCode::Deletion}.pack());
            deletion_reads.cigar_ops.push_back(fastgatk::io::CigarOp{10, match}.pack());
            deletion_reads.cigar_offsets.push_back(static_cast<std::uint32_t>(deletion_reads.cigar_ops.size()));
            deletion_reads.flags.push_back(0);
        }
        const auto deletion_result = fastgatk::calling::run(
            deletion_reads, {cigar_reference}, options);
        found = false;
        for (const auto& candidate : deletion_result.candidates) {
            if (candidate.position == 9 &&
                candidate.reference_allele == cigar_reference.substr(9, 2) &&
                candidate.alternate_allele == cigar_reference.substr(9, 1)) {
                found = true;
                break;
            }
        }
        if (!found) throw std::runtime_error("CIGAR deletion candidate was not assembled");

        // Mixed CIGAR evidence must be represented from the retained graph
        // haplotype, rather than leaking its primitive input I/D operations
        // past EventMap. SequenceMatch/SequenceMismatch are legal BAM
        // operators and must advance both cursors like M.
        fastgatk::io::ReadBatch mixed_reads;
        mixed_reads.offsets.push_back(0);
        mixed_reads.cigar_offsets.push_back(0);
        const std::string mixed_reference = cigar_reference;
        const std::string mixed_read = mixed_reference.substr(0, 5) + "C" + "AA" +
            mixed_reference.substr(9, 5);
        for (std::size_t record = 0; record < 4; ++record) {
            mixed_reads.bases.insert(mixed_reads.bases.end(), mixed_read.begin(), mixed_read.end());
            mixed_reads.qualities.insert(mixed_reads.qualities.end(), mixed_read.size(), 30);
            mixed_reads.offsets.push_back(static_cast<std::uint32_t>(mixed_reads.bases.size()));
            mixed_reads.positions.push_back(0);
            mixed_reads.tids.push_back(0);
            mixed_reads.mapq.push_back(60);
            mixed_reads.cigar_ops.push_back(fastgatk::io::CigarOp{2,
                fastgatk::io::CigarOpCode::SequenceMatch}.pack());
            mixed_reads.cigar_ops.push_back(fastgatk::io::CigarOp{3,
                fastgatk::io::CigarOpCode::Match}.pack());
            mixed_reads.cigar_ops.push_back(fastgatk::io::CigarOp{1,
                fastgatk::io::CigarOpCode::Insertion}.pack());
            mixed_reads.cigar_ops.push_back(fastgatk::io::CigarOp{2,
                fastgatk::io::CigarOpCode::SequenceMismatch}.pack());
            mixed_reads.cigar_ops.push_back(fastgatk::io::CigarOp{2,
                fastgatk::io::CigarOpCode::Deletion}.pack());
            mixed_reads.cigar_ops.push_back(fastgatk::io::CigarOp{5,
                fastgatk::io::CigarOpCode::Match}.pack());
            mixed_reads.cigar_offsets.push_back(
                static_cast<std::uint32_t>(mixed_reads.cigar_ops.size()));
            mixed_reads.flags.push_back(0);
        }
        const auto mixed_result = fastgatk::calling::run(
            mixed_reads, {mixed_reference}, options);
        bool mixed_insertion = false;
        bool mixed_deletion = false;
        for (const auto& candidate : mixed_result.candidates) {
            mixed_insertion = mixed_insertion ||
                (candidate.position == 8 && candidate.reference_allele == "G" &&
                 candidate.alternate_allele == "GCAA");
            mixed_deletion = mixed_deletion ||
                (candidate.position == 4 && candidate.reference_allele == "CTGTG" &&
                 candidate.alternate_allele == "C");
        }
        if (!mixed_insertion || !mixed_deletion)
            throw std::runtime_error("mixed CIGAR edits were not normalized and retained");

        // Two nearby indels on the same read must remain in the retained
        // source-to-sink graph haplotype. EventMap may emit their canonical
        // primitive alleles instead of preserving the input CIGAR block.
        fastgatk::io::ReadBatch compound_reads;
        compound_reads.offsets.push_back(0);
        compound_reads.cigar_offsets.push_back(0);
        const std::string compound_reference = cigar_reference;
        for (std::size_t record = 0; record < 4; ++record) {
            const auto sequence = compound_reference.substr(0, 10) + "GG" +
                compound_reference.substr(10, 3) + compound_reference.substr(14, 10);
            compound_reads.bases.insert(compound_reads.bases.end(), sequence.begin(), sequence.end());
            compound_reads.qualities.insert(compound_reads.qualities.end(), sequence.size(), 30);
            compound_reads.offsets.push_back(static_cast<std::uint32_t>(compound_reads.bases.size()));
            compound_reads.positions.push_back(0);
            compound_reads.tids.push_back(0);
            compound_reads.mapq.push_back(60);
            compound_reads.cigar_ops.push_back(fastgatk::io::CigarOp{10, match}.pack());
            compound_reads.cigar_ops.push_back(fastgatk::io::CigarOp{2, insertion}.pack());
            compound_reads.cigar_ops.push_back(fastgatk::io::CigarOp{3, match}.pack());
            compound_reads.cigar_ops.push_back(fastgatk::io::CigarOp{1,
                fastgatk::io::CigarOpCode::Deletion}.pack());
            compound_reads.cigar_ops.push_back(fastgatk::io::CigarOp{10, match}.pack());
            compound_reads.cigar_offsets.push_back(
                static_cast<std::uint32_t>(compound_reads.cigar_ops.size()));
            compound_reads.flags.push_back(0);
        }
        const auto compound_result = fastgatk::calling::run(
            compound_reads, {compound_reference}, options);
        bool compound_insertion = false;
        bool compound_deletion = false;
        for (const auto& candidate : compound_result.candidates) {
            compound_insertion = compound_insertion ||
                (candidate.position == 9 && candidate.reference_allele == "T" &&
                 candidate.alternate_allele == "TGG");
            compound_deletion = compound_deletion ||
                (candidate.position == 12 && candidate.reference_allele == "AC" &&
                 candidate.alternate_allele == "A");
        }
        if (!compound_insertion || !compound_deletion)
            throw std::runtime_error("compound mixed-CIGAR haplotype was not assembled");

        // Internal soft-clip plus insertion is retained as one compound
        // allele when both events occur on the same short haplotype.  The
        // clipped sequence is not emitted as a speculative standalone event.
        fastgatk::io::ReadBatch compound_softclip_reads;
        compound_softclip_reads.offsets.push_back(0);
        compound_softclip_reads.cigar_offsets.push_back(0);
        const std::string compound_softclip_reference = cigar_reference;
        for (std::size_t record = 0; record < 4; ++record) {
            const std::string sequence = compound_softclip_reference.substr(0, 10) +
                "GG" + "TC" + compound_softclip_reference.substr(10, 10);
            compound_softclip_reads.bases.insert(compound_softclip_reads.bases.end(),
                                                  sequence.begin(), sequence.end());
            compound_softclip_reads.qualities.insert(
                compound_softclip_reads.qualities.end(), sequence.size(), 30);
            compound_softclip_reads.offsets.push_back(
                static_cast<std::uint32_t>(compound_softclip_reads.bases.size()));
            compound_softclip_reads.positions.push_back(0);
            compound_softclip_reads.tids.push_back(0);
            compound_softclip_reads.mapq.push_back(60);
            compound_softclip_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{10, match}.pack());
            compound_softclip_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{2, insertion}.pack());
            compound_softclip_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{2, fastgatk::io::CigarOpCode::SoftClip}.pack());
            compound_softclip_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{10, match}.pack());
            compound_softclip_reads.cigar_offsets.push_back(
                static_cast<std::uint32_t>(compound_softclip_reads.cigar_ops.size()));
            compound_softclip_reads.flags.push_back(0);
        }
        const auto compound_softclip_result = fastgatk::calling::run(
            compound_softclip_reads, {compound_softclip_reference}, options);
        bool compound_softclip_found = false;
        for (const auto& candidate : compound_softclip_result.candidates) {
            compound_softclip_found = compound_softclip_found ||
                (candidate.position == 9 && candidate.reference_allele == "T" &&
                 candidate.alternate_allele == "TGGTC");
        }
        if (!compound_softclip_found || compound_softclip_result.compound_cigar_candidates == 0)
            throw std::runtime_error("compound soft-clip haplotype was not assembled");

        // GATK hard-clips terminal soft-clipped bases before it threads reads
        // into ReadThreadingAssembler.  They can contribute activity, but a
        // terminal clip alone is not an HC EventMap insertion.
        fastgatk::io::ReadBatch softclip_reads;
        softclip_reads.offsets.push_back(0);
        softclip_reads.cigar_offsets.push_back(0);
        for (std::size_t record = 0; record < 4; ++record) {
            const std::string sequence = "CG" + std::string(20, 'A');
            softclip_reads.bases.insert(softclip_reads.bases.end(), sequence.begin(), sequence.end());
            softclip_reads.qualities.insert(softclip_reads.qualities.end(), sequence.size(), 30);
            softclip_reads.offsets.push_back(static_cast<std::uint32_t>(softclip_reads.bases.size()));
            softclip_reads.positions.push_back(1);
            softclip_reads.tids.push_back(0);
            softclip_reads.mapq.push_back(60);
            softclip_reads.cigar_ops.push_back(fastgatk::io::CigarOp{2,
                fastgatk::io::CigarOpCode::SoftClip}.pack());
            softclip_reads.cigar_ops.push_back(fastgatk::io::CigarOp{20, match}.pack());
            softclip_reads.cigar_offsets.push_back(
                static_cast<std::uint32_t>(softclip_reads.cigar_ops.size()));
            softclip_reads.flags.push_back(0);
        }
        const auto softclip_result = fastgatk::calling::run(
            softclip_reads, {std::string(40, 'A')}, options);
        if (!softclip_result.candidates.empty())
            throw std::runtime_error("terminal soft clip escaped GATK EventMap assembly");

        fastgatk::io::ReadBatch reverse_softclip_reads;
        reverse_softclip_reads.offsets.push_back(0);
        reverse_softclip_reads.cigar_offsets.push_back(0);
        for (std::size_t record = 0; record < 4; ++record) {
            // The stored clipped sequence is already in alignment orientation;
            // reverse strand changes the insertion anchor, not its bases.
            const std::string sequence = std::string(20, 'A') + "GT";
            reverse_softclip_reads.bases.insert(
                reverse_softclip_reads.bases.end(), sequence.begin(), sequence.end());
            reverse_softclip_reads.qualities.insert(
                reverse_softclip_reads.qualities.end(), sequence.size(), 30);
            reverse_softclip_reads.offsets.push_back(
                static_cast<std::uint32_t>(reverse_softclip_reads.bases.size()));
            reverse_softclip_reads.positions.push_back(1);
            reverse_softclip_reads.tids.push_back(0);
            reverse_softclip_reads.mapq.push_back(60);
            reverse_softclip_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{20, match}.pack());
            reverse_softclip_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{2, fastgatk::io::CigarOpCode::SoftClip}.pack());
            reverse_softclip_reads.cigar_offsets.push_back(
                static_cast<std::uint32_t>(reverse_softclip_reads.cigar_ops.size()));
            reverse_softclip_reads.flags.push_back(0x10U);
        }
        const auto reverse_softclip_result = fastgatk::calling::run(
            reverse_softclip_reads, {std::string(40, 'A')}, options);
        if (!reverse_softclip_result.candidates.empty())
            throw std::runtime_error("reverse terminal soft clip escaped GATK EventMap assembly");

        // An internal soft clip has no direct CIGAR I operation.  The read
        // sequence still enters the bounded k-mer graph, whose materialized
        // haplotype is aligned back to the reference and promoted to the
        // same left-anchored insertion candidate.
        fastgatk::io::ReadBatch internal_softclip_reads;
        internal_softclip_reads.offsets.push_back(0);
        internal_softclip_reads.cigar_offsets.push_back(0);
        // Use a unique reference context.  An all-A fixture collapses its
        // source and sink into one cyclic k-mer vertex, which has no GATK
        // assembly anchor for a soft-clipped branch and must be discarded.
        const std::string internal_softclip_reference =
            "GCTACGATCAGTCCGATAGCTTGACCTAGTGCATGACGA";
        for (std::size_t record = 0; record < 4; ++record) {
            const std::string sequence = internal_softclip_reference.substr(0, 10) + "CG" +
                internal_softclip_reference.substr(10, 10);
            internal_softclip_reads.bases.insert(internal_softclip_reads.bases.end(),
                                                  sequence.begin(), sequence.end());
            internal_softclip_reads.qualities.insert(
                internal_softclip_reads.qualities.end(), sequence.size(), 30);
            internal_softclip_reads.offsets.push_back(
                static_cast<std::uint32_t>(internal_softclip_reads.bases.size()));
            internal_softclip_reads.positions.push_back(0);
            internal_softclip_reads.tids.push_back(0);
            internal_softclip_reads.mapq.push_back(60);
            internal_softclip_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{10, match}.pack());
            internal_softclip_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{2, fastgatk::io::CigarOpCode::SoftClip}.pack());
            internal_softclip_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{10, match}.pack());
            internal_softclip_reads.cigar_offsets.push_back(
                static_cast<std::uint32_t>(internal_softclip_reads.cigar_ops.size()));
            internal_softclip_reads.flags.push_back(0);
        }
        const auto internal_softclip_result = fastgatk::calling::run(
            internal_softclip_reads, {internal_softclip_reference}, options);
        found = false;
        for (const auto& candidate : internal_softclip_result.candidates) {
            if (candidate.reference_allele == "A" &&
                candidate.alternate_allele == "ACG") {
                found = true;
                break;
            }
        }
        if (!found || internal_softclip_result.graph_variant_candidates == 0)
            throw std::runtime_error("internal soft-clip graph candidate was not assembled");

        // A CIGAR reference skip (N) is not an EventMap deletion.  Production
        // HaplotypeCaller rejects it with WellformedReadFilter; the permissive
        // direct-library entry point must likewise avoid promoting it.
        fastgatk::io::ReadBatch graph_deletion_reads;
        graph_deletion_reads.offsets.push_back(0);
        graph_deletion_reads.cigar_offsets.push_back(0);
        const std::string graph_reference = cigar_reference;
        for (std::size_t record = 0; record < 4; ++record) {
            const auto sequence = graph_reference.substr(0, 20) + graph_reference.substr(22, 20);
            graph_deletion_reads.bases.insert(graph_deletion_reads.bases.end(),
                                              sequence.begin(), sequence.end());
            graph_deletion_reads.qualities.insert(graph_deletion_reads.qualities.end(),
                                                   sequence.size(), 30);
            graph_deletion_reads.offsets.push_back(
                static_cast<std::uint32_t>(graph_deletion_reads.bases.size()));
            graph_deletion_reads.positions.push_back(0);
            graph_deletion_reads.tids.push_back(0);
            graph_deletion_reads.mapq.push_back(60);
            graph_deletion_reads.cigar_ops.push_back(fastgatk::io::CigarOp{20,
                fastgatk::io::CigarOpCode::Match}.pack());
            graph_deletion_reads.cigar_ops.push_back(fastgatk::io::CigarOp{2,
                fastgatk::io::CigarOpCode::ReferenceSkip}.pack());
            graph_deletion_reads.cigar_ops.push_back(fastgatk::io::CigarOp{20,
                fastgatk::io::CigarOpCode::Match}.pack());
            graph_deletion_reads.cigar_offsets.push_back(
                static_cast<std::uint32_t>(graph_deletion_reads.cigar_ops.size()));
            graph_deletion_reads.flags.push_back(0);
        }
        const auto graph_deletion_result = fastgatk::calling::run(
            graph_deletion_reads, {graph_reference}, options);
        if (!graph_deletion_result.candidates.empty())
            throw std::runtime_error("reference skip was promoted to an HC EventMap deletion");

        // A single locus can carry multiple alternate bases.  Keep both
        // candidates through the likelihood path so the VCF writer can emit a
        // true multiallelic record instead of silently dropping the weaker ALT.
        fastgatk::io::ReadBatch multiallelic_reads;
        multiallelic_reads.offsets.push_back(0);
        const std::string multiallelic_reference =
            "TACGCTCCAAGATGCAGCAAACTCACCTAACGTATGTTGCGACTATCTGGCCAGTGTTAAGTCCTAGCAGAGTTTCAATT";
        for (std::size_t record = 0; record < 8; ++record) {
            const char alternate = record < 4 ? 'C' : 'G';
            const auto sequence = multiallelic_reference.substr(0, 20) + alternate +
                multiallelic_reference.substr(21, 20);
            multiallelic_reads.bases.insert(multiallelic_reads.bases.end(), sequence.begin(), sequence.end());
            multiallelic_reads.qualities.insert(multiallelic_reads.qualities.end(), sequence.size(), 30);
            multiallelic_reads.offsets.push_back(static_cast<std::uint32_t>(multiallelic_reads.bases.size()));
            multiallelic_reads.positions.push_back(0);
            multiallelic_reads.tids.push_back(0);
            multiallelic_reads.mapq.push_back(60);
            multiallelic_reads.flags.push_back(0);
        }
        const auto multiallelic_result = fastgatk::calling::run(
            multiallelic_reads, {multiallelic_reference}, options);
        std::size_t multiallelic_candidates = 0;
        for (const auto& candidate : multiallelic_result.candidates)
            if (candidate.position == 20 && candidate.reference == 'A' &&
                (candidate.alternate == 'C' || candidate.alternate == 'G'))
                ++multiallelic_candidates;
        if (multiallelic_candidates != 2)
            throw std::runtime_error("multiple alternate alleles at one locus were not preserved");
        bool snp_prior_checked = false;
        // GenotypePriorCalculator.assumingHW receives the per-ALT
        // heterozygosity in log10 space; it does not divide a multiallelic
        // SNP by three at this layer.
        const double expected_snp_het = std::log10(1.0e-3) - std::log10(3.0);
        for (std::size_t index = 0; index < multiallelic_result.candidates.size(); ++index) {
            const auto& candidate = multiallelic_result.candidates[index];
            if (candidate.reference_allele.empty() && candidate.alternate_allele.empty() &&
                index < multiallelic_result.candidate_prior_het.size()) {
                if (std::abs(multiallelic_result.candidate_prior_het[index] - expected_snp_het) > 1e-9)
                    throw std::runtime_error("SNP assumingHW prior is not GATK-compatible");
                snp_prior_checked = true;
                break;
            }
        }
        if (!snp_prior_checked) throw std::runtime_error("SNP prior regression candidate missing");

        // ReferenceConfidenceModel skips bases whose quality is <= HC's
        // min-base-quality, while ReadThreadingAssembler keeps Q10 bases
        // (baseIsUsableForAssembly uses >=).  Exercise both source behaviors:
        // without another active signal this Q10-only locus is inactive, but
        // --force-active must still thread it and retain EventMap provenance.
        fastgatk::io::ReadBatch graph_snp_reads;
        graph_snp_reads.offsets.push_back(0);
        graph_snp_reads.cigar_offsets.push_back(0);
        const std::string graph_snp_reference = multiallelic_reference;
        for (std::size_t record = 0; record < 4; ++record) {
            auto sequence = graph_snp_reference.substr(0, 41);
            sequence[20] = 'C';
            graph_snp_reads.bases.insert(graph_snp_reads.bases.end(), sequence.begin(), sequence.end());
            graph_snp_reads.qualities.insert(graph_snp_reads.qualities.end(), sequence.size(), 30);
            graph_snp_reads.qualities[graph_snp_reads.qualities.size() - 21] = 10;
            graph_snp_reads.offsets.push_back(static_cast<std::uint32_t>(graph_snp_reads.bases.size()));
            graph_snp_reads.positions.push_back(0);
            graph_snp_reads.tids.push_back(0);
            graph_snp_reads.mapq.push_back(60);
            graph_snp_reads.flags.push_back(0);
            graph_snp_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{41, match}.pack());
            graph_snp_reads.cigar_offsets.push_back(
                static_cast<std::uint32_t>(graph_snp_reads.cigar_ops.size()));
        }
        auto graph_snp_options = options;
        graph_snp_options.min_base_quality = 10;
        const auto inactive_q10_snp_result = fastgatk::calling::run(
            graph_snp_reads, {graph_snp_reference}, graph_snp_options);
        if (inactive_q10_snp_result.assembly_regions != 0 ||
            !inactive_q10_snp_result.candidates.empty())
            throw std::runtime_error("Q10-only HC activity was not filtered before assembly");
        graph_snp_options.force_active = true;
        const auto graph_snp_result = fastgatk::calling::run(
            graph_snp_reads, {graph_snp_reference}, graph_snp_options);
        bool graph_snp_found = false;
        bool graph_snp_is_graph_derived = false;
        for (const auto& candidate : graph_snp_result.candidates) {
            const bool is_graph_snp = candidate.position == 20 && candidate.reference == 'A' &&
                candidate.alternate == 'C';
            graph_snp_found = graph_snp_found || is_graph_snp;
            graph_snp_is_graph_derived = graph_snp_is_graph_derived ||
                (is_graph_snp && candidate.graph_derived);
        }
        // `graph_snp_candidates` counts only candidates newly appended by
        // EventMap.  A CIGAR-provisional candidate with the same allele is
        // deliberately merged and marked graph-derived instead, so asserting
        // that counter here would require a duplicate candidate.
        if (!graph_snp_found || !graph_snp_is_graph_derived ||
            graph_snp_result.pairhmm_graph_snp_posterior_pairs == 0)
            throw std::runtime_error("graph-supported SNP mismatch was not materialized: paths=" +
                                     std::to_string(graph_snp_result.graph_haplotype_paths) +
                                     " candidates=" + std::to_string(graph_snp_result.candidate_sites) +
                                     " graph_snp=" + std::to_string(graph_snp_result.graph_snp_candidates) +
                                     " graph_snp_posterior=" + std::to_string(
                                         graph_snp_result.pairhmm_graph_snp_posterior_pairs) +
                                     " graph_variant=" + std::to_string(graph_snp_result.graph_variant_candidates) +
                                     " pairhmm_graph=" + std::to_string(graph_snp_result.pairhmm_graph_haplotypes) +
                                     " pairhmm_pairs=" + std::to_string(graph_snp_result.pairhmm_pairs) +
                                     " seqs=" + std::to_string(graph_snp_result.graph_haplotype_sequences.size()) +
                                     " nonref=" + std::to_string(graph_snp_result.graph_haplotype_path_has_non_reference_edge.size()) +
                                     " support=" + (graph_snp_result.graph_haplotype_path_support.empty() ?
                                         std::string("none") : std::to_string(graph_snp_result.graph_haplotype_path_support.back())) +
                                     " altspan=" + (graph_snp_result.graph_haplotype_path_alt_read_starts.empty() ?
                                         std::string("none") : std::to_string(graph_snp_result.graph_haplotype_path_alt_read_starts.back()) +
                                         "-" + std::to_string(graph_snp_result.graph_haplotype_path_alt_read_ends.back())) +
                                     " seq=" + (graph_snp_result.graph_haplotype_sequences.empty() ?
                                         std::string("none") : graph_snp_result.graph_haplotype_sequences.back()) +
                                     " paths_detail=" + [&]() {
                                         std::string detail;
                                         for (std::size_t path = 0; path < graph_snp_result.graph_haplotype_sequences.size(); ++path) {
                                             if (!detail.empty()) detail += ";";
                                             detail += std::to_string(path) + ":tid=" +
                                                 std::to_string(path < graph_snp_result.graph_haplotype_path_tids.size() ?
                                                     graph_snp_result.graph_haplotype_path_tids[path] : -99) +
                                                 ",start=" + std::to_string(path < graph_snp_result.graph_haplotype_path_starts.size() ?
                                                     graph_snp_result.graph_haplotype_path_starts[path] : -99) +
                                                 ",end=" + std::to_string(path < graph_snp_result.graph_haplotype_path_ends.size() ?
                                                     graph_snp_result.graph_haplotype_path_ends[path] : -99) +
                                                 ",seq=" + graph_snp_result.graph_haplotype_sequences[path];
                                         }
                                         return detail;
                                     }());
        // EventMap is a per-haplotype Host state, not merely the flattened
        // candidate union.  The same Q10 graph branch must retain the SNP in
        // at least one reference-connected haplotype EventMap.
        const auto graph_event_present = [](const auto& result, const std::string& ref,
                                            const std::string& alt, const std::int32_t position) {
            return std::any_of(result.graph_haplotype_event_maps.begin(),
                               result.graph_haplotype_event_maps.end(), [&](const auto& event_map) {
                return std::any_of(event_map.events.begin(), event_map.events.end(),
                    [&](const auto& event) {
                        return event.position == position && event.reference == ref &&
                            event.alternate == alt && event.graph_derived;
                    });
            });
        };
        if (!graph_event_present(graph_snp_result, "A", "C", 20) ||
            graph_snp_result.graph_haplotype_event_count == 0)
            throw std::runtime_error("per-haplotype graph EventMap lost SNP provenance");

        // A phased pair of graph mismatches is emitted as one same-length MNP
        // when max-mnp-distance is enabled.
        fastgatk::io::ReadBatch graph_mnp_reads;
        graph_mnp_reads.offsets.push_back(0);
        graph_mnp_reads.cigar_offsets.push_back(0);
        for (std::size_t record = 0; record < 8; ++record) {
            auto sequence = graph_snp_reference.substr(0, 41);
            sequence[20] = 'C';
            sequence[21] = 'T';
            graph_mnp_reads.bases.insert(graph_mnp_reads.bases.end(), sequence.begin(), sequence.end());
            graph_mnp_reads.qualities.insert(graph_mnp_reads.qualities.end(), sequence.size(), 30);
            if (record < 4) {
                graph_mnp_reads.qualities[graph_mnp_reads.qualities.size() - 21] = 10;
                graph_mnp_reads.qualities[graph_mnp_reads.qualities.size() - 20] = 10;
            }
            graph_mnp_reads.offsets.push_back(static_cast<std::uint32_t>(graph_mnp_reads.bases.size()));
            graph_mnp_reads.positions.push_back(0);
            graph_mnp_reads.tids.push_back(0);
            graph_mnp_reads.mapq.push_back(60);
            graph_mnp_reads.flags.push_back(0);
            graph_mnp_reads.cigar_ops.push_back(
                fastgatk::io::CigarOp{41, match}.pack());
            graph_mnp_reads.cigar_offsets.push_back(
                static_cast<std::uint32_t>(graph_mnp_reads.cigar_ops.size()));
        }
        auto graph_mnp_options = graph_snp_options;
        graph_mnp_options.max_mnp_distance = 1;
        const auto graph_mnp_result = fastgatk::calling::run(
            graph_mnp_reads, {graph_snp_reference}, graph_mnp_options);
        bool graph_mnp_found = false;
        bool graph_mnp_component_snp = false;
        for (const auto& candidate : graph_mnp_result.candidates) {
            graph_mnp_found = graph_mnp_found ||
                (candidate.position == 20 && candidate.reference_allele == "AC" &&
                 candidate.alternate_allele == "CT");
            graph_mnp_component_snp = graph_mnp_component_snp ||
                (candidate.position == 20 && candidate.reference_allele.empty() &&
                 candidate.alternate == 'C') ||
                (candidate.position == 21 && candidate.reference_allele.empty() &&
                 candidate.alternate == 'T');
        }
        if (!graph_mnp_found || graph_mnp_result.graph_mnp_candidates == 0 || graph_mnp_component_snp)
            throw std::runtime_error("graph-supported MNP mismatch was not materialized");
        if (!graph_event_present(graph_mnp_result, "AC", "CT", 20))
            throw std::runtime_error("per-haplotype graph EventMap did not retain MNP provenance");

        // EventMap.makeCompoundEvents combines events with the same start;
        // they must not escape as a SNP plus a second ALT.  The explicit
        // CIGAR avoids depending on a particular SW tie-break placement:
        // A>C followed by an insertion after that A is one A>CGG event.
        auto graph_compound_event_haplotype = graph_snp_reference;
        graph_compound_event_haplotype[20] = 'C';
        graph_compound_event_haplotype.insert(21, "GG");
        const auto graph_compound_event_map =
            fastgatk::calling::event_map_from_haplotype_cigar(
                graph_snp_reference, graph_compound_event_haplotype, 0, 0,
                "21M2I59M", 0);
        const auto compound_event_present = std::any_of(
            graph_compound_event_map.events.begin(), graph_compound_event_map.events.end(),
            [](const auto& event) {
                return event.position == 20 && event.reference == "A" &&
                event.alternate == "CGG" && event.graph_derived;
            });
        if (!compound_event_present || graph_compound_event_map.events.size() != 1)
            throw std::runtime_error("per-haplotype EventMap did not combine SNP/insertion event");
        // EventMap stores offsets in the complete reference coordinate frame,
        // not local positions in a padded assembly slice. Keep the public
        // Host helper honest at a nonzero CIGAR start; graph provenance uses
        // the same transition after its Kokkos-scored SW traceback.
        const std::string shifted_graph_snp_reference =
            std::string(100, 'N') + graph_snp_reference;
        const auto shifted_compound_event_map =
            fastgatk::calling::event_map_from_haplotype_cigar(
                shifted_graph_snp_reference, graph_compound_event_haplotype, 7, 100,
                "21M2I59M", 0);
        const auto shifted_compound_event_present = std::any_of(
            shifted_compound_event_map.events.begin(), shifted_compound_event_map.events.end(),
            [](const auto& event) {
                return event.tid == 7 && event.position == 120 && event.reference == "A" &&
                event.alternate == "CGG";
            });
        if (!shifted_compound_event_present || shifted_compound_event_map.events.size() != 1)
            throw std::runtime_error("per-haplotype EventMap lost complete-reference coordinate");
        // Softclip CIGARs span the complete haplotype.  The SW alignment's
        // read_start is the first aligned base, not the CIGAR walk cursor;
        // a leading soft clip must therefore be consumed exactly once.
        const auto leading_softclip_event_map =
            fastgatk::calling::event_map_from_haplotype_cigar(
                graph_snp_reference, std::string{"GG"} + graph_snp_reference, 0, 0,
                "2S80M", 0);
        if (!leading_softclip_event_map.events.empty())
            throw std::runtime_error("leading soft clip created a spurious EventMap event");
        auto graph_snp_deletion_haplotype = graph_snp_reference;
        graph_snp_deletion_haplotype[20] = 'C';
        graph_snp_deletion_haplotype.erase(21, 2);
        const auto graph_snp_deletion_event_map =
            fastgatk::calling::event_map_from_haplotype_cigar(
                graph_snp_reference, graph_snp_deletion_haplotype, 0, 0,
                "21M2D57M", 0);
        const auto snp_deletion_event_present = std::any_of(
            graph_snp_deletion_event_map.events.begin(),
            graph_snp_deletion_event_map.events.end(), [](const auto& event) {
                return event.position == 20 && event.reference == "ACT" &&
                event.alternate == "C" && event.graph_derived;
            });
        if (!snp_deletion_event_present || graph_snp_deletion_event_map.events.size() != 1)
            throw std::runtime_error("per-haplotype EventMap did not combine SNP/deletion event");
        auto graph_indel_compound_haplotype = graph_snp_reference.substr(0, 21);
        graph_indel_compound_haplotype += "TT";
        graph_indel_compound_haplotype += graph_snp_reference.substr(23);
        const auto graph_indel_compound_event_map =
            fastgatk::calling::event_map_from_haplotype_cigar(
                graph_snp_reference, graph_indel_compound_haplotype, 0, 0,
                "21M2I2D57M", 0);
        const auto indel_compound_event_present = std::any_of(
            graph_indel_compound_event_map.events.begin(),
            graph_indel_compound_event_map.events.end(), [](const auto& event) {
                return event.position == 20 && event.reference == "ACT" &&
                event.alternate == "ATT" && event.graph_derived;
            });
        if (!indel_compound_event_present || graph_indel_compound_event_map.events.size() != 1)
            throw std::runtime_error("per-haplotype EventMap did not combine insertion/deletion event");

        // AssemblyRegion boundaries are part of the PairHMM contract: two
        // distant active loci must not be combined into one contig-wide local
        // haplotype set. The synthetic pair is deliberately > max_region_size
        // apart while each locus has enough support to become active.
        fastgatk::io::ReadBatch separated_region_reads;
        separated_region_reads.offsets.push_back(0);
        separated_region_reads.cigar_offsets.push_back(0);
        const std::string separated_reference = unique_five_mer_reference(512);
        for (const auto [variant_position, alternate] :
             std::vector<std::pair<std::size_t, char>>{{42, 'T'}, {401, 'G'}}) {
            const auto start = variant_position - 10;
            for (std::size_t record = 0; record < 8; ++record) {
                auto sequence = separated_reference.substr(start, 21);
                if (record >= 4) sequence[10] = alternate;
                separated_region_reads.bases.insert(separated_region_reads.bases.end(), sequence.begin(), sequence.end());
                separated_region_reads.qualities.insert(separated_region_reads.qualities.end(), sequence.size(), 30);
                separated_region_reads.offsets.push_back(static_cast<std::uint32_t>(separated_region_reads.bases.size()));
                separated_region_reads.positions.push_back(static_cast<std::int32_t>(start));
                separated_region_reads.tids.push_back(0);
                separated_region_reads.mapq.push_back(60);
                separated_region_reads.flags.push_back(0);
                separated_region_reads.cigar_ops.push_back(
                    fastgatk::io::CigarOp{21, match}.pack());
                separated_region_reads.cigar_offsets.push_back(
                    static_cast<std::uint32_t>(separated_region_reads.cigar_ops.size()));
            }
        }
        auto separated_region_options = options;
        separated_region_options.graph_allow_non_unique_kmers_in_ref = false;
        const auto separated_result = fastgatk::calling::run(
            separated_region_reads, {separated_reference}, separated_region_options);
        if (separated_result.assembly_regions < 2 ||
            separated_result.pairhmm_assembly_region_groups != 2 ||
            separated_result.pairhmm_unassigned_candidates != 0 ||
            !separated_result.pairhmm_assembly_region_partitioned ||
            separated_result.assembly_cross_region_candidates != 0)
            throw std::runtime_error("PairHMM did not preserve separated AssemblyRegion boundaries: regions=" +
                                     std::to_string(separated_result.assembly_regions) +
                                     " groups=" + std::to_string(separated_result.pairhmm_assembly_region_groups) +
                                     " unassigned=" + std::to_string(separated_result.pairhmm_unassigned_candidates) +
                                     " assembly_cross=" + std::to_string(separated_result.assembly_cross_region_candidates));

        // AssemblyRegionWalker preserves independently popped regions even
        // when their padding overlaps.  The padding is only an assembly
        // context, not a request to coalesce the source callRegion()
        // transactions.  Keeping this explicit prevents a tempting but
        // incorrect whole-window PairHMM optimization.
        fastgatk::io::ReadBatch nearby_region_reads;
        nearby_region_reads.offsets.push_back(0);
        nearby_region_reads.cigar_offsets.push_back(0);
        const std::string nearby_reference = unique_five_mer_reference(256);
        for (const auto [variant_position, alternate] :
             std::vector<std::pair<std::size_t, char>>{{42, 'C'}, {142, 'G'}}) {
            const auto start = variant_position - 10;
            for (std::size_t record = 0; record < 8; ++record) {
                auto sequence = nearby_reference.substr(start, 21);
                if (record >= 4) sequence[10] = alternate;
                nearby_region_reads.bases.insert(
                    nearby_region_reads.bases.end(), sequence.begin(), sequence.end());
                nearby_region_reads.qualities.insert(nearby_region_reads.qualities.end(), sequence.size(), 30);
                nearby_region_reads.offsets.push_back(
                    static_cast<std::uint32_t>(nearby_region_reads.bases.size()));
                nearby_region_reads.positions.push_back(static_cast<std::int32_t>(start));
                nearby_region_reads.tids.push_back(0);
                nearby_region_reads.mapq.push_back(60);
                nearby_region_reads.flags.push_back(0);
                nearby_region_reads.cigar_ops.push_back(
                    fastgatk::io::CigarOp{21, match}.pack());
                nearby_region_reads.cigar_offsets.push_back(
                    static_cast<std::uint32_t>(nearby_region_reads.cigar_ops.size()));
            }
        }
        auto nearby_region_options = separated_region_options;
        nearby_region_options.min_assembly_region_size = 20;
        nearby_region_options.max_assembly_region_size = 50;
        nearby_region_options.max_probability_propagation_distance = 10;
        nearby_region_options.assembly_region_padding = 128;
        const auto nearby_result = fastgatk::calling::run(
            nearby_region_reads, {nearby_reference}, nearby_region_options);
        if (nearby_result.assembly_regions < 2 ||
            nearby_result.assembly_region_union_count != nearby_result.assembly_regions ||
            nearby_result.assembly_region_union_merges != 0 ||
            nearby_result.pairhmm_assembly_region_groups != nearby_result.assembly_regions)
            throw std::runtime_error("overlapping AssemblyRegion halos were incorrectly coalesced: raw=" +
                                     std::to_string(nearby_result.assembly_regions) +
                                     " union=" + std::to_string(nearby_result.assembly_region_union_count) +
                                     " merges=" + std::to_string(nearby_result.assembly_region_union_merges) +
                                     " groups=" + std::to_string(nearby_result.pairhmm_assembly_region_groups));

        // Complete raw flow tags select GATK FlowBased's dedicated direct
        // flow-space Kokkos engine (not the distinct FlowBasedHMM recurrence);
        // this guards the end-to-end Host decoder -> Kokkos model-selection
        // path while ordinary batches above continue to use regular mode.
        fastgatk::io::ReadBatch flow_reads;
        flow_reads.offsets.push_back(0);
        flow_reads.cigar_offsets.push_back(0);
        flow_reads.flow_tp_offsets.push_back(0);
        flow_reads.flow_t0_offsets.push_back(0);
        flow_reads.flow_order_offsets.push_back(0);
        const std::string flow_reference = multiallelic_reference;
        for (std::size_t record = 0; record < 4; ++record) {
            // Keep a real reference-connected graph branch; the direct Flow
            // likelihood is
            // selected after AssemblyRegion trimming, not from raw flow tags
            // on an otherwise unassemblable read.  The 41-base read also
            // exceeds HC's ten-base post-trimming minimum.
            auto sequence = flow_reference.substr(0, 41);
            sequence[20] = 'C';
            flow_reads.bases.insert(flow_reads.bases.end(), sequence.begin(), sequence.end());
            flow_reads.qualities.insert(flow_reads.qualities.end(), sequence.size(), 30);
            flow_reads.offsets.push_back(static_cast<std::uint32_t>(flow_reads.bases.size()));
            flow_reads.positions.push_back(0);
            flow_reads.tids.push_back(0);
            flow_reads.mapq.push_back(60);
            flow_reads.flags.push_back(0);
            flow_reads.cigar_ops.push_back(fastgatk::io::CigarOp{
                static_cast<std::uint32_t>(sequence.size()), match}.pack());
            flow_reads.cigar_offsets.push_back(static_cast<std::uint32_t>(flow_reads.cigar_ops.size()));
            flow_reads.flow_tp.insert(flow_reads.flow_tp.end(), sequence.size(), 0);
            flow_reads.flow_tp_offsets.push_back(static_cast<std::uint32_t>(flow_reads.flow_tp.size()));
            flow_reads.flow_t0_offsets.push_back(static_cast<std::uint32_t>(flow_reads.flow_t0_phred.size()));
            flow_reads.flow_orders.insert(flow_reads.flow_orders.end(), {'A', 'C', 'G', 'T'});
            flow_reads.flow_order_offsets.push_back(static_cast<std::uint32_t>(flow_reads.flow_orders.size()));
            flow_reads.flow_max_hmer.push_back(12);
        }
        auto flow_options = options;
        flow_options.min_alt_support = 2;
        const auto flow_result = fastgatk::calling::run(
            flow_reads, {flow_reference}, flow_options);
        if (!flow_result.pairhmm_used || flow_result.pairhmm_error_model != "flow" ||
            flow_result.pairhmm_flow_reads == 0 || flow_result.pairhmm_flow_haplotypes == 0 ||
            flow_result.pairhmm_flow_clipping_fallbacks != 0)
            throw std::runtime_error("complete flow tags did not select FlowBased direct likelihood: used=" +
                std::to_string(flow_result.pairhmm_used) + " model=" + flow_result.pairhmm_error_model +
                " reads=" + std::to_string(flow_result.pairhmm_flow_reads) +
                " haplotypes=" + std::to_string(flow_result.pairhmm_flow_haplotypes) +
                " clipped=" + std::to_string(flow_result.pairhmm_flow_haplotypes_clipped) +
                " fallbacks=" + std::to_string(flow_result.pairhmm_flow_clipping_fallbacks) +
                " pairs=" + std::to_string(flow_result.pairhmm_pairs));
        if (result.read_haplotype_cigar_filtered_pairs > result.read_haplotype_cigar_pairs ||
            result.read_haplotype_softclip_filtered_pairs >
                result.read_haplotype_cigar_filtered_pairs)
            throw std::runtime_error("read-haplotype CIGAR evidence gate counters are invalid");
        std::cout << "{\"status\":\"pass\",\"insertion_sites\":"
                  << result.candidate_sites << ",\"deletion_sites\":"
                  << deletion_result.candidate_sites << ",\"softclip_sites\":"
                  << softclip_result.candidate_sites << ",\"reverse_softclip_sites\":"
                  << reverse_softclip_result.candidate_sites << ",\"multiallelic_candidates\":"
                  << multiallelic_candidates << ",\"reverse_cigar_insertion_alt\":\""
                  << reverse_insertion_alt << "\",\"reverse_cigar_insertion_position\":"
                  << (reverse_insertion_found ? 5 : -1) << ",\"mixed_cigar_candidates\":"
                  << mixed_result.candidate_sites << ",\"graph_variant_candidates\":"
                  << internal_softclip_result.graph_variant_candidates
                  << ",\"graph_snp_candidates\":"
                  << graph_snp_result.graph_snp_candidates
                  << ",\"graph_snp_posterior_pairs\":"
                  << graph_snp_result.pairhmm_graph_snp_posterior_pairs
                  << ",\"graph_mnp_candidates\":"
                  << graph_mnp_result.graph_mnp_candidates
                  << ",\"compound_cigar_candidates\":"
                  << compound_result.compound_cigar_candidates
                  << ",\"compound_softclip_candidates\":"
                  << compound_softclip_result.compound_cigar_candidates
                  << ",\"read_haplotype_cigar_filtered_pairs\":"
                  << result.read_haplotype_cigar_filtered_pairs
                  << ",\"read_haplotype_softclip_filtered_pairs\":"
                  << result.read_haplotype_softclip_filtered_pairs
                  << ",\"graph_deletion_candidates\":"
                  << graph_deletion_result.graph_deletion_candidates << "}\n";
        Kokkos::finalize();
        return 0;
    } catch (...) {
        Kokkos::finalize();
        throw;
    }
}
