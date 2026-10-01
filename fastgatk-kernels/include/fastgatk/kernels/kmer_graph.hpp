#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace fastgatk::kernels {

struct KmerGraphInput {
    std::vector<std::uint8_t> bases;
    std::vector<std::uint32_t> offsets;
    // Optional read contig IDs.  If omitted, all reads use contig 0.
    std::vector<std::int32_t> tids;
    // Optional reference contigs flattened with the same offset convention as
    // `bases`.  Reference k-mers anchor graph connectivity and are never
    // treated as read support when reporting `input_kmers`.
    std::vector<std::uint8_t> reference_bases;
    std::vector<std::uint32_t> reference_offsets;
    // Optional coordinate metadata for each reference segment.  A segment is
    // still flattened in `reference_bases`, but this lets Host consumers map
    // materialized paths back to the originating BAM contig/window.
    std::vector<std::int32_t> reference_tids;
    std::vector<std::int32_t> reference_starts;
    // Optional half-open reference window ends, paired with reference_starts.
    // Keeping the end avoids treating a truncated/dead-end graph path as a
    // large deletion when it is aligned back to the reference scaffold.
    std::vector<std::int32_t> reference_ends;
    // Optional BAM/SAM flags for reads.  K-mer threading preserves the stored
    // BAM base order, matching GATKRead.getBases(); flags remain available to
    // callers for alignment-aware policy outside the graph payload.
    // This field is appended to preserve source compatibility with existing
    // aggregate initializers that stop at reference_ends.
    std::vector<std::uint16_t> flags;
    // Optional sample ordinal for each read.  Missing values are treated as
    // sample 0, matching the single-sample HC/Mutect2 default.  The field is
    // appended so older aggregate initializers remain source-compatible.
    std::vector<std::int32_t> sample_ids;
};

struct KmerGraphOptions {
    std::uint32_t k = 21;
    std::uint32_t min_count = 1;
    // Minimum number of read observations required to retain a non-reference
    // branch. This is the native equivalent of GATK's --min-pruning; the
    // reference scaffold remains available independently of this threshold.
    std::uint32_t min_pruning = 1;
    // Path-count protection is explicit so the same graph policy can be
    // tuned for desktop, node-local and accelerator memory budgets.  GATK's
    // K-best search terminates at reference sink vertices, not after a fixed
    // number of bases; max_depth == 0 therefore preserves that source-to-sink
    // default.  A positive value bounds only the unrejoined non-reference
    // branch.  It must never cap the shared reference prefix/suffix, which
    // would turn every AssemblyRegion longer than the guard into zero
    // haplotypes rather than a bounded alternate traversal.
    std::size_t max_paths = 64;
    std::size_t max_depth = 0;
    // GATK ReadThreadingAssembler-shaped dangling-branch policy.  A branch
    // that leaves the reference scaffold and terminates without rejoining it
    // must span at least this many materialized bases to be retained.  Zero
    // disables the length floor (useful for low-coverage compatibility
    // probes).  `recover_all_dangling_branches` mirrors GATK's *connected*
    // graph walk: it permits forks while choosing the heaviest graph edge.
    // It does not rescue a disconnected component by aligning it globally to
    // the reference; ReadThreadingAssembler removes those components after
    // dangling-end recovery.
    std::uint32_t min_dangling_branch_length = 4;
    bool recover_all_dangling_branches = false;
    // When an alternate cycle is detected, GATK retries the graph with a
    // larger k-mer (normally +10) unless this switch is set.  The retry is
    // Host orchestration; each individual graph/count stage still uses the
    // same Kokkos API.
    bool dont_increase_kmer_sizes_for_cycles = false;
    // Minimum exact read/reference matches required before a connected
    // dangling-end CIGAR can rejoin the graph. This mirrors GATK's
    // minMatchingBasesToDanglingEndRecovery. The GATK 4 default is -1,
    // which selects its legacy mismatch-bounded dangling-head merge rule;
    // non-negative values select the exact-prefix rule.
    // Appended to preserve existing aggregate initializer compatibility.
    std::int32_t min_dangling_matching_bases = -1;
    // GATK skips a k-mer attempt when the reference haplotype contains a
    // repeated k-mer unless this compatibility switch is enabled.  The
    // low-level API remains permissive by default for source compatibility;
    // HC/Mutect2 set the option explicitly to GATK's default (false).
    bool allow_non_unique_kmers_in_ref = true;
    // Number of samples that must independently reach min_pruning for a
    // branch to survive.  This is GATK MultiSampleEdge's pruning semantics;
    // a missing sample map is interpreted as one sample.
    std::uint32_t num_pruning_samples = 1;
    // Optional adaptive chain-pruning policy.  The defaults preserve the
    // historical fixed min_pruning path; callers can enable the GATK-shaped
    // probabilistic gate explicitly when they provide multi-sample evidence.
    bool use_adaptive_pruning = false;
    double initial_error_rate_for_pruning = 0.001;
    // ReadThreadingAssemblerArgumentCollection stores these in natural-log
    // units: MathUtils.log10ToLog(1) and MathUtils.log10ToLog(4).  They are
    // intentionally not the convenient-looking 2/3 values: the seeding
    // gate at ln(10^4) determines whether a one-read chain is removed before
    // dangling-end recovery.
    double pruning_log_odds_threshold = 2.30258509299404568402;
    double pruning_seeding_log_odds_threshold = 9.21034037197618273607;
    std::uint32_t max_unpruned_variants = 100;
    // Experimental GATK linked-de-Bruijn mode deliberately keeps the
    // sequence graph unsimplified; the raw bounded path topology remains
    // available for Host recovery and telemetry.
    bool disable_seqgraph_simplification = false;
    // Linked-de-Bruijn mode retains junction-tree-like edge connectivity.
    // When artificial recovery is enabled, paths for high-support graph
    // edges that were not visited by the bounded traversal are synthesized
    // from the best completed path, matching GATK's uncovered-edge recovery.
    bool linked_de_bruijn_graph = false;
    bool disable_artificial_haplotype_recovery = false;
    // Optional ordered --kmer-size requests from the GATK CLI.  This is
    // appended so all existing aggregate initializers retain their meaning;
    // an empty list means use the scalar `k` field and preserves the low-level
    // API's historical behavior.
    std::vector<std::uint32_t> requested_kmer_sizes;
};

// Host-side SeqGraph rewrite input.  A path is represented as one contiguous
// sequence vertex; the rewriter then materializes shared prefixes/suffixes and
// the same diamond/tail topology that GATK's SeqGraph simplifier operates on.
// Keeping this contract independent of k-mer node IDs lets the exact rewrite
// rules be reused by CPU, OpenMP and accelerator graph front-ends.
struct SeqGraphPath {
    std::string sequence;
    std::int32_t tid = -1;
    std::int32_t start = -1;
    std::int32_t end = -1;
    bool reference = false;
};

struct SeqGraphSimplificationResult {
    std::size_t initial_nodes = 0;
    std::size_t initial_edges = 0;
    std::size_t final_nodes = 0;
    std::size_t final_edges = 0;
    std::size_t linear_chain_merges = 0;
    std::size_t diamond_merges = 0;
    std::size_t tail_merges = 0;
    std::size_t suffix_splits = 0;
    std::size_t suffix_merges = 0;
    // The transformed graph's path vertices.  These are deduplicated and
    // sorted, so callers can use them as a deterministic haplotype ordering
    // without exposing mutable graph node IDs.
    std::vector<std::string> path_sequences;
};

// Weighted graph edge used by the GATK ChainPruner-compatible primitive.
// `sample_multiplicities` contains one entry per observed sample (zero-count
// samples may be omitted); the pruning multiplicity is the Nth largest entry,
// where N is `KmerGraphOptions::num_pruning_samples`.
struct SeqGraphEdge {
    std::uint32_t from = 0;
    std::uint32_t to = 0;
    std::uint32_t multiplicity = 0;
    bool reference = false;
    std::vector<std::uint32_t> sample_multiplicities;
};

struct SeqGraphPruningResult {
    std::vector<std::uint8_t> keep_edges;
    std::size_t chain_count = 0;
    std::size_t pruned_chain_count = 0;
    std::size_t adaptive_pruned_chain_count = 0;
    double estimated_error_rate = 0.0;
};

// GATK's AdaptiveChainPruner uses Mutect2Engine.logLikelihoodRatio for the
// constant-error edge model. Exposing this scalar helper keeps the numerical
// contract testable without coupling callers to graph traversal details.
// The implementation follows Commons Math/GATK's digamma and binomial-log
// branches so Host, Serial and accelerator front-ends share one definition.
double seqgraph_constant_error_log_likelihood_ratio(
    std::uint32_t n_ref,
    std::uint32_t n_alt,
    double error_probability);

struct KmerGraphResult {
    std::size_t input_reads = 0;
    std::size_t input_kmers = 0;
    std::size_t nodes = 0;
    std::size_t edges = 0;
    std::size_t branching_nodes = 0;
    std::size_t reference_nodes = 0;
    std::size_t reference_edges = 0;
    std::size_t reference_connected_nodes = 0;
    std::size_t dangling_nodes = 0;
    std::size_t pruned_nodes = 0;
    std::size_t dangling_branch_paths = 0;
    std::size_t dangling_branch_bases = 0;
    std::size_t dangling_recovered_paths = 0;
    std::size_t dangling_recovered_bases = 0;
    // Number/bases of linked-de-Bruijn artificial haplotypes synthesized for
    // pivotal edges that were not covered by the bounded path traversal.
    std::size_t artificial_haplotype_recovery_paths = 0;
    std::size_t artificial_haplotype_recovery_bases = 0;
    // SeqGraph-style non-branching path compression performed after k-mer
    // pruning. `nodes/edges` remain the raw de Bruijn topology; these fields
    // expose the deterministic maximal-path topology consumed by traversal
    // ordering and path materialization telemetry.
    std::size_t seqgraph_nodes = 0;
    std::size_t seqgraph_edges = 0;
    std::size_t seqgraph_linear_chain_merges = 0;
    std::size_t seqgraph_diamond_merges = 0;
    std::size_t seqgraph_tail_merges = 0;
    std::size_t seqgraph_suffix_splits = 0;
    std::size_t seqgraph_suffix_merges = 0;
    std::size_t adaptive_pruned_nodes = 0;
    std::uint32_t kmer_size = 0;
    std::size_t kmer_iterations = 1;
    // k-mer sizes whose graph attempt GATK's createGraph accepted (its
    // "Using kmer size of N in read threading assembler" log).
    std::vector<std::uint32_t> kmer_sizes_used;
    bool has_non_reference_cycles = false;
    // K-mers observed more than once within any individual read/reference
    // sequence, matching ReadThreadingGraph.getNonUniqueKmers() semantics.
    // Sorted decoded strings keep the dynamic-k selection decision stable.
    std::vector<std::string> non_unique_kmers;
    // Number of repeated k-mers found in the reference scaffold only.  This
    // is separate from the combined read/reference `non_unique_kmers` set so
    // Host orchestration can implement GATK's skip-and-increase-k policy.
    std::size_t reference_non_unique_kmers = 0;
    bool reference_kmer_rejected = false;
    std::size_t reference_path_count = 0;
    std::size_t haplotype_path_count = 0;
    // Bounded, deterministic sequence materialization for the first paths.
    // Each sequence is a Host-owned k-mer path (one base appended per edge),
    // suitable for SW/PairHMM consumers without exposing graph node IDs.
    std::vector<std::string> haplotype_path_sequences;
    std::vector<std::int32_t> haplotype_path_tids;
    std::vector<std::int32_t> haplotype_path_starts;
    std::vector<std::int32_t> haplotype_path_ends;
    // Maximum read-edge support observed along each materialized path.  The
    // reference scaffold contributes no support, so callers can distinguish
    // a pure reference path from a read-supported alternate branch.
    std::vector<std::uint32_t> haplotype_path_support;
    // GraphBasedKBestHaplotypeFinder's accumulated log10 edge-probability
    // score for each materialized path.  This is the score carried by GATK's
    // Haplotype object and is deliberately distinct from read×haplotype
    // PairHMM likelihoods and from the later haplotype-to-reference SW score.
    // Host uses it when HaplotypeCaller must reduce an oversized allele map
    // for --max-genotype-count before marginalizing read likelihoods.
    std::vector<double> haplotype_path_scores;
    // True when traversal used at least one edge absent from the reference
    // scaffold.  Reference-only cycle expansions are retained for telemetry
    // but must not be interpreted as alternate haplotypes by callers.
    std::vector<std::uint8_t> haplotype_path_has_non_reference_edge;
    // Read-derived branch span (half-open, in the materialized sequence).
    // When a non-reference dangling-tail recovery edge enters a compressed
    // reference vertex, the span ends at that joining base rather than
    // absorbing the vertex's full downstream reference suffix.  This lets
    // Host projection ignore artificial leading/trailing gaps introduced by
    // bounded path halos.
    std::vector<std::uint32_t> haplotype_path_alt_read_starts;
    std::vector<std::uint32_t> haplotype_path_alt_read_ends;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
    bool used = false;
};

// Build a deterministic compact de Bruijn graph. Host owns k-mer encoding,
// sorting and stable node IDs; Kokkos performs occurrence/degree counting so
// the same graph stage runs on Serial/OpenMP/CUDA/HIP execution spaces.
KmerGraphResult build_kmer_graph_kokkos(
    const KmerGraphInput& input, KmerGraphOptions options = {});

// Apply the deterministic Host SeqGraph simplifier to materialized paths.
// This is intentionally a pure API: no Kokkos state is required, and the
// resulting topology is suitable for callers that want to run their own
// downstream path/likelihood kernel over the simplified graph.
SeqGraphSimplificationResult simplify_seqgraph_paths(
    const std::vector<SeqGraphPath>& paths,
    std::size_t min_common_tail_bases = 10);

// Apply GATK ChainPruner/AdaptiveChainPruner semantics to a weighted directed
// graph.  The function is Host-side because chain discovery is variable
// topology; all multiplicity vectors and the resulting keep mask are plain
// data that can be consumed by a Kokkos traversal kernel.
SeqGraphPruningResult prune_seqgraph_chains(
    std::size_t node_count,
    const std::vector<SeqGraphEdge>& edges,
    const KmerGraphOptions& options,
    const std::vector<std::string>* node_sequences = nullptr);

}  // namespace fastgatk::kernels
