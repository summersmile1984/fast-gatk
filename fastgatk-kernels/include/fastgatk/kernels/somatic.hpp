#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace fastgatk::kernels {

// Per-candidate, per-read log10 likelihoods produced by the PairHMM stage.
// The matrices use candidate-major contiguous rows.  A value <= -1e299 is
// treated as missing so this API is independent of HTSlib sentinels.
struct SomaticLikelihoodResult {
    std::vector<double> tlod;
    // In Mutect2 reference-confidence mode this carries the TLOD of the
    // synthetic <NON_REF> allele at each concrete EventMap locus.  Its
    // likelihood row is constructed from the per-fragment next-best concrete
    // allele likelihoods, exactly as AlleleLikelihoods.addNonReferenceAllele
    // does before SomaticGenotypingEngine evaluates its all-alleles-vs-
    // without-allele evidence.  The matrix reduction itself remains the
    // regular Kokkos multiallelic kernel.
    std::vector<double> non_reference_tlod;
    // Posterior mean allele fraction of the same synthetic <NON_REF> row.
    // Concrete ALT fractions in reference-confidence mode are evaluated in
    // that shared REF/ALT/<NON_REF> model as well, rather than by the normal
    // biallelic convenience result.
    std::vector<double> non_reference_allele_fraction;
    // GATK SomaticGenotypingEngine::diploidAltLogOdds: log10 likelihood of
    // hom-ref minus the ref/alt heterozygous likelihood. Mutect2 writes this
    // as NLOD and applies --normal-lod in the same log10 domain.
    std::vector<double> normal_log10_odds;
    std::vector<double> best_allele_fraction;
    std::vector<double> reference_log10_likelihood;
    std::vector<double> best_log10_likelihood;
    std::vector<std::size_t> informative_reads;
    // Number of evidence units supplied to the Kokkos reduction. For
    // Mutect2 this is the post-Host fragment/read-name grouping cardinality;
    // for callers that do not group evidence it equals the input read count.
    std::size_t evidence_groups = 0;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

// Deterministic posterior decomposition for a tumor candidate.  The dense
// per-read likelihood arithmetic remains Kokkos-owned; Host code supplies
// prior/filter policy and maps the posterior fields to VCF annotations.
// Probabilities are ordinary [0,1] values and the three posterior vectors sum
// to one up to floating-point normalization error.
struct SomaticPosteriorResult {
    std::vector<double> somatic_probability;
    std::vector<double> germline_probability;
    std::vector<double> artifact_probability;
    std::vector<double> contamination_adjusted_allele_fraction;
    std::vector<double> orientation_bias_probability;
    std::vector<double> somatic_log10_evidence;
    std::vector<double> germline_log10_evidence;
    std::vector<double> artifact_log10_evidence;
    std::vector<std::size_t> informative_reads;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

// Evaluate a deterministic somatic mixture likelihood for each candidate.
// The null model is all-reference; the alternative uses the same two-state
// Dirichlet variational posterior/evidence form as GATK's
// SomaticLikelihoodsEngine.  minimum_allele_fraction=0 gives GATK's default
// [1,1] pseudocounts; a positive value uses the Beta(1+epsilon,1) prior that
// GATK derives from minAF.  The legacy allele_fraction_grid argument is
// retained for ABI/source compatibility and validation, but the evidence path
// is no longer a max-over-grid shortcut.
// A read with one missing allele likelihood still contributes through its
// finite allele, matching NaturalLogUtils' -Infinity handling.
// Host code owns prior/filter policy, while Kokkos owns the dense per-read
// mixture arithmetic and stable candidate ordering.
SomaticLikelihoodResult calculate_somatic_likelihood_kokkos(
    const std::vector<double>& reference_read_likelihoods,
    const std::vector<double>& alternate_read_likelihoods,
    std::size_t candidate_count,
    std::size_t read_count,
    std::size_t allele_fraction_grid = 101,
    double minimum_allele_fraction = 0.0);

// Sparse equivalent of calculate_somatic_likelihood_kokkos.  Host code owns
// the irregular candidate-to-fragment membership and supplies CSR row
// offsets; Kokkos still owns every per-fragment mixture update, evidence
// calculation, and NLOD reduction.  `candidate_offsets` has
// candidate_count + 1 entries and partitions the paired likelihood arrays.
// Cells absent from a candidate row are semantically equivalent to a dense
// row with both allele likelihoods missing, so they do not affect evidence.
// `evidence_group_count` is the cardinality of the original grouped evidence
// collection and is retained for diagnostics/metrics only.
SomaticLikelihoodResult calculate_somatic_likelihood_sparse_kokkos(
    const std::vector<double>& reference_group_likelihoods,
    const std::vector<double>& alternate_group_likelihoods,
    const std::vector<std::uint32_t>& candidate_offsets,
    std::size_t candidate_count,
    std::size_t evidence_group_count,
    std::size_t allele_fraction_grid = 101,
    double minimum_allele_fraction = 0.0);

// Evaluate GATK SomaticLikelihoodsEngine-style all-alleles-vs-without-alt
// log evidence for one locus.  `log10_likelihoods` is an allele-major,
// row-contiguous matrix with `allele_count` rows and `read_count` columns;
// row zero is REF and rows 1..N are concrete ALT alleles (an optional
// <NON_REF> row may be supplied as the final row by a caller using reference
// confidence).  The returned vectors have `allele_count - 1` entries, one
// for each non-REF row, and use the same row order as the input.  By default
// every evidence evaluation uses flat Dirichlet pseudocounts [1,...,1],
// matching GATK's minAF=0 somatic model; a positive minimum_allele_fraction
// changes each ALT pseudocount using GATK's Beta prior.  Missing/invalid values (<= -1e299, infinities
// or NaN) are ignored per allele; a read with at least one valid allele remains
// evidence.  Stable row/column order is retained across Kokkos execution spaces.
// The optional minimum_allele_fraction applies the same ALT pseudocount prior
// to every non-REF row as GATK's SomaticGenotypingEngine.
SomaticLikelihoodResult calculate_somatic_multiallelic_likelihood_kokkos(
    const std::vector<double>& log10_likelihoods,
    std::size_t allele_count,
    std::size_t read_count,
    double minimum_allele_fraction = 0.0);

// Evaluate a compact somatic posterior model over candidate-major tumor and
// normal likelihood matrices.  The tumor matrix has candidate_count rows and
// `read_count` columns; a non-empty normal matrix may have a different number
// of columns (its width is inferred as normal_matrix.size()/candidate_count),
// matching GATK's independent tumor/normal AlleleLikelihoods evidence counts.
// Empty normal matrices mean that no normal sample is available.  `f1r2`/`r1f2`
// are alternate-support counts used by the orientation artifact term.
// Contamination is a fraction in [0,1); priors must be positive and are
// normalized in log10 space on device.
SomaticPosteriorResult calculate_somatic_posterior_kokkos(
    const std::vector<double>& tumor_reference_read_likelihoods,
    const std::vector<double>& tumor_alternate_read_likelihoods,
    const std::vector<double>& normal_reference_read_likelihoods,
    const std::vector<double>& normal_alternate_read_likelihoods,
    const std::vector<std::uint32_t>& f1r2,
    const std::vector<std::uint32_t>& r1f2,
    std::size_t candidate_count,
    std::size_t read_count,
    double contamination = 0.0,
    double somatic_prior = 1.0e-4,
    double germline_prior = 1.0e-3,
    double artifact_prior = 1.0e-2,
    std::size_t allele_fraction_grid = 101);

}  // namespace fastgatk::kernels
